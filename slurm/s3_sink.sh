#!/bin/bash
# Copy the output and checkpoints of the runs under <root> to S3 while they run, and free local space. Run on a
# login node (never in a job), for example in tmux. Every 30 minutes, for each run directory <root>/<run>:
#   - every file not modified in the last 30 minutes is uploaded to $PREFIX/<run>/ if it changed since its last upload,
#     with its md5 in the object's metadata, and verified (size and md5) before it counts as uploaded;
#   - local copies are removed only of verified files that have not changed since: closed output parts (a later
#     _part<n> of the same file exists) and checkpoints older than the newest two iterations.
# If the AWS credentials fail (the SSO token expires after 8-12 hours) nothing is uploaded or removed and the log says
# so on every pass until `~/.local/bin/aws sso login --profile aeolus --use-device-code` is run again.
#
#   bash slurm/s3_sink.sh <root>          # stop: touch <root>/.s3_sink_stop
root=${1:?usage: s3_sink.sh <root>}
aws=~/.local/bin/aws
profile=aeolus
prefix=s3://aeolus-dev/users/greg@aeolus.earth/delta/oceananigans-paper
bucket=aeolus-dev
key_prefix=users/greg@aeolus.earth/delta/oceananigans-paper
state=$root/.s3_sink_state  # verified uploads: <path> <size> <mtime> <md5>
log=$root/s3_sink.log
quiet_minutes=30
touch "$state"

say() { echo "$(date '+%F %T') $*" | tee -a "$log"; }

authenticated() {
    $aws sts get-caller-identity --profile $profile > /dev/null 2>&1 && return 0
    say "!!!!!!!! AWS AUTHENTICATION FAILED: nothing uploaded or removed. Run: ~/.local/bin/aws sso login --profile aeolus --use-device-code"
    return 1
}

verified() {  # path: uploaded, verified and unchanged since
    local f=$1 last
    last=$(grep -F "$f " "$state" | tail -1 | awk '{print $2, $3}')
    [ -n "$last" ] && [ "$last" = "$(stat -c '%s %Y' "$f")" ]
}

upload() {  # run, path
    local run=$1 f=$2 name size mtime md5 remote
    name=$(basename "$f")
    size=$(stat -c %s "$f"); mtime=$(stat -c %Y "$f")
    grep -q "^$f $size $mtime " "$state" && return 0
    md5=$(md5sum "$f" | cut -d' ' -f1)
    if ! $aws s3 cp "$f" "$prefix/$run/$name" --profile $profile --metadata md5=$md5 --only-show-errors; then
        say "upload FAILED: $run/$name"; return 1
    fi
    remote=$($aws s3api head-object --bucket $bucket --key "$key_prefix/$run/$name" --profile $profile \
             --query '[ContentLength, Metadata.md5]' --output text 2>/dev/null)
    if [ "$remote" = "$size	$md5" ] && [ "$(stat -c '%s %Y' "$f")" = "$size $mtime" ]; then
        echo "$f $size $mtime $md5" >> "$state"
        say "uploaded $run/$name ($((size / 1000000)) MB)"
    else
        say "VERIFY MISMATCH for $run/$name (remote: $remote, local: $size $md5); kept"
    fi
}

removable() {  # path: a closed output part or an old checkpoint
    local f=$1 dir name stem part iteration newest
    dir=$(dirname "$f"); name=$(basename "$f")
    if [[ $name =~ ^(.*)_part([0-9]+)\.jld2$ ]]; then
        stem=${BASH_REMATCH[1]}; part=${BASH_REMATCH[2]}
        ls "$dir" | grep -qE "^${stem}_part$((part + 1))\.jld2$"
    elif [[ $name =~ ^checkpoint.*_iteration([0-9]+)\.jld2$ ]]; then
        iteration=${BASH_REMATCH[1]}
        newest=$(ls "$dir" | sed -nE 's/^checkpoint.*_iteration([0-9]+)\.jld2$/\1/p' | sort -n | uniq | tail -2 | head -1)
        [ -n "$newest" ] && [ "$iteration" -lt "$newest" ]
    else
        return 1
    fi
}

say "s3 sink started for $root -> $prefix"
until [ -f "$root/.s3_sink_stop" ]; do
    if authenticated; then
        for dir in "$root"/*/; do
            run=$(basename "$dir")
            for f in $(find "$dir" -maxdepth 1 -type f ! -name '*.partial' -mmin +$quiet_minutes | sort); do
                upload "$run" "$f"
            done
            for f in $(find "$dir" -maxdepth 1 -type f -name '*.jld2' | sort); do
                if removable "$f" && verified "$f" && authenticated; then
                    rm -f -- "$f" && say "removed local $run/$(basename "$f") (verified in S3)"
                fi
            done
        done
    fi
    sleep 1800
done
say "s3 sink stopped ($root/.s3_sink_stop)"
