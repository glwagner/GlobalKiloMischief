# Stage the inputs a run segment reads onto the NVMe allocation and point NumericalEarth there (source this file
# after env.sh, on the compute node). Cold reads from /work/hdd made the model build take up to 27 minutes.
#
#   source slurm/stage_inputs.sh <stop_days>
#
# The mirror holds real copies of the GLORYS initial condition, the 2015 JRA55 runoff, the bathymetry cache, and
# the ERA5 hours of the segment (from 12 days before `stop_days` to 2 days after, plus the first three days of
# 2015, which every model build reads). Every other file is a symlink to /work/hdd, so a forcing series that checks
# that all its files exist finds them. Downloads still go to /work/hdd.

hdd_data=/work/hdd/bhcr/glwagner/numericalearth_data
mirror=/work/nvme/bhcr/glwagner/numericalearth_mirror
stage_start=$(date +%s)
mkdir -p $mirror/ERA5 $mirror/GLORYS $mirror/JRA55 $mirror/field_cache

# Datasets the runs do not read stay on /work/hdd.
for dataset in $(ls $hdd_data); do
    case $dataset in ERA5|GLORYS|JRA55|field_cache) ;; *) [ -e $mirror/$dataset ] || ln -s $hdd_data/$dataset $mirror/$dataset ;; esac
done

# Small or read-at-every-build inputs are copied once (the repeat-year JRA55 files are not used).
for dataset in GLORYS field_cache; do
    for file in $(ls $hdd_data/$dataset); do
        [ -f $mirror/$dataset/$file ] && [ ! -L $mirror/$dataset/$file ] || cp $hdd_data/$dataset/$file $mirror/$dataset/$file
    done
done
for file in $(ls $hdd_data/JRA55); do
    case $file in
        RYF.*) [ -e $mirror/JRA55/$file ] || ln -s $hdd_data/JRA55/$file $mirror/JRA55/$file ;;
        *) [ -f $mirror/JRA55/$file ] && [ ! -L $mirror/JRA55/$file ] || cp $hdd_data/JRA55/$file $mirror/JRA55/$file ;;
    esac
done

# ERA5: the segment's hours as copies, every other hour as a symlink.
first_date=$(date -d "2015-01-01 + $(( $1 - 12 )) days" +%Y-%m-%d)
last_date=$(date -d "2015-01-01 + $(( $1 + 2 )) days" +%Y-%m-%d)
in_window() {  # ERA5 file name
    local d=$(echo $1 | sed 's/.*_\(2015-..-..\)T.*/\1/')
    [[ "$d" < "2015-01-04" ]] || { [[ ! "$d" < "$first_date" ]] && [[ ! "$d" > "$last_date" ]]; }
}
to_copy=$(mktemp -p $mirror)
for file in $(ls $hdd_data/ERA5 | grep "HourlySingleLevel_2015-..-..T..\.nc$"); do
    if in_window $file; then
        [ -L $mirror/ERA5/$file ] || [ ! -e $mirror/ERA5/$file ] && echo $file >> $to_copy
    elif [ ! -L $mirror/ERA5/$file ]; then
        ln -sfn $hdd_data/ERA5/$file $mirror/ERA5/$file  # creates the link, or turns an old copy back into one
    fi
done
xargs -a $to_copy -P 16 -I{} cp --remove-destination $hdd_data/ERA5/{} $mirror/ERA5/{}
echo "staged $(wc -l < $to_copy) ERA5 files ($first_date to $last_date) in $(( $(date +%s) - stage_start )) s; mirror holds $(du -sh --exclude='*.nc.partial' $mirror 2>/dev/null | cut -f1)"
mv $to_copy $mirror/.last_staged

export NUMERICALEARTH_DATA_DIRECTORY=$mirror
