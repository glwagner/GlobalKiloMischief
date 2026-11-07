# Fetch the ERA5 hourly forcing of a run from NCAR's public S3 mirror (no CDS account or queue), one month at a time,
# printing how long each month takes. Light enough for a login node: one time chunk of one variable at a time.
#
#   START_DATE=2021-06-30 julia --project=. fetch_era5_ncar.jl [days=91]

using NumericalEarth
using Dates
using Downloads: download

start_date = DateTime(get(ENV, "START_DATE", "2015-01-01"))
end_date = start_date + Day(parse(Int, get(ARGS, 1, "91")))
dataset = ERA5HourlySingleLevel(archive = NCARArchive())
names = (:eastward_velocity, :northward_velocity, :temperature, :dewpoint_temperature, :mean_sea_level_pressure,
         :total_precipitation, :downwelling_shortwave_radiation, :downwelling_longwave_radiation)

month_start = start_date
while month_start < end_date
    month_end = min(end_date, DateTime(Dates.year(month_start), Dates.month(month_start), Dates.daysinmonth(month_start), 23))
    elapsed = @elapsed for name in names
        download(Metadata(name; dataset, dates = (month_start, month_end)))
    end
    @info "ERA5 from NCAR, $(Date(month_start)) to $(Date(month_end)): $(round(elapsed / 60, digits = 1)) minutes"
    global month_start = month_end + Hour(1)
end
