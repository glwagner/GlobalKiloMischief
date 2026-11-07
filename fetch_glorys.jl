# Download the GLORYS initial state of a run (temperature, salinity, sea ice thickness and concentration on START_DATE).
#   START_DATE=2021-06-30 julia --project=. fetch_glorys.jl

using NumericalEarth
using CopernicusMarine  # loads the GLORYS download extension
using Dates
using Downloads: download

start_date = DateTime(get(ENV, "START_DATE", "2015-01-01"))
download(MetadataSet(:temperature, :salinity, :sea_ice_thickness, :sea_ice_concentration;
                     dataset = GLORYSDaily(), date = start_date))
