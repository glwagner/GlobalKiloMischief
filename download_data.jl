# Download the ERA5, JRA55, GLORYS, and ETOPO files for a forcing window without building a model.
# Runs on the login node (CPU only), so GPU jobs find their data cached.
#
#   julia --project download_data.jl [days=30]

using NumericalEarth
using CDSAPI, CopernicusMarine
using Dates

start_date = DateTime(2015, 1, 1)
end_date = start_date + Day(parse(Int, get(ARGS, 1, "30")) + 1)

download(Metadatum(:bottom_height; dataset = ETOPO2022()))

download(MetadataSet(:temperature, :salinity, :sea_ice_thickness, :sea_ice_concentration;
                     dataset = GLORYSDaily(), date = start_date))

download(MetadataSet(:river_freshwater_flux, :iceberg_freshwater_flux; dataset = MultiYearJRA55(), start_date, end_date))

download(MetadataSet(:eastward_velocity, :northward_velocity, :temperature, :dewpoint_temperature,
                     :surface_pressure, :total_precipitation,
                     :downwelling_shortwave_radiation, :downwelling_longwave_radiation;
                     dataset = ERA5HourlySingleLevel(), start_date, end_date))
