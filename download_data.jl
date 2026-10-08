# Download the ERA5, JRA55, GLORYS, and ETOPO files for a forcing window and cache the regridded bathymetry
# for both resolutions, without building a model. Runs on the login node (CPU only), so GPU jobs find
# everything cached and the ranks of a distributed run never regrid concurrently.
#
#   julia --project download_data.jl [days=30]

using CDSAPI, CopernicusMarine
using Downloads: download
include("setup.jl")

end_date = start_date + Day(parse(Int, get(ARGS, 1, "30")) + 1)

for cells_per_degree in (6, 12)
    global_bottom_height(360cells_per_degree, 170cells_per_degree; halo = (7, 7, 4))
end

download(MetadataSet(:river_freshwater_flux, :iceberg_freshwater_flux; dataset = MultiYearJRA55(), start_date, end_date))

download(MetadataSet(:temperature, :salinity, :sea_ice_thickness, :sea_ice_concentration;
                     dataset = GLORYSDaily(), date = start_date))

download(MetadataSet(:eastward_velocity, :northward_velocity, :temperature, :dewpoint_temperature,
                     :surface_pressure, :total_precipitation,
                     :downwelling_shortwave_radiation, :downwelling_longwave_radiation;
                     dataset = ERA5HourlySingleLevel(), start_date, end_date))
