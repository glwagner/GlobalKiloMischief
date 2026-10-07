# Global 1/6ᵒ ocean--sea ice simulation on a single GPU. Each run picks up from the latest checkpoint,
# so run again with a larger `stop_days` to extend it.
#
#   julia --project sixth_degree.jl [stop_days=30]

include("setup.jl")

stop_days = parse(Int, get(ARGS, 1, "30"))
max_Δt = 10minutes

grid = global_grid(GPU(); cells_per_degree = 6)
model = ocean_sea_ice_model(grid; max_Δt, end_date = start_date + Day(stop_days + 1))
simulation = global_simulation(model; name = "sixth_degree", Δt = 1minute, max_Δt, stop_time = stop_days * days)

run!(simulation; pickup = true)
