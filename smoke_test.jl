# Smoke test: build the full 1/6ᵒ configuration on one GPU, downloading and caching every dataset it needs
# (ETOPO, GLORYS, ERA5, JRA55), then take a few hundred time steps.
#
#   julia --project smoke_test.jl [iterations=300]

include("setup.jl")

iterations = parse(Int, get(ARGS, 1, "300"))
max_Δt = 20minutes

grid = global_grid(GPU(); cells_per_degree = 6)
@info "Built grid" grid

model = ocean_sea_ice_model(grid; max_Δt, end_date = start_date + Day(3))
simulation = global_simulation(model; name = "smoke_test", Δt = 1minute, max_Δt, stop_time = 2days)
simulation.stop_iteration = iterations

run!(simulation)
