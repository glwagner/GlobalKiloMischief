# Print the array shapes of the fields `set_from_checkpoint!` reads, without reading the arrays.
#   julia --project=. regrid/checkpoint_shapes.jl <checkpoint.jld2>
using JLD2

names = ["ocean/model/tracers/T", "ocean/model/velocities/u", "ocean/model/velocities/v",
         "ocean/model/free_surface/displacement", "sea_ice/model/ice_thickness", "sea_ice/model/velocities/u"]

jldopen(ARGS[1], "r") do file
    for name in names
        dataset = JLD2.get_dataset(file, "simulation/model/" * name * "/data")
        println(name, "  ", Int.(reverse(dataset.dataspace.dimensions)))
    end
end
