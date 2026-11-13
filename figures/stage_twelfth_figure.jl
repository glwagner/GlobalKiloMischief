# Everything global_simulation_globes.jl reads for the twelfth-degree figure, stitched across ranks into one Float32
# file, so the figure can be remade without the run output.
#   julia --project=<figure env> figures/stage_twelfth_figure.jl <run directory> <output.jld2>

include(joinpath(@__DIR__, "global_simulation_globes.jl"))

directory, output_path = ARGS[1], ARGS[2]
ocean_iteration, time = last(snapshots(directory, "ocean_surface"))
ice_iteration = Dict(t => i for (i, t) in snapshots(directory, "sea_ice"))[time]

uo = stitch(directory, "ocean_surface", "u", ocean_iteration; face_x = true)
vo = stitch(directory, "ocean_surface", "v", ocean_iteration; face_y = true)
ui = stitch(directory, "sea_ice", "u", ice_iteration; face_x = true)
vi = stitch(directory, "sea_ice", "v", ice_iteration; face_y = true)
ℵi = stitch(directory, "sea_ice", "ℵ", ice_iteration)
Nx, Ny = size(ℵi)
grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), z = (0, 1), halo = (7, 7, 1))
λ = Array(λnodes(grid, Center(), Center(), Center()))
φ = Array(φnodes(grid, Center(), Center(), Center()))
land = (uo[1:end-1, :] .== 0) .& (uo[2:end, :] .== 0) .& (vo[:, 1:end-1] .== 0) .& (vo[:, 2:end] .== 0)

f32(a) = Float32.(a)
jldsave(output_path; u_ocean = f32(uo), v_ocean = f32(vo), u_ice = f32(ui), v_ice = f32(vi), ice_concentration = f32(ℵi),
        longitude = f32(λ), latitude = f32(φ), land, time_seconds = time, snapshot = timestamp(time),
        notes = "u at (Face, Center), Nx+1 × Ny; v at (Center, Face), Nx × Ny+1; others at cell centers. " *
                "TripolarGrid 1/12ᵒ, 80ᵒS to the north poles; `land` as in speeds(); start 2015-01-01.")
@info "Wrote $output_path: $(timestamp(time)); sizes u $(size(uo)), v $(size(vo)), centers $(size(ℵi))"
