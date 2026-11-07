using JLD2
using Oceananigans.Architectures: on_architecture

"""
    set_from_checkpoint!(model, directory; iteration, Nz_old = 100, depth = 6000, rank = 0)

Overwrite the ocean and sea ice state of the coupled `model` with checkpoint `iteration` in `directory`, written by a
run on the same horizontal grid and partition but with `Nz_old` levels of the same `ExponentialDiscretization`.
Each ocean column is interpolated linearly in z between the old cell centers, using only the old wet cells (S > 0);
new cells below the deepest old wet cell take its value, cells above the top old center take the top value. Columns
that were dry in the old run keep the state `model` already has. The sea ice fields are copied. The clock is not
touched, so the state starts at the model's own `start_date`.
"""
function set_from_checkpoint!(model, directory; iteration, Nz_old = 100, depth = 6000, rank = 0)
    path = joinpath(directory, "checkpoint_rank$(rank)_iteration$(iteration).jld2")
    ocean = model.ocean.model
    sea_ice = model.sea_ice.model
    z_old = reference_centers(Nz_old, depth)
    z_new = reference_centers(size(ocean.grid, 3), depth)

    jldopen(path, "r") do file
        old(name) = file["simulation/model/" * name * "/data"]
        S = old("ocean/model/tracers/S")
        Hz = (size(S, 3) - Nz_old) ÷ 2
        wet = S[8:end-7, 8:end-7, Hz+1:Hz+Nz_old] .> 0
        @info "Rank $rank: old checkpoint T $(size(S)), vertical halo $Hz, $(count(any(wet, dims = 3))) wet columns"

        for (name, field) in (("ocean/model/tracers/T", ocean.tracers.T),
                              ("ocean/model/tracers/S", ocean.tracers.S),
                              ("ocean/model/tracers/e", ocean.tracers.e),
                              ("ocean/model/velocities/u", ocean.velocities.u),
                              ("ocean/model/velocities/v", ocean.velocities.v))
            data = name == "ocean/model/tracers/S" ? S : old(name)
            new = Array(interior(field))
            Nx, Ny = size(new, 1), size(new, 2)
            columns = data[8:7+Nx, 8:7+Ny, Hz+1:Hz+Nz_old]
            Threads.@threads for j in 1:Ny
                for i in 1:Nx
                    interpolate_column!(view(new, i, j, :), view(columns, i, j, :),
                                        view(wet, min(i, size(wet, 1)), min(j, size(wet, 2)), :), z_old, z_new)
                end
            end
            interior(field) .= on_architecture(architecture(field), new)
        end

        surface_fields = (("ocean/model/free_surface/displacement", ocean.free_surface.displacement),
                          ("sea_ice/model/ice_thickness", sea_ice.ice_thickness),
                          ("sea_ice/model/ice_concentration", sea_ice.ice_concentration),
                          ("sea_ice/model/snow_thickness", sea_ice.snow_thickness),
                          ("sea_ice/model/velocities/u", sea_ice.velocities.u),
                          ("sea_ice/model/velocities/v", sea_ice.velocities.v))
        for (name, field) in surface_fields
            haskey(file, "simulation/model/" * name * "/data") || continue
            data = old(name)
            Nx, Ny = size(interior(field), 1), size(interior(field), 2)
            k = (size(data, 3) + 1) ÷ 2  # the single level, between its vertical halos
            interior(field) .= on_architecture(architecture(field), reshape(data[8:7+Nx, 8:7+Ny, k], Nx, Ny, 1))
        end
    end

    foreach(fill_halo_regions!, (ocean.tracers.T, ocean.tracers.S, ocean.tracers.e, ocean.velocities.u,
                                 ocean.velocities.v, ocean.free_surface.displacement, sea_ice.ice_thickness,
                                 sea_ice.ice_concentration, sea_ice.velocities.u, sea_ice.velocities.v))
    return nothing
end

function reference_centers(Nz, depth)
    grid = RectilinearGrid(CPU(); size = Nz, z = ExponentialDiscretization(Nz, -depth, 0; scale = 1100),
                           topology = (Flat, Flat, Bounded))
    return collect(znodes(grid, Center()))
end

function interpolate_column!(new, old, wet, z_old, z_new)
    kb = findfirst(wet)
    kb === nothing && return nothing
    for (n, z) in enumerate(z_new)
        new[n] = if z ≤ z_old[kb]
            old[kb]
        elseif z ≥ z_old[end]
            old[end]
        else
            k = searchsortedlast(z_old, z)
            w = (z - z_old[k]) / (z_old[k+1] - z_old[k])
            (1 - w) * old[k] + w * old[k+1]
        end
    end
    return nothing
end
