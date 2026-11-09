# Subgrid topographic roughness for a scale-aware wave drag, from ETOPO2022 (60 arc seconds) on 1ᵒ tiles.
#
# Johnston, Shakespeare & Constantinou (2025, arXiv:2503.12845, Eq. 4) give the propagating wave drag on a tidal flow U
# over topography of spectrum |ĥ|² as F = √((N² - ω²)(ω² - f²)) / (4π|ω|) U ∫ |ĥ(K)|² K² dK (isotropic, per unit
# density). Per unit area, with P(k, l) the spectral density of h such that ∫∫ P dk dl = ⟨h²⟩, this is
# F = S N √(1 - f²/ω²) U / 2 for N ≫ ω, S = ∫∫ P K dk dl, which is Jayne & St. Laurent's (2001) κ h² made spectral.
# Only topography finer than the model resolves should drag the barotropic tide (the model generates internal tides at
# resolved scales itself), so S is summed over K > 2π / λc with λc = 2 √(Δx Δy) of the model grid at each latitude.
# On each tile: land is filled with the mean ocean depth, a plane is removed, a Hann window applied and P normalized so
# that its sum is the variance of the detrended depth. Tiles shallower than 100 m on average, or less than half ocean,
# get S = 0. Roughness at wavelengths below ~4 km (twice the ETOPO spacing) is missing.
#
#   julia --project=. tides/subgrid_wave_drag.jl   # writes tides/subgrid_roughness.jld2: S for 1/6ᵒ, 1/12ᵒ and all K

using NCDatasets
using FFTW
using JLD2
using Printf
using Statistics

const R = 6371e3
const points = 60  # 1ᵒ tiles of 1' points

path = joinpath(ENV["NUMERICALEARTH_DATA_DIRECTORY"], "ETOPO", "ETOPO_2022_v1_60s_N90W180_surface.nc")
z = NCDataset(ds -> Float32.(ds["z"][:, :]), path)  # (longitude, latitude), -180 to 180, -90 to 90
@info "Read ETOPO: $(size(z))"

window = [sinpi((i - 0.5) / points)^2 * sinpi((j - 0.5) / points)^2 for i in 1:points, j in 1:points]
x = repeat((1:points) .- (points + 1) / 2, 1, points)
y = repeat(((1:points) .- (points + 1) / 2)', points, 1)
plane = hcat(ones(points^2), vec(x), vec(y))
frequencies = fftfreq(points)  # cycles per point

cutoffs = (sixth = 6, twelfth = 12, all = nothing)
S = Dict(name => zeros(360, 180) for name in keys(cutoffs))

Threads.@threads for J in 1:180
    φ = -90 + J - 0.5
    Δx = R * cosd(φ) * deg2rad(1 / 60)
    Δy = R * deg2rad(1 / 60)
    k = 2π .* frequencies ./ Δx
    l = 2π .* frequencies ./ Δy
    K = [sqrt(a^2 + b^2) for a in k, b in l]
    for I in 1:360
        h = Float64.(z[(I - 1) * points + 1:I * points, (J - 1) * points + 1:J * points])
        ocean = h .< 0
        (count(ocean) < points^2 / 2 || mean(h[ocean]) > -100) && continue
        h[.!ocean] .= mean(h[ocean])
        h .-= reshape(plane * (plane \ vec(h)), points, points)
        P = abs2.(fft(h .* window)) ./ (points^2 * sum(abs2, window))
        for (name, cells_per_degree) in pairs(cutoffs)
            Kc = if isnothing(cells_per_degree)
                0.0
            else
                model_Δx = R * cosd(φ) * deg2rad(1 / cells_per_degree)
                model_Δy = R * deg2rad(1 / cells_per_degree)
                2π / (2 * sqrt(model_Δx * model_Δy))
            end
            S[name][I, J] = sum(P[n] * K[n] for n in eachindex(P) if K[n] > Kc)
        end
    end
end

for name in keys(cutoffs)
    s = S[name][S[name] .> 0]
    @info @sprintf("S %-8s: %d tiles, median %.2f m, 90th percentile %.2f m, max %.1f m", name, length(s), median(s),
                   quantile(s, 0.9), maximum(s))
end
jldsave(joinpath(@__DIR__, "subgrid_roughness.jld2"); S_sixth = S[:sixth], S_twelfth = S[:twelfth], S_all = S[:all],
        longitude = collect(-179.5:1:179.5), latitude = collect(-89.5:1:89.5))
