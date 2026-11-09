# TPXO9 M2 and TICON-3 M2 gauges on the 0.36ᵒ tripolar grid for the tide-differentiability session, exactly as
# score_matrix.jl uses them (included with no runs to score), in Greenwich phase and in the model's raw-time convention
# (η = Re Z e^{iωt}, t from 2015-01-01) for a fit window centred on day 45, as export_eot20.jl does for EOT20.
include(joinpath(@__DIR__, "score_matrix.jl"))
N = 125.04452 - 1934.136261jc(DateTime(2015, 1, 1) + Day(45))
u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
to_model(Z) = @. abs(Z) * f * cis(angle(Z) + deg2rad(V₀ + u))   # inverse of greenwich()
compared = masks[:TPXO9]
deep = compared .& (depth .> 1000) .& (abs.(φ) .< 66)
jldsave(joinpath(@__DIR__, "tpxo9_m2_on_1000.jld2"); Z_greenwich = TPXO, Z_model_convention = to_model(TPXO), compared, deep, area,
        λ, φ, V₀, u, f, note = "TPXO9 nearest p04arc9 cell within 0.5ᵒ, as score_matrix.jl; deep = compared & native depth > 1000 m & |φ| < 66ᵒ")
jldsave(joinpath(@__DIR__, "ticon3_m2_gauges_on_1000.jld2"); Z_greenwich = gZ, Z_model_convention = to_model(gZ), λ = gλ, φ = gφ,
        nearest_cell = nearest, valid = gvalid, deep_200m = gdeep, Nx = size(depth, 1),
        note = "TICON-3 coastal M2 gauges; nearest_cell is a linear index into the 1000×472 grid (0 = no wet cell within 1.5 cells); " *
               "score = 100 sqrt(mean |Z_model[nearest_cell] - Z|² / 2) over valid (all) or deep_200m (model depth > 200 m), in cm")
println("TPXO9: $(count(compared)) compared, $(count(deep)) deep; gauges: $(count(gvalid)) valid, $(count(gdeep)) > 200 m")
