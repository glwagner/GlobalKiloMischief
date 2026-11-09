# EOT20 M2 on the 0.36ᵒ tripolar grid for the tide-differentiability session: bilinear as in score_matrix.jl, in Greenwich
# phase (Z_G = A e^{-iG}) and converted to the model's raw-time convention (η = Re Z e^{iωt}, t from 2015-01-01) for a fit
# window centred on day 45, plus the deep and comparison masks used for scoring.
using JLD2, Dates
const NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"), "NCDatasets"))
here = @__DIR__
depth, λ, φ, area = load(joinpath(here, "..", "diag", "grid_1000.jld2"), "depth", "λ", "φ", "area")
mom6depth = load(joinpath(here, "..", "..", "GlobalKiloMischief", "tides", "mom6_p04arc9_on_1000.jld2"), "depth")
eot = NCDatasets.NCDataset(joinpath(here, "..", "..", "egbert_ray", "data", "EOT20", "ocean_tides", "M2_ocean_eot20.nc"))
lo, la = Float64.(eot["lon"][:]), Float64.(eot["lat"][:])
EZ = map((a, g) -> ismissing(a) || ismissing(g) ? complex(NaN) : a / 100 * cis(-deg2rad(g)), eot["amplitude"][:, :], eot["phase"][:, :]); close(eot)
function bilinear(x0, y0)
    x = mod(x0 - lo[1], 360) / (lo[2] - lo[1]) + 1; y = (y0 - la[1]) / (la[2] - la[1]) + 1
    i = clamp(floor(Int, x), 1, length(lo) - 1); j = clamp(floor(Int, y), 1, length(la) - 1); a, b = x - i, y - j
    (1 - a) * (1 - b) * EZ[i, j] + a * (1 - b) * EZ[i+1, j] + (1 - a) * b * EZ[i, j+1] + a * b * EZ[i+1, j+1]
end
Z_greenwich = bilinear.(λ, φ)
jc(d) = (Dates.datetime2julian(d) - 2451545.0) / 36525
T₀ = jc(DateTime(2015, 1, 1)); V₀ = mod(2(280.4664567 + 36000.76983T₀) - 2(218.3164477 + 481267.88123421T₀), 360)
N = 125.04452 - 1934.136261jc(DateTime(2015, 1, 1) + Day(45))
u = -2.14sind(N); f = 1.0004 - 0.0373cosd(N) + 0.0002cosd(2N)
Z_model_convention = @. abs(Z_greenwich) * f * cis(angle(Z_greenwich) + deg2rad(V₀ + u))   # inverse of score_matrix.jl greenwich()
wet = (depth .> 0) .& (mom6depth .> 0)
compared = wet .& isfinite.(Z_greenwich)
deep = compared .& (depth .> 1000) .& (abs.(φ) .< 66)
jldsave(joinpath(here, "eot20_m2_on_1000.jld2"); Z_greenwich, Z_model_convention, compared, deep, area, λ, φ, V₀, u, f,
        note = "deep uses the native (ETOPO) depth > 1000 m and wet in both depths, exactly as score_matrix.jl; score = sqrt(Σ area |ΔZ|² / Σ area) / √2 over the mask, as score_matrix.jl")
println("saved: $(count(compared)) compared, $(count(deep)) deep cells")
