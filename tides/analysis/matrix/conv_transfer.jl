# Transfer function of ConvolutionLoadingOperator variants on the 0.36ᵒ aquaplanet TripolarGrid (all cells wet), vs the
# analytic gradient of the degree-n response with Wang PREM load Love numbers (as diag/sal_decomp.jl part A), the
# relative residual, the time per application, and an adjoint check of loading_gradient_transpose!.
#   julia --project=sal_env analysis/matrix/conv_transfer.jl
using NumericalEarth, Oceananigans, CUDA, Printf, Random
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ
const O = NumericalEarth.Oceans
Nx, Ny = 1000, 472; g = 9.80665; R = 6371e3
love_w = let rows = [parse.(Float64, split(l)) for l in eachline(joinpath(@__DIR__, "..", "..", "egbert_ray", "data", "LoveNumbers", "wang-prem-lln.txt")) if !startswith(strip(l), '#') && !isempty(strip(l))]
    rows = filter(r -> r[1] ≤ 400, rows); [1 + r[3] - r[2] for r in rows]; end
γ(L, n) = 3 * 1035 / 5517 * L / (2n + 1)
grid = TripolarGrid(GPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
cg = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo = (7, 7, 4), z = (-6000, 0))
λa = [deg2rad(λnode(i, j, 1, cg, Center(), Center(), Center())) for i in 1:Nx, j in 1:Ny]
φa = [deg2rad(φnode(i, j, 1, cg, Center(), Center(), Center())) for i in 1:Nx, j in 1:Ny]
area = [Azᶜᶜᶜ(i, j, 1, cg) for i in 1:Nx, j in 1:Ny]
function Pnm(n, m, x)
    p = sqrt(1 / 2); c = sqrt(1 - x^2)
    for k in 1:m; p *= sqrt((2k + 1) / 2k) * c; end
    m == n && return p
    p2, p1 = 0.0, p
    for l in m+1:n
        a = sqrt((4l^2 - 1) / (l^2 - m^2)); b = sqrt(((l - 1)^2 - m^2) / (4 * (l - 1)^2 - 1))
        p2, p1 = p1, a * (x * p1 - b * p2)
    end
    p1
end
variants = [("f32 (default)", (; FT = Float32)),
            ("f32 dedup", (; FT = Float32, remove_duplicates = true)),
            ("f32 dedup+subtract", (; FT = Float32, remove_duplicates = true, subtract_uniform = true)),
            ("f64 dedup", (; FT = Float64, remove_duplicates = true)),
            ("f64 (default)", (; FT = Float64))]
p = Field{Center, Center, Nothing}(grid); a¹ = Field{Center, Center, Nothing}(grid); a² = Field{Center, Center, Nothing}(grid)
b¹ = Field{Center, Center, Nothing}(grid); b² = Field{Center, Center, Nothing}(grid); q = Field{Center, Center, Nothing}(grid)
println("variant\tn\tm\ttransfer\trel_residual\trel_residual_|φ|<85")
for (label, kw) in variants
    t = @elapsed op = O.ConvolutionLoadingOperator(grid; components = :grid, kw...)
    @printf("# %s: built in %.1f s\n", label, t)
    for n in (1, 2, 3, 5, 10, 40)
        m = n ÷ 2; δ = 1e-4
        Y = [Pnm(n, m, sin(φ)) * cos(m * λ) for (λ, φ) in zip(λa, φa)]
        set!(p, reshape(g .* Y, Nx, Ny, 1))
        # expected along-grid components: finite differences of the expected potential at the operator's cell positions
        O.loading_gradient!(a¹, a², op, p)
        E, N = Array(interior(a¹, :, :, 1)), Array(interior(a², :, :, 1))
        bad = count(.!isfinite.(E) .| .!isfinite.(N)); E = ifelse.(isfinite.(E), E, 0.0); N = ifelse.(isfinite.(N), N, 0.0)
        # compare via the potential's analytic gradient projected on e¹/e² is awkward; use the geographic operator's
        # convention: grid components equal geographic ones away from the fold (|φ| < 55ᵒ there the grid is lat-lon)
        Yφ = [(Pnm(n, m, sin(φ + δ)) - Pnm(n, m, sin(φ - δ))) / 2δ * cos(m * λ) for (λ, φ) in zip(λa, φa)]
        Yλ = [-m * Pnm(n, m, sin(φ)) * sin(m * λ) / cos(φ) for (λ, φ) in zip(λa, φa)]
        Ee = -γ(love_w[n+1], n) * g / R .* Yλ; Ne = -γ(love_w[n+1], n) * g / R .* Yφ
        w = ifelse.(φa .< deg2rad(50), area, 0.0)   # lat-lon part of the tripolar grid
        T = sum(w .* (E .* Ee .+ N .* Ne)) / sum(w .* (Ee .^ 2 .+ Ne .^ 2))
        res = sqrt(sum(w .* ((E .- T .* Ee) .^ 2 .+ (N .- T .* Ne) .^ 2)) / sum(w .* (Ee .^ 2 .+ Ne .^ 2))) / abs(T)
        # whole-sphere residual magnitude: rms of the output over all cells relative to rms over the lat-lon part
        big = maximum(abs, E)
        @printf("%s\t%d\t%d\t%.4f\t%.4f\tmax|E|/rms %.1f, non-finite %d\n", label, n, m, T, res, big / sqrt(sum(w .* Ee .^ 2) / sum(w)), bad)
        flush(stdout)
    end
    CUDA.synchronize(); t = @elapsed (for _ in 1:3; O.loading_gradient!(a¹, a², op, p); end; CUDA.synchronize())
    Random.seed!(1); set!(p, randn(Nx, Ny, 1)); set!(b¹, randn(Nx, Ny, 1)); set!(b², randn(Nx, Ny, 1))
    O.loading_gradient!(a¹, a², op, p); O.loading_gradient_transpose!(q, op, b¹, b²)
    fin(x) = ifelse.(isfinite.(x), x, 0.0)
    lhs = sum(fin(Array(interior(a¹))) .* Array(interior(b¹)) .+ fin(Array(interior(a²))) .* Array(interior(b²)))
    rhs = sum(Array(interior(p)) .* fin(Array(interior(q))))
    @printf("# %s: %.3f s per application; adjoint check <Ap,b>/<p,Aᵀb> - 1 = %.2e\n", label, t / 3, lhs / rhs - 1)
    flush(stdout)
end
