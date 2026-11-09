# Offline attribution of SAL differences, for one finished *_harmonics run (GPU job, M2 fields only):
#  A (once, PART_A=true): aquaplanet transfer function of each discrete operator against the analytic answer, Y_n^m input.
#  B: SAL force at velocity points from every method applied to the SAME M2 η: size, complex projection on the PREM
#     convolution, coastal vs deep; symmetric-kernel check Im∫Z* Ψ[Z]; spectra by degree.
#  C (if the run saved Zu, Zv): M2 energy budget with the discrete operators, including Coriolis and pressure work.
# julia -t 4 --project=sal_env analysis/diag/sal_decomp.jl CUDA <run name without .jld2> <χ>
using NumericalEarth, Oceananigans, Oceananigans.Units, Tidejinks, SpeedyTransforms, ConservativeRegridding
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnode, φnode, peripheral_node
using Oceananigans.Operators
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using JLD2, Printf, Dates, LinearAlgebra
import SPICE
using CUDA
const O = NumericalEarth.Oceans
root = joinpath(@__DIR__, "..", "..", "GlobalKiloMischief", "tides")
run_name, χ = ARGS[2], parse(Float64, ARGS[3])
g = 9.80665; R = 6371e3; ρ = 1035.0
Nx, Ny = 1000, 472; halo = (7, 7, 4); Δt = 600.0
love_w = let rows = [parse.(Float64, split(l)) for l in eachline(joinpath(root, "..", "..", "egbert_ray", "data", "LoveNumbers", "wang-prem-lln.txt")) if !startswith(strip(l), '#') && !isempty(strip(l))]
    rows = filter(r -> r[1] ≤ 400, rows); [1 + r[3] - r[2] for r in rows]; end
love_g = [ (r = parse.(Float64, split(l)); 1 + r[3] - r[2]) for l in eachline(joinpath(root, "..", "..", "egbert_ray", "data", "LoveNumbers", "gegout-lln.txt")) if !startswith(strip(l), '#') && !isempty(strip(l))][1:401]
fit(n) = O.one_plus_k_minus_h(O.AsymptoticLoadLoveNumbers(), n)
γ(L, n) = 3 * 1035 / 5517 * L / (2n + 1)
println("# n  1+k-h: Wang  Gegout  fit  fit/3  ratio fit/3 to Wang;  γ_n Wang, γ_n fit/3")
for n in (1, 2, 3, 4, 5, 8, 10, 15, 20, 30, 40, 80, 200)
    @printf("# %3d  %.3f %.3f %.3f %.3f  %.3f   %.4f %.4f\n", n, love_w[n+1], love_g[n+1], fit(n), fit(n) / 3, fit(n) / 3 / love_w[n+1], γ(love_w[n+1], n), γ(fit(n) / 3, n))
end
flush(stdout)

arch = GPU()
if get(ENV, "PART_A", "false") == "true"
    gridA = TripolarGrid(arch; size = (Nx, Ny, 1), halo, z = (-6000, 0))
    cA = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-6000, 0))
    λa = [deg2rad(λnode(i, j, 1, cA, Center(), Center(), Center())) for i in 1:Nx, j in 1:Ny]
    φa = [deg2rad(φnode(i, j, 1, cA, Center(), Center(), Center())) for i in 1:Nx, j in 1:Ny]
    area = [Azᶜᶜᶜ(i, j, 1, cA) for i in 1:Nx, j in 1:Ny]
    function Pnm(n, m, x)
        p = sqrt(1 / 2) ; c = sqrt(1 - x^2)
        for k in 1:m; p *= sqrt((2k + 1) / 2k) * c; end
        m == n && return p
        p2, p1 = 0.0, p
        for l in m+1:n
            a = sqrt((4l^2 - 1) / (l^2 - m^2)); b = sqrt(((l - 1)^2 - m^2) / (4 * (l - 1)^2 - 1))
            p2, p1 = p1, a * (x * p1 - b * p2)
        end
        return p1
    end
    ops = get(ENV, "PART_A_ONLY", "false") == "true" ? (prem32 = ConvolutionLoadingOperator(gridA; FT = Float32), prem64 = ConvolutionLoadingOperator(gridA; FT = Float64), paper64 = ConvolutionLoadingOperator(gridA; FT = Float64, love_numbers = AsymptoticLoadLoveNumbers(), self_term = false, earth_density = 3 * 5517)) : (prem = ConvolutionLoadingOperator(gridA; FT = Float32), prem_noself = ConvolutionLoadingOperator(gridA; FT = Float32, self_term = false),
           paper = ConvolutionLoadingOperator(gridA; FT = Float32, love_numbers = AsymptoticLoadLoveNumbers(), self_term = false, earth_density = 3 * 5517))
    sh = (n40g24 = SphericalHarmonicSelfAttractionAndLoading(gridA; truncation = 40, nlat_half = 24), n40g128 = SphericalHarmonicSelfAttractionAndLoading(gridA; truncation = 40, nlat_half = 128))
    p = Field{Center, Center, Nothing}(gridA); a¹ = Field{Center, Center, Nothing}(gridA); a² = Field{Center, Center, Nothing}(gridA)
    println("n\tm\top\texpected_love_ratio\ttransfer_vs_WangPREM\trel_residual")
    for n in (1, 2, 3, 5, 10, 40)
        m = n ÷ 2
        δ = 1e-4
        Y = [Pnm(n, m, sin(φ)) * cos(m * λ) for (λ, φ) in zip(λa, φa)]
        Yφ = [(Pnm(n, m, sin(φ + δ)) - Pnm(n, m, sin(φ - δ))) / 2δ * cos(m * λ) for (λ, φ) in zip(λa, φa)]
        Yλ = [abs(φ) < deg2rad(85) ? -m * Pnm(n, m, sin(φ)) * sin(m * λ) / cos(φ) : 0.0 for (λ, φ) in zip(λa, φa)]
        pole = abs.(φa) .< deg2rad(85); E₀ = ifelse.(isfinite.(area) .& pole, area, 0.0)
        set!(p, reshape(g .* Y, Nx, Ny, 1))
        for (name, op) in pairs(ops)
            loading_gradient!(a¹, a², op, p)
            E, N = Array(interior(a¹, :, :, 1)), Array(interior(a², :, :, 1)); E = ifelse.(isfinite.(E), E, 0.0); N = ifelse.(isfinite.(N), N, 0.0)
            Ee = -γ(love_w[n+1], n) * g / R .* Yλ; Ne = -γ(love_w[n+1], n) * g / R .* Yφ
            T = sum(E₀ .* (E .* Ee .+ N .* Ne)) / sum(E₀ .* (Ee .^ 2 .+ Ne .^ 2))
            res = sqrt(sum(E₀ .* ((E .- T .* Ee) .^ 2 .+ (N .- T .* Ne) .^ 2)) / sum(E₀ .* (Ee .^ 2 .+ Ne .^ 2))) / abs(T)
            @printf("%d\t%d\t%s\t%.3f\t%.4f\t%.4f\n", n, m, name, name in (:paper, :paper64) ? fit(n) / 3 / love_w[n+1] : 1.0, T, res)
        end
        if n ≤ 40
            for (name, s) in pairs(sh)
                O.apply_loading_operator!(a¹, s.operator, p)
                Ψ = Array(interior(a¹, :, :, 1)); Ψe = -γ(love_g[n+1], n) * g .* Y
                T = sum(area .* Ψ .* Ψe) / sum(area .* Ψe .^ 2)
                @printf("%d\t%d\t%s\t1.000\t%.4f\t%.4f\n", n, m, name, T, sqrt(sum(area .* (Ψ .- T .* Ψe) .^ 2) / sum(area .* Ψe .^ 2)) / abs(T))
            end
        end
        flush(stdout)
    end
    get(ENV, "PART_A_ONLY", "false") == "true" && exit(0)
end

cd(root); SPICE.furnsh("kernels.txt")
cpu_grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), halo, z = (-6000, 0))
bottom_height = regrid_bathymetry(cpu_grid; minimum_depth = 20, interpolation_passes = 3, major_basins = 1)
let zb = interior(bottom_height, :, :, 1)
    for pass in 1:10
        H = -Array(zb); n = 0
        for j in 2:Ny-1, i in 1:Nx
            H[i, j] > 0 || continue
            nb = filter(>(0), (H[mod1(i-1, Nx), j], H[mod1(i+1, Nx), j], H[i, j-1], H[i, j+1]))
            if !isempty(nb) && all(>(5 * H[i, j]), nb); zb[i, j] = 0; n += 1; end
        end
        n == 0 && break
    end
end
function make(a)
    ug = TripolarGrid(a; size = (Nx, Ny, 1), halo, z = (-6000, 0))
    zb = Field{Center, Center, Nothing}(ug); set!(zb, Array(interior(bottom_height))); fill_halo_regions!(zb)
    return ImmersedBoundaryGrid(ug, PartialCellBottom(zb; minimum_fractional_cell_height = 1e-3); active_cells_map = true)
end
grid = make(arch); gridc = make(CPU())
c(i, j, f, l) = f(i, j, 1, gridc, l...)
λ = [mod(λnode(i, j, 1, gridc, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, gridc, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
H = -Array(interior(bottom_height, :, :, 1)); wet = H .> 0; area = [Azᶜᶜᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny]
wetu = [!peripheral_node(i, j, 1, gridc, Face(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
wetv = [!peripheral_node(i, j, 1, gridc, Center(), Face(), Center()) for i in 1:Nx, j in 1:Ny]
Vu = [Azᶠᶜᶜ(i, j, 1, gridc) * Δzᶠᶜᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny] .* wetu
Vv = [Azᶜᶠᶜ(i, j, 1, gridc) * Δzᶜᶠᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny] .* wetv
Au = [Azᶠᶜᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny] .* wetu; Av = [Azᶜᶠᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny] .* wetv
coast = falses(Nx, Ny)  # within 2 cells of land, equatorward of 66ᵒ
for j in 3:Ny-2, i in 1:Nx
    wet[i, j] || continue
    any(!wet[mod1(i + a, Nx), j + b] for a in -2:2, b in -2:2) && abs(φ[i, j]) < 66 && (coast[i, j] = true)
end
masks = (deep = (H .> 1000) .& (abs.(φ) .< 66) .& .!coast, coast = coast, all = wet)
shift_u(m) = m .| m[[mod1(i - 1, Nx) for i in 1:Nx], :]; shift_v(m) = m .| m[:, [max(j - 1, 1) for j in 1:Ny]]

h = load(joinpath("output", run_name * ".jld2"))
Z = ComplexF64.(h["Z"]["M2"]) .* wet
println("# wet cells ", count(wet), "; run ", run_name, "; samples ", h["samples"])

# --- spherical harmonics by direct quadrature on the model grid, as MOM6 does
function analysis(f, L)
    a = [zeros(ComplexF64, L + 1, L + 1) for _ in 1:Threads.maxthreadid()]; b = deepcopy(a)
    cells = findall(wet)
    Threads.@threads :static for q in cells
        A, B = a[Threads.threadid()], b[Threads.threadid()]
        x = sind(φ[q]); y = cosd(φ[q]); w = f[q] * area[q] / R^2
        pmm = sqrt(1 / 4π)
        for m in 0:L
            m > 0 && (pmm *= sqrt((2m + 1) / 2m) * y)
            s = m == 0 ? 1.0 : sqrt(2.0); cc, ss = w * s * cosd(m * λ[q]), w * s * sind(m * λ[q])
            p2, p1 = 0.0, pmm
            A[m+1, m+1] += p1 * cc; B[m+1, m+1] += p1 * ss
            for l in m+1:L
                e = sqrt((4l^2 - 1) / (l^2 - m^2)); d = sqrt(((l - 1)^2 - m^2) / (4 * (l - 1)^2 - 1))
                p2, p1 = p1, e * (x * p1 - d * p2)
                A[l+1, m+1] += p1 * cc; B[l+1, m+1] += p1 * ss
            end
        end
    end
    return sum(a), sum(b)
end
function synthesis(A, B, weights, L)
    out = zeros(ComplexF64, Nx, Ny)
    Threads.@threads :static for q in findall(wet)
        x = sind(φ[q]); y = cosd(φ[q]); pmm = sqrt(1 / 4π); S = 0.0 + 0im
        for m in 0:min(L, size(A, 2) - 1)
            m > 0 && (pmm *= sqrt((2m + 1) / 2m) * y)
            s = m == 0 ? 1.0 : sqrt(2.0); cc, ss = s * cosd(m * λ[q]), s * sind(m * λ[q])
            p2, p1 = 0.0, pmm
            S += weights[m+1] * p1 * (A[m+1, m+1] * cc + B[m+1, m+1] * ss)
            for l in m+1:L
                e = sqrt((4l^2 - 1) / (l^2 - m^2)); d = sqrt(((l - 1)^2 - m^2) / (4 * (l - 1)^2 - 1))
                p2, p1 = p1, e * (x * p1 - d * p2)
                S += weights[l+1] * p1 * (A[l+1, m+1] * cc + B[l+1, m+1] * ss)
            end
        end
        out[q] = S
    end
    return out
end
Lmax = parse(Int, get(ENV, "LMAX", "200"))
@time A, B = analysis(Z, Lmax)
power = [sum(abs2, @view(A[n+1, :])) + sum(abs2, @view(B[n+1, :])) for n in 0:Lmax]
weights(L, nmax) = [n > nmax ? 0.0 : γ(L(n), n) for n in 0:Lmax]
Wang(n) = love_w[n+1]; Geg(n) = love_g[n+1]; Paper(n) = fit(n) / 3
println("n\tpower_Z_frac\tpsi_Wang\tpsi_Paper_over_Wang\tcum_diff_frac")
tot = sum(power[2:end])
dcum = 0.0; dall = sum(power[n+1] * (γ(Paper(n), n) - γ(Wang(n), n))^2 for n in 1:Lmax)
for n in 1:Lmax
    global dcum += power[n+1] * (γ(Paper(n), n) - γ(Wang(n), n))^2
    (n ≤ 10 || n % 10 == 0) && @printf("%d\t%.4f\t%.4f\t%.3f\t%.3f\n", n, power[n+1] / tot, γ(Wang(n), n), γ(Paper(n), n) / γ(Wang(n), n), dcum / dall)
end
# Effective scalar β: Ψ/(gη) projected on η over the ocean, for each kernel
Ψs = Dict{String, Matrix{ComplexF64}}()
for (name, L, nmax) in (("SHquad n40 Gegout", Geg, 40), ("SHquad n40 Wang", Wang, 40), ("SHquad n200 Wang", Wang, Lmax), ("SHquad n40 paper/3", Paper, 40), ("SHquad n200 paper/3", Paper, Lmax))
    Ψs[name] = -g .* synthesis(A, B, weights(L, nmax), Lmax)
    ip = sum(conj.(Z) .* Ψs[name] .* area)
    # Ψ = - β g η: β_eff = - <Z, Ψ> / g <Z, Z>. Imag ≠ 0 would mean net work by a symmetric kernel.
    @printf("# %s: β_eff = %.4f, Im part %.2e; over deep %.4f\n", name, -real(ip) / (g * sum(abs2.(Z) .* area)), imag(ip) / (g * sum(abs2.(Z) .* area)),
            -real(sum((conj.(Z) .* Ψs[name] .* area)[masks.deep])) / (g * sum((abs2.(Z) .* area)[masks.deep])))
end
flush(stdout)
Ψs["scalar β=0.1"] = -0.1 * g .* Z

# --- model operators on the same η
pfield(x) = (f = Field{Center, Center, Nothing}(grid); set!(f, reshape(x, Nx, Ny, 1)); fill_halo_regions!(f); f)
cc(x) = (f = Field{Center, Center, Nothing}(gridc); set!(f, reshape(x, Nx, Ny, 1)); fill_halo_regions!(f); f)
function grad_faces(Ψ) # -∂Ψ at faces as for the model's potential forcing
    Fu = zeros(ComplexF64, Nx, Ny); Fv = zeros(ComplexF64, Nx, Ny)
    for (part, unit) in ((real, 1.0 + 0im), (imag, 0.0 + 1im))
        f = cc(part.(Ψ))
        for j in 1:Ny, i in 1:Nx
            wetu[i, j] && (Fu[i, j] -= unit * ∂xᶠᶜᶜ(i, j, 1, gridc, f))
            wetv[i, j] && (Fv[i, j] -= unit * ∂yᶜᶠᶜ(i, j, 1, gridc, f))
        end
    end
    return Fu, Fv
end
force = Dict{String, Any}()
for (k, v) in Ψs; force[k] = grad_faces(v); end
force["pressure g∇η"] = grad_faces(g .* Z)
sh40 = SphericalHarmonicSelfAttractionAndLoading(grid; truncation = 40, nlat_half = 24)
shp = Field{Center, Center, Nothing}(grid)
ΨSH = zeros(ComplexF64, Nx, Ny)
for (part, unit) in ((real, 1.0 + 0im), (imag, 0.0 + 1im))
    O.apply_loading_operator!(shp, sh40.operator, pfield(g .* part.(Z)))
    ΨSH .+= unit .* Array(interior(shp, :, :, 1)) .* wet
end
Ψs["model SH n40 (Gaussian N24)"] = ΨSH
force["model SH n40 (Gaussian N24)"] = grad_faces(ΨSH)
@printf("# model SH symmetric check: Im<Z,Ψ>/<Z,gZ> = %.2e, β_eff = %.4f\n", imag(sum(conj.(Z) .* ΨSH .* area)) / (g * sum(abs2.(Z) .* area)), -real(sum(conj.(Z) .* ΨSH .* area)) / (g * sum(abs2.(Z) .* area)))
function conv_force(sal)
    Fu = zeros(ComplexF64, Nx, Ny); Fv = zeros(ComplexF64, Nx, Ny)
    for (part, unit) in ((real, 1.0 + 0im), (imag, 0.0 + 1im))
        loading_gradient!(sal.gradient..., sal.operator, pfield(g .* part.(Z)))
        a = cc(Array(interior(sal.gradient[1], :, :, 1)) .* wet); b = cc(Array(interior(sal.gradient[2], :, :, 1)) .* wet)
        for j in 1:Ny, i in 1:Nx
            wetu[i, j] && (Fu[i, j] -= unit * (a[i-1, j, 1] + a[i, j, 1]) / 2)
            wetv[i, j] && (Fv[i, j] -= unit * (b[i, j-1, 1] + b[i, j, 1]) / 2)
        end
    end
    return Fu, Fv
end
force["conv PREM"] = conv_force(ConvolutionSelfAttractionAndLoading(grid; FT = Float32))
force["conv PREM no self term"] = conv_force(ConvolutionSelfAttractionAndLoading(grid; FT = Float32, self_term = false))
force["conv paper"] = conv_force(ConvolutionSelfAttractionAndLoading(grid; FT = Float32, love_numbers = AsymptoticLoadLoveNumbers(), self_term = false, earth_density = 3 * 5517))
Fref = force["conv PREM"]; Fp = force["pressure g∇η"]
ip(a, b, mu, mv) = sum((conj.(a[1]) .* b[1] .* Au)[mu]) + sum((conj.(a[2]) .* b[2] .* Av)[mv])
println("op\tmask\trms_F(1e-5m/s2)\trms/rms(g∇η)\talpha_vs_convPREM\tphase_deg\trms_diff_from_convPREM\tdiff/rms_ref")
for (k, F) in sort(collect(force), by = first), (mn, m) in pairs(masks)
    mu, mv = shift_u(m), shift_v(m)
    nn = real(ip(F, F, mu, mv)); nr = real(ip(Fref, Fref, mu, mv)); α = ip(Fref, F, mu, mv) / nr; A₀ = sum(Au[mu]) + sum(Av[mv])
    d = (F[1] .- Fref[1], F[2] .- Fref[2])
    @printf("%s\t%s\t%.3f\t%.4f\t%.4f\t%.2f\t%.3f\t%.3f\n", k, mn, 1e5sqrt(nn / A₀), sqrt(nn / real(ip(Fp, Fp, mu, mv))), abs(α), rad2deg(angle(α)), 1e5sqrt(real(ip(d, d, mu, mv)) / A₀), sqrt(real(ip(d, d, mu, mv)) / nr))
end
flush(stdout)
jldsave(joinpath(@__DIR__, "sal_decomp_forces_" * replace(run_name, "barotropic_tides_0.36deg_dt10min_" => "", "_harmonics" => "") * ".jld2");
        Fu = Dict(k => ComplexF32.(v[1]) for (k, v) in force if k in ("conv PREM", "conv paper", "SHquad n40 Wang", "SHquad n200 Wang")),
        Fv = Dict(k => ComplexF32.(v[2]) for (k, v) in force if k in ("conv PREM", "conv paper", "SHquad n40 Wang", "SHquad n200 Wang")), λ = Float32.(λ), φ = Float32.(φ), H = Float32.(H))

# --- M2 energy budget with the discrete operators
haskey(h, "Zu") || exit(0)
const K_names = ["M2", "N2", "S2", "K1", "O1"]
Zuall = [ComplexF64.(h["Zu"][k][1:Nx, 1:Ny]) for k in K_names]; Zvall = [ComplexF64.(h["Zv"][k][1:Nx, 1:Ny]) for k in K_names]
Zu, Zv = Zuall[1] .* wetu, Zvall[1] .* wetv
W(F, mu = trues(Nx, Ny), mv = trues(Nx, Ny)) = 1e-12 * ρ / 2 * (sum((real.(conj.(Zu) .* F[1]) .* Vu)[mu]) + sum((real.(conj.(Zv) .* F[2]) .* Vv)[mv]))
ω = 2π ./ (3600 .* [12.4206012, 12.65834751, 12.0, 23.9344696, 25.8193417])
ts = h["first_day"] * 86400 .+ 3600.0 .* (0:h["samples"]-1)
D = hcat(ones(length(ts)), [cos(w * t) for t in ts, w in ω], [sin(w * t) for t in ts, w in ω]); Pinv = pinv(D)
Φ = Tidejinks.TidalPotential(gridc, (DateTime(2015, 1, 1), DateTime(2015, 1, 1) + Day(62)); reference_frame = "IAU_EARTH")
aΦ = zeros(Nx, Ny); bΦ = zeros(Nx, Ny)
Threads.@threads for q in findall(wet)
    for (n, t) in enumerate(ts)
        v = Φ(λ[q], φ[q], t); aΦ[q] += Pinv[2, n] * v; bΦ[q] += Pinv[7, n] * v
    end
end
Ftide = grad_faces(aΦ .- im .* bΦ)
# Quadratic drag: u(t) reconstructed from the five fitted constituents, as in the implicit drag -c_d |u| u / Δz; M2-projected
cd₀ = 0.003
Dzu = [Δzᶠᶜᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny]; Dzv = [Δzᶜᶠᶜ(i, j, 1, gridc) for i in 1:Nx, j in 1:Ny]
sum_u = [zeros(Threads.maxthreadid()) for _ in 1:2]
acu = [zeros(ComplexF64, Nx, Ny) for _ in 1:Threads.maxthreadid()]; acv = deepcopy(acu)
rate_ip1 = [mod1(i + 1, Nx) for i in 1:Nx]
total_cubed = zeros(Threads.maxthreadid())
Threads.@threads :static for n in eachindex(ts)
    t = ts[n]; id = Threads.threadid()
    u = sum(real.(Zuall[k] .* cis(ω[k] * t)) for k in 1:5) .* wetu; v = sum(real.(Zvall[k] .* cis(ω[k] * t)) for k in 1:5) .* wetv
    wm = (Pinv[2, n] - im * Pinv[7, n])
    for j in 2:Ny-1, i in 1:Nx
        ip1, im1 = rate_ip1[i], mod1(i - 1, Nx)
        if wetu[i, j]
            vm = (v[im1, j] + v[i, j] + v[im1, j+1] + v[i, j+1]) / 4
            acu[id][i, j] += wm * (-cd₀ * sqrt(u[i, j]^2 + vm^2) * u[i, j] / Dzu[i, j])
            total_cubed[id] += cd₀ * sqrt(u[i, j]^2 + vm^2) * u[i, j]^2 * Au[i, j] * ρ / length(ts)
        end
        if wetv[i, j]
            um = (u[i, j] + u[ip1, j] + u[i, j-1] + u[ip1, j-1]) / 4
            acv[id][i, j] += wm * (-cd₀ * sqrt(v[i, j]^2 + um^2) * v[i, j] / Dzv[i, j])
            total_cubed[id] += cd₀ * sqrt(v[i, j]^2 + um^2) * v[i, j]^2 * Av[i, j] * ρ / length(ts)
        end
    end
end
Fbot = (sum(acu), sum(acv))
C = reshape(reinterpret(Float64, read("wavedrag_halfdeg_C.bin")), 720, 312)
tile(λ, φ) = (j = floor(Int, (φ + 78) / 0.5) + 1; 1 ≤ j ≤ 312 ? C[clamp(floor(Int, mod(λ, 360) / 0.5) + 1, 1, 720), j] : 0.0)
rate = [H[q] > 0 ? min(χ * tile(λ[q], φ[q]) / H[q], 1 / 2Δt) : 0.0 for q in CartesianIndices(H)]
rf = cc(rate)
Fwave = ([-(rf[i-1, j, 1] + rf[i, j, 1]) / 2 * Zu[i, j] for i in 1:Nx, j in 1:Ny], [-(rf[i, j-1, 1] + rf[i, j, 1]) / 2 * Zv[i, j] for i in 1:Nx, j in 1:Ny])
coriolis = HydrostaticSphericalCoriolis()
function cor(sc)
    co = HydrostaticSphericalCoriolis(scheme = sc)
    Fu = zeros(ComplexF64, Nx, Ny); Fv = zeros(ComplexF64, Nx, Ny)
    for (part, unit) in ((real, 1.0 + 0im), (imag, 0.0 + 1im))
        uf = XFaceField(gridc); vf = YFaceField(gridc)
        interior(uf, :, :, 1) .= part.(Zu); interior(vf, 1:Nx, 1:Ny, 1) .= part.(Zv); fill_halo_regions!((uf, vf))
        for j in 2:Ny-1, i in 1:Nx
            wetu[i, j] && (Fu[i, j] -= unit * x_f_cross_U(i, j, 1, gridc, co, (u = uf, v = vf, w = nothing)))
            wetv[i, j] && (Fv[i, j] -= unit * y_f_cross_U(i, j, 1, gridc, co, (u = uf, v = vf, w = nothing)))
        end
    end
    return Fu, Fv
end
# Viscous work by the Laplacian ν ∇²u on the (nearly orthogonal) grid: - ν ∫ H |∇u|² estimated from face differences
Fp2 = (-Fp[1] .* -1, -Fp[2] .* -1)  # force −g∇η is grad_faces(gZ)
println("# M2 budget, TW (+ = into the M2 tide, or by − sign a sink), run ", run_name)
for (k, v) in sort(collect(force), by = first); @printf("# work %-32s %8.4f   deep %8.4f coast %8.4f\n", k, W(v), W(v, shift_u(masks.deep), shift_v(masks.deep)), W(v, shift_u(masks.coast), shift_v(masks.coast))); end
@printf("# work tidal potential %.4f\n# work bottom drag (M2-projected) %.4f; all-constituent mean ρ c_d|u|³ %.4f TW\n# work wave drag %.4f\n", W(Ftide), W(Fbot), 1e-12 * sum(total_cubed), W(Fwave))
let a = cor(Oceananigans.Advection.EnstrophyConserving()), b = cor(Oceananigans.Advection.EnergyConserving()), c = cor(Oceananigans.Coriolis.TriadScheme())
    @printf("# Coriolis rms force enst %.4e, rms(enst-ener) %.4e, rms(enst-triad) %.4e; work triad %.4f\n", sqrt(sum(abs2, a[1]) + sum(abs2, a[2])), sqrt(sum(abs2, a[1] .- b[1]) + sum(abs2, a[2] .- b[2])), sqrt(sum(abs2, a[1] .- c[1]) + sum(abs2, a[2] .- c[2])), W(c)); end
@printf("# work Coriolis EnstrophyConserving(run) %.4f; EnergyConserving %.4f\n", W(cor(Oceananigans.Advection.EnstrophyConserving())), W(cor(Oceananigans.Advection.EnergyConserving())))
