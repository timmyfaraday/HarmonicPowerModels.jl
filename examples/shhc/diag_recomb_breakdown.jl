################################################################################
# Diagnostic: per-harmonic recombination vs chance-penalty breakdown at the      #
# binding bus, and Monte-Carlo validation of the Cantelli/Gaussian operating     #
# points (true violation probability vs ε). SOC-only (fast; no QCQP).            #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Distributions, Random
const HPM = HarmonicPowerModels
const PMs = PowerModels
include(joinpath(@__DIR__, "pce_params.jl"))

ipopt_ma27 = optimizer_with_attributes(Ipopt.Optimizer,
                "linear_solver" => "ma27", "hsllib" => HSL_jll.libhsl_path,
                "print_level" => 0, "tol" => 1e-10, "constr_viol_tol" => 1e-10)
mosek_soc = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

data = PMs.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc.m"))
H = [1,5,7,11,13]; harm = [5,7,11,13]
α̂, β̂ = PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]
ε = PCE_PARAMS["epsilon"]
base_ranges = Dict(h => PCE_PARAMS["harmonics"][h] for h in harm)
N = 4
all_units = sort(parse.(Int, collect(keys(data["load"]))))
units = all_units[1:N]
get_I(sol,h) = (s=sol["solution"]; haskey(s,"nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h])
λc = HPM.cantelli_lambda(ε); λg = HPM.gaussian_lambda(ε)

hd0 = HPM.replicate(deepcopy(data), H=H)
Z, bus_pos = HPM.reduced_harmonic_impedance(hd0)
ubus  = Dict(d => hd0["nw"]["1"]["load"]["$d"]["load_bus"] for d in units)
nstar = first(b for b in (ubus[d] for d in units) if haskey(bus_pos,b))
ang = HPM.constructive_worst_case_angles(Z, bus_pos, hd0, nstar)
unc = HPM.recentre_angle_ranges(base_ranges, ang, α̂, β̂)

function solve_soc(bound)
    p = Dict{String,Any}(deepcopy(PCE_PARAMS)...)
    p["harmonics"]=unc; p["harmonic_units"]=units; p["bound"]=bound
    hd = HPM.replicate(deepcopy(data), H=H)
    sol = HPM.solve_shhc(hd, sHHC_SOC, mosek_soc; pce_params=p, hpf_optimizer=ipopt_ma27)
    return sol, hd
end
sol_c, hd_c = solve_soc("cantelli")
sol_g, _    = solve_soc("gaussian")

pce = hd_c["pce"]; alpha_V = pce["alpha_V"]; psi = pce["psi_norms"]; K = pce["K"]
fund = hd_c["nw"]["1"]
bus_idx = Dict(b["bus_i"]=>i for (i,(_,b)) in enumerate(sort(collect(fund["bus"]); by=x->x[2]["bus_i"])))
n3 = bus_idx[nstar]
v1 = fund["bus"]["$nstar"]["vm"]
np3 = bus_pos[nstar]

println("="^78)
@printf("IEEE-30 N=%d, binding bus %d (IHD).  V1=%.4f  λ_cant=%.3f λ_gauss=%.3f\n", N, nstar, v1, λc, λg)

# ── per-harmonic recombination vs penalty (per unit hosting current, I_h=1) ──────
println("\nPER-HARMONIC at bus $nstar (values are |V|²/I_h², so I_h-independent):")
@printf("  %-4s %-12s %-12s %-12s %-10s %-10s %-10s %-10s\n",
        "h","Wdet(worst)","E[W]","recomb%","σ(W)","λcσ/Wdet%","λgσ/Wdet%","ihd lim/I²")
for h in harm
    A = alpha_V[h]
    Wdet = (sum(abs(Z[h][np3, bus_pos[ubus[d]]]) for d in units if haskey(bus_pos,ubus[d])))^2
    EW = A[n3,1]
    σ  = sqrt(sum(psi[k+1]*A[n3,idx]^2 for (idx,k) in enumerate(K) if k!=0; init=0.0))
    ihdmax = hd_c["nw"]["$h"]["bus"]["$nstar"]["ihdmax"]
    @printf("  %-4d %-12.4e %-12.4e %-+12.3f %-10.3e %-10.3f %-10.3f %-10.4e\n",
            h, Wdet, EW, (1-EW/Wdet)*100, σ, λc*σ/Wdet*100, λg*σ/Wdet*100, (ihdmax*v1)^2)
end
println("  recomb% = mean reduction from worst case; λσ/Wdet% = chance penalty.")
println("  Uplift sign per h ≈ recomb% − (λσ/Wdet%) scaled; Cantelli penalty ≫ Gaussian.")

# ── Monte-Carlo true violation at each operating point ──────────────────────────
Random.seed!(1); ns = 400_000; Beta = Distributions.Beta(α̂,β̂)
function mc(g)
    out = Dict{Int,Float64}()
    for h in harm
        lim=(hd_c["nw"]["$h"]["bus"]["$nstar"]["ihdmax"]*v1)^2; Ih=sqrt(g[h]); v=0
        zs=Dict(d=>(haskey(bus_pos,ubus[d]) ? Z[h][np3,bus_pos[ubus[d]]] : 0im) for d in units)
        lo=Dict(d=>unc[h][d][1] for d in units); wd=Dict(d=>unc[h][d][2]-unc[h][d][1] for d in units)
        for _ in 1:ns
            Vr=0.0;Vi=0.0
            for d in units
                φ=lo[d]+wd[d]*rand(Beta); zr,zi=real(zs[d]),imag(zs[d]); c,s=cos(φ),sin(φ)
                Vr+=Ih*(zr*c-zi*s); Vi+=Ih*(zi*c+zr*s)
            end
            (Vr^2+Vi^2)>lim && (v+=1)
        end
        out[h]=v/ns
    end
    out
end
gc=Dict(h=>get_I(sol_c,h)^2 for h in harm); gg=Dict(h=>get_I(sol_g,h)^2 for h in harm)
mc_c=mc(gc); mc_g=mc(gg)
@printf("\nMonte-Carlo true P[IHD violated] at bus %d (target ε=%.2f), n=%d:\n", nstar, ε, ns)
@printf("  %-4s  %-14s  %-14s\n","h","Cantelli op.pt","Gaussian op.pt")
for h in harm; @printf("  %-4d  %-14.4f  %-14.4f\n", h, mc_c[h], mc_g[h]); end
@printf("\nHosting totals: SOC-Cantelli=%.5e  SOC-Gaussian=%.5e  (Gauss/Cant %+.2f%%)\n",
        sum(get_I(sol_c,h) for h in harm), sum(get_I(sol_g,h) for h in harm),
        (sum(get_I(sol_g,h) for h in harm)/sum(get_I(sol_c,h) for h in harm)-1)*100)
println("="^78)
