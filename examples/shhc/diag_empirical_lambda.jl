################################################################################
# Diagnostic: empirical (quantile-calibrated) λ*(ε) for the binding chance        #
# constraints. For each constrained aggregate S, sample its output distribution   #
# and compute λ* = (q_{1−ε}(S) − E[S]) / σ(S). Compare with Cantelli & Gaussian.  #
# For IHD (scalar-linear in g) λ* is g-independent; for THD/RMS it is evaluated    #
# at the Gaussian-SOC solution (one iteration).                                    #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Distributions, Random, Statistics
const HPM = HarmonicPowerModels
const PMs = PowerModels
include(joinpath(@__DIR__, "pce_params.jl"))

const NET  = get(ENV, "SHHC_NET", "industrial")
const NUN  = parse(Int, get(ENV, "SHHC_N", "4"))
const M    = parse(Int, get(ENV, "SHHC_MC", "100000"))
datafile = NET == "ieee30" ? "30bus_network_shhc.m" : "industrial_network_shhc.m"

ipopt = optimizer_with_attributes(Ipopt.Optimizer,"linear_solver"=>"ma27",
            "hsllib"=>HSL_jll.libhsl_path,"print_level"=>0,"tol"=>1e-10,"constr_viol_tol"=>1e-10)
mosek = optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true)

data = PMs.parse_file(joinpath(@__DIR__,"data",datafile))
H=[1,5,7,11,13]; harm=[5,7,11,13]
α̂,β̂=PCE_PARAMS["alpha_hat"],PCE_PARAMS["beta_hat"]; ε=PCE_PARAMS["epsilon"]
base=Dict(h=>PCE_PARAMS["harmonics"][h] for h in harm)
λc=HPM.cantelli_lambda(ε); λg=HPM.gaussian_lambda(ε)
all_units=sort(parse.(Int,collect(keys(data["load"]))))
units = NET=="ieee30" ? all_units[1:NUN] : all_units
get_I(s,h)=(x=s["solution"]; haskey(x,"nw") ? x["nw"]["$h"]["load"]["1"]["cmd"] : x[h])

# physical-mean Gaussian SOC solution (for g and limits)
p=deepcopy(PCE_PARAMS); p["harmonics"]=Dict{Int,Any}(h=>base[h] for h in harm)
p["harmonic_units"]=units; p["bound"]="gaussian"
hd=HPM.replicate(deepcopy(data),H=H)
sol=HPM.solve_shhc(hd,sHHC_SOC,mosek;pce_params=p,hpf_optimizer=ipopt)
g=Dict(h=>get_I(sol,h)^2 for h in harm)

Z,bp=HPM.reduced_harmonic_impedance(hd)
fund=hd["nw"]["1"]
v1=Dict(b["bus_i"]=>b["vm"] for (_,b) in fund["bus"])
ubus=Dict(d=>hd["nw"]["1"]["load"]["$d"]["load_bus"] for d in units)

# sample per-harmonic squared response a_h(ω)[n] = |Σ_d Z[h][n,d] e^{jφ}|²  (per I²)
Random.seed!(1); Bdist=Distributions.Beta(α̂,β̂)
nonref=sort(collect(keys(bp)))
# accumulate samples of a_h[n] only for buses that have an IHD/THD limit
λstar(x) = (quantile(x,1-ε)-mean(x))/std(x)

# storage: per (h, busid) vector of a samples; and per busid THD aggregate S
A=Dict((h,n)=>Vector{Float64}(undef,M) for h in harm, n in nonref)
for m in 1:M
    ω=Dict(d=>rand(Bdist) for d in units)
    for h in harm
        # voltage response per non-ref bus
        for n in nonref
            Vr=0.0;Vi=0.0
            for d in units
                bb=ubus[d]; haskey(bp,bb) || continue
                z=Z[h][bp[n],bp[bb]]; φ=base[h][1]+(base[h][2]-base[h][1])*ω[d]
                zr,zi=real(z),imag(z); c,s=cos(φ),sin(φ)
                Vr+=zr*c-zi*s; Vi+=zi*c+zr*s
            end
            A[(h,n)][m]=Vr*Vr+Vi*Vi
        end
    end
end

println("="^74)
@printf("Empirical λ* on %s (units=%d, M=%d).  Cantelli λ=%.3f  Gaussian λ=%.3f\n",
        NET, length(units), M, λc, λg)

# IHD per (bus,h): λ* g-independent. report the binding (max utilisation) ones
println("\nIHD (per bus,harmonic) — λ* = (q95−E)/σ of |V_{n,h}|²:")
@printf("  %-6s %-4s %-12s %-10s %-10s %-10s\n","bus","h","E[a]","σ[a]","λ*","util(λg)")
rows=[]
for h in harm, n in nonref
    bs="$n"; haskey(hd["nw"]["$h"]["bus"],bs) || continue
    bus=hd["nw"]["$h"]["bus"][bs]; haskey(bus,"ihdmax")||continue
    a=A[(h,n)]; Ea=mean(a); σa=std(a); σa>0||continue
    lim=(bus["ihdmax"]*v1[n])^2
    util=(Ea+λg*σa)*g[h]/lim
    push!(rows,(n,h,Ea,σa,λstar(a),util))
end
for r in sort(rows;by=x->-x[6])[1:min(6,length(rows))]
    @printf("  %-6d %-4d %-12.4e %-10.3e %-10.4f %-10.4f\n", r...)
end

# THD per bus: S=Σ_h a_h g_h, λ* at the Gaussian g
println("\nTHD (per bus, cross-harmonic at Gaussian g) — λ* = (q95−E)/σ of Σ_h|V|²:")
@printf("  %-6s %-12s %-10s %-10s %-10s\n","bus","E[S]","σ[S]","λ*","util(λg)")
trows=[]
for (bs,fb) in fund["bus"]
    n=fb["bus_i"]; haskey(bp,n)||continue; haskey(fb,"thdmax")||continue
    S=zeros(M); for h in harm; S .+= A[(h,n)].*g[h]; end
    Es=mean(S);σs=std(S);σs>0||continue
    lim=(fb["thdmax"]*v1[n])^2
    push!(trows,(n,Es,σs,λstar(S),(Es+λg*σs)/lim))
end
for r in sort(trows;by=x->-x[5])[1:min(5,length(trows))]
    @printf("  %-6d %-12.4e %-10.3e %-10.4f %-10.4f\n", r...)
end
println("\n  λ* ≪ Gaussian(1.645) ≪ Cantelli(4.36) ⇒ bounded/left-skewed output;")
println("  using λ* makes E[S]+λ*σ = q95(S) exactly ⇒ hosting = q05 (5% violation).")
println("="^74)
