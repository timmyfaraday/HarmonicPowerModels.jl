################################################################################
# Diagnostic: locate every reference (dHHC@E[φ], dHHC worst-case, sHHC Cantelli/  #
# Gaussian) inside the MC hosting distribution, and report the IMPLIED VIOLATION  #
# probability of committing each hosting level:                                   #
#     viol(h) = P[H(ω) < h] = fraction of realisations that cannot support h.     #
# For a correct (1−ε) chance constraint the sHHC should give viol ≈ ε (=5%),      #
# i.e. sit at the 5th percentile of HOSTING (= 95th percentile of DISTORTION).    #
# Industrial network (clean: dHHC↔SOC gap ~0.25%).                                #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Distributions, Random, Statistics
const HPM = HarmonicPowerModels
const PMs = PowerModels
include(joinpath(@__DIR__, "pce_params.jl"))

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "linear_solver"=>"ma27",
            "hsllib"=>HSL_jll.libhsl_path, "print_level"=>0,
            "tol"=>1e-10, "constr_viol_tol"=>1e-10)
mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET"=>true)

data = PMs.parse_file(joinpath(@__DIR__, "data", "industrial_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]
α̂,β̂ = PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]; Eω=α̂/(α̂+β̂); ε=PCE_PARAMS["epsilon"]
base = Dict(h=>PCE_PARAMS["harmonics"][h] for h in harm)
units = sort(parse.(Int, collect(keys(data["load"])))); uset=Set(units)
get_I(s,h)=(x=s["solution"]; haskey(x,"nw") ? x["nw"]["$h"]["load"]["1"]["cmd"] : x[h])
N_MC = 4000

function dhhc(φ)   # φ: h=>unit=>rad
    hd=HPM.replicate(deepcopy(data),H=H)
    for h in harm, k in keys(hd["nw"]["$h"]["load"])
        d=parse(Int,k); d in uset && (hd["nw"]["$h"]["load"][k]["ref_angle"]=rad2deg(φ[h][d]))
    end
    HPM.solve_hhc(hd,dHHC_SOC,mosek,ipopt)["objective"]
end
function soc(bound)
    p=deepcopy(PCE_PARAMS); p["harmonics"]=Dict{Int,Any}(h=>base[h] for h in harm)
    p["harmonic_units"]=units; p["bound"]=bound
    HPM.solve_shhc(HPM.replicate(deepcopy(data),H=H),sHHC_SOC,mosek;pce_params=p,hpf_optimizer=ipopt)["objective_value"]
end

# references
Eφ   = Dict(h=>Dict(d=>base[h][1]+(base[h][2]-base[h][1])*Eω for d in units) for h in harm)   # mean angle (coherent)
H_Eφ = dhhc(Eφ)
Z,bp = HPM.reduced_harmonic_impedance(HPM.replicate(deepcopy(data),H=H))
ubus = unique(HPM.replicate(deepcopy(data),H=H)["nw"]["1"]["load"]["$d"]["load_bus"] for d in units)
H_wc = minimum(let ang=HPM.constructive_worst_case_angles(Z,bp,HPM.replicate(deepcopy(data),H=H),ns)
                   dhhc(Dict(h=>Dict(d=>ang[h][d] for d in units) for h in harm)) end
               for ns in ubus if haskey(bp,ns))
H_c, H_g = soc("cantelli"), soc("gaussian")

# MC
Random.seed!(1); B=Distributions.Beta(α̂,β̂); obj=zeros(N_MC)
for m in 1:N_MC
    ω=Dict(d=>rand(B) for d in units)
    φ=Dict(h=>Dict(d=>base[h][1]+(base[h][2]-base[h][1])*ω[d] for d in units) for h in harm)
    obj[m]=dhhc(φ)
end
sort!(obj)
pct(h) = 100*count(<(h),obj)/N_MC                       # implied violation % = P[H(ω)<h]
qz(p)  = quantile(obj,p)

@printf("\nIndustrial — MC dHHC hosting distribution (N=%d)\n", N_MC)
@printf("  min=%.6e  q05=%.6e  median=%.6e  q95=%.6e  max=%.6e  mean=%.6e\n",
        minimum(obj), qz(.05), qz(.5), qz(.95), maximum(obj), mean(obj))
println("\n  reference            hosting        MC-percentile   implied violation P[H<h]")
for (nm,h) in (("dHHC worst-case",H_wc),("dHHC @E[φ] (mean)",H_Eφ),
               ("sHHC Cantelli",H_c),("sHHC Gaussian",H_g),
               ("MC q05",qz(.05)),("MC median",qz(.5)),("MC q95",qz(.95)))
    @printf("  %-20s %.6e   %6.1f%%         %6.1f%%\n", nm, h, pct(h), pct(h))
end
println("\n  A valid 95% chance constraint ⇒ implied violation ≈ 5% ⇒ sits at hosting-q05")
println("  (= distortion-q95). 'Near q95 of hosting' would mean ~95% violation.")
