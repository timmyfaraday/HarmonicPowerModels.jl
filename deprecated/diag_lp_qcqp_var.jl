# diag_lp_qcqp_var.jl — compare LP vs QCQP |V|²/I² PCE (mean & variance) per
# harmonic at a bus, uncertain case. Are the variance coefficients equal?
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf
const PMs=PowerModels; const HPM=HarmonicPowerModels
include(joinpath(HPM.BASE_DIR,"examples","shhc","pce_params.jl"))
ipopt=optimizer_with_attributes(Ipopt.Optimizer,"linear_solver"=>"ma27","hsllib"=>HSL_jll.libhsl_path,"print_level"=>0)
mosek=optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true)
data=PMs.parse_file(joinpath(HPM.BASE_DIR,"examples","shhc","data","industrial_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]

# LP (SOCP) — caches pce["alpha_V"][h] = |V|²/I² PCE (Nn × P), matrix bus index
hl=HPM.replicate(deepcopy(data),H=H)
sl=HPM.solve_shhc(hl, sHHC_LP, mosek; pce_params=deepcopy(PCE_PARAMS), hpf_optimizer=ipopt)
pce=hl["pce"]; psi=pce["psi_norms"]; aV=pce["alpha_V"]
fund=hl["nw"]["1"]; bus_ids=sort([b["bus_i"] for (_,b) in fund["bus"]]); bidx=Dict(id=>i for (i,id) in enumerate(bus_ids))

# QCQP — w_pce per bus (|V|²), divide by I_h² to normalise
hq=HPM.replicate(deepcopy(data),H=H)
sq=HPM.solve_shhc(hq, sHHC_QCQP, ipopt; pce_params=deepcopy(PCE_PARAMS))
getI(sol,h)=sol["solution"]["nw"]["$h"]["load"]["1"]["cmd"]

var_of(coef) = sum(psi[k+1]*coef[k+1]^2 for k in 1:length(coef)-1; init=0.0)

println("bus 8: |V|²/I² mean (k=0) and variance — LP vs QCQP, per harmonic")
@printf("%-3s | %-13s %-13s %-8s | %-13s %-13s %-8s\n","h","mean_LP","mean_QCQP","Δ%","var_LP","var_QCQP","Δ%")
n8 = bidx[8]
for h in harm
    aLP = aV[h][n8, :]                     # length-P
    Ih  = getI(sq,h)
    wq  = sq["solution"]["nw"]["$h"]["bus"]["8"]["w_pce"]   # Dict k=>|V|²
    ks  = sort(collect(keys(wq)))
    aQ  = [wq[k]/Ih^2 for k in ks]         # |V|²/I² PCE (normalised), length-P
    mL, mQ = aLP[1], aQ[1]
    vL, vQ = var_of(aLP), var_of(aQ)
    @printf("%-3d | %-13.6e %-13.6e %+-8.2f | %-13.6e %-13.6e %+-8.2f\n",
            h, mL, mQ, 100*(mL/mQ-1), vL, vQ, 100*(vL/vQ-1))
end
