# diag_qcqp_thd.jl — measure per-bus THD voltage & limit in dHHC vs QCQP (det E[φ])
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf
const PMs = PowerModels; const HPM = HarmonicPowerModels
include(joinpath(HPM.BASE_DIR, "examples", "shhc", "pce_params.jl"))
ipopt = optimizer_with_attributes(Ipopt.Optimizer,"linear_solver"=>"ma27","hsllib"=>HSL_jll.libhsl_path,"print_level"=>0)
mosek = optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true)
data = PMs.parse_file(joinpath(HPM.BASE_DIR,"examples","shhc","data","industrial_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]
α̂,β̂=PCE_PARAMS["alpha_hat"],PCE_PARAMS["beta_hat"]; Eω=α̂/(α̂+β̂)
Eφ(h)=(hb=PCE_PARAMS["harmonics"][h]; hb.v_min+(hb.v_max-hb.v_min)*Eω)
detp()=(p=deepcopy(PCE_PARAMS); p["harmonics"]=Dict{Int,Any}(h=>(Eφ(h),Eφ(h)) for h in harm); p)

hd=HPM.replicate(deepcopy(data),H=H); for h in harm,(_,l) in hd["nw"]["$h"]["load"]; l["ref_angle"]=rad2deg(Eφ(h)); end
sd = HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt)
sq = HPM.solve_shhc(HPM.replicate(deepcopy(data),H=H), sHHC_QCQP, ipopt; pce_params=detp())

# dHHC bus voltage at (nw,bus)
dV2(sol,h,i) = (b=sol["solution"]["nw"]["$h"]["bus"]["$i"]; b["vr"]^2 + b["vi"]^2)
# QCQP lifted |V|² k=0 if present, else from vr_pce/vi_pce[0]
function qV2(sol,h,i)
    b=sol["solution"]["nw"]["$h"]["bus"]["$i"]
    if haskey(b,"w_pce"); w=b["w_pce"]; return w[first(sort(collect(keys(w))))]
    elseif haskey(b,"vr_pce"); k0=first(sort(collect(keys(b["vr_pce"])))); return b["vr_pce"][k0]^2+b["vi_pce"][k0]^2
    else; return NaN; end
end
buses = sort(collect(keys(sd["solution"]["nw"]["1"]["bus"])), by=x->parse(Int,x))
fund = HPM.replicate(deepcopy(data),H=H)  # for thdmax (raw); use solution vm
println("per-bus THD²=Σ_h|V_h|²  vs limit=(thdmax·V1)²   [dHHC | QCQP]")
@printf("%-5s %-7s %-13s %-13s %-13s %-8s %-8s\n","bus","thdmax","V1_d","THD2_d","THD2_q","bind_d","bind_q")
for i in buses
    bd = sd["solution"]["nw"]["1"]["bus"][i]
    v1d = sqrt(get(bd,"vr",NaN)^2+get(bd,"vi",NaN)^2)
    thd = try PMs.parse_file(joinpath(HPM.BASE_DIR,"examples","shhc","data","industrial_network_shhc.m")); catch; end
    # thdmax from data bus_harmonics standard — read from hd fundamental bus
    tm = get(hd["nw"]["1"]["bus"][i], "thdmax", NaN)
    v1 = get(hd["nw"]["1"]["bus"][i], "vm", NaN)
    thd2d = sum(dV2(sd,h,i) for h in harm)
    thd2q = sum(qV2(sq,h,i) for h in harm)
    lim = (tm*v1)^2
    @printf("%-5s %-7.4f %-13.6e %-13.6e %-13.6e %-8s %-8s\n", i, tm, v1, thd2d, thd2q,
            thd2d/lim>0.999 ? "YES" : "no", thd2q/lim>0.999 ? "YES" : "no")
end
