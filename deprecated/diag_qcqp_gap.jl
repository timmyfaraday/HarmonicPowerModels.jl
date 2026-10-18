# diag_qcqp_gap.jl — why is the deterministic sHHC_QCQP more conservative than
# dHHC_SOC / sHHC_LP at the point-mass angle E[φ]? (throwaway; run from main env)
using HarmonicPowerModels, PowerModels
using Ipopt, HiGHS, MosekTools, HSL_jll
using JuMP, Printf
const PMs = PowerModels; const HPM = HarmonicPowerModels
include(joinpath(HPM.BASE_DIR, "examples", "shhc", "pce_params.jl"))

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "linear_solver"=>"ma27",
            "hsllib"=>HSL_jll.libhsl_path, "print_level"=>0)
mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET"=>true)
highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)

data = PMs.parse_file(joinpath(HPM.BASE_DIR,"examples","shhc","data","industrial_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]
α̂,β̂ = PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]; Eω=α̂/(α̂+β̂)
Eφ(h)=(hb=PCE_PARAMS["harmonics"][h]; hb.v_min+(hb.v_max-hb.v_min)*Eω)
detp() = (p=deepcopy(PCE_PARAMS); p["harmonics"]=Dict{Int,Any}(h=>(Eφ(h),Eφ(h)) for h in harm); p)
getI(sol,h)=(s=sol["solution"]; haskey(s,"nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h])

# baselines
hd=HPM.replicate(deepcopy(data),H=H); for h in harm,(_,l) in hd["nw"]["$h"]["load"]; l["ref_angle"]=rad2deg(Eφ(h)); end
sol_d = HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt)
sol_l = HPM.solve_shhc(HPM.replicate(deepcopy(data),H=H), sHHC_LP, highs; pce_params=detp(), hpf_optimizer=ipopt)

qcqp() = HPM.solve_shhc(HPM.replicate(deepcopy(data),H=H), sHHC_QCQP, ipopt; pce_params=detp())

println("baseline dHHC vs LP vs QCQP (deterministic E[φ]):")
sol_q0 = qcqp()
@printf("%-3s %-14s %-14s %-14s %-10s\n","h","dHHC","LP","QCQP","QCQP/dHHC%")
for h in harm
    @printf("%-3d %-14.6e %-14.6e %-14.6e %+.2f\n",h,getI(sol_d,h),getI(sol_l,h),getI(sol_q0,h),(getI(sol_q0,h)/getI(sol_d,h)-1)*100)
end

# cumulative disabling of QCQP chance-constraint families
noop_irms()=@eval HPM constraint_current_rms_chance(pm::sHHC_QCQP,b::Int;nw_ids,K=0:2)=nothing
noop_vrms()=@eval HPM constraint_voltage_rms_chance(pm::sHHC_QCQP,i::Int;nw_ids,K=0:2)=nothing
noop_thd() =@eval HPM constraint_voltage_thd_chance(pm::sHHC_QCQP,i::Int;harmonics,K=0:2)=nothing

println("\ncumulative disable → QCQP/dHHC % per harmonic:")
function row(label)
    s=qcqp()
    @printf("%-26s ",label)
    for h in harm; @printf("%+7.2f ",(getI(s,h)/getI(sol_d,h)-1)*100); end
    println()
end
row("all chance on")
noop_irms(); row("− current RMS")
noop_vrms(); row("− current+voltage RMS")
noop_thd();  row("− RMS+THD (IHD only)")
println("\n(≈0 means QCQP matches dHHC. harm order: ", harm, ")")
