# diag_ieee30_qcqp.jl — root-cause the IEEE-30 sHHC_QCQP infeasibility (throwaway)
using HarmonicPowerModels, PowerModels
using Ipopt
import JuMP
using Printf
const PMs = PowerModels
const HPM = HarmonicPowerModels
include(joinpath(HPM.BASE_DIR, "examples", "shhc", "pce_params.jl"))
const SILENT = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")

data = PMs.parse_file(joinpath(HPM.BASE_DIR, "examples", "shhc", "data", "30bus_network_shhc.m"))
H = [1, 5, 7, 11, 13]

# bus types in raw data
type3 = sort([parse(Int,k) for (k,b) in data["bus"] if b["bus_type"]==3])
println("raw data: bus_type==3 buses = ", type3, "  (", length(type3), " of ", length(data["bus"]), ")")

hdata = HPM.replicate(data, H=H)

# ref_buses per network (PowerModels ref) — need a pm; instead inspect bus_type per nw
for nw in ["1","5"]
    t3 = sort([parse(Int,k) for (k,b) in hdata["nw"][nw]["bus"] if b["bus_type"]==3])
    println("nw=$nw: bus_type==3 = ", t3)
end

# 1. fundamental HPF
println("\n── fundamental HPF ──")
hpf_status = "?"
try
    HPM.update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, SILENT)
    global hpf_status = "ok"
catch e
    global hpf_status = "ERROR: " * first(split(sprint(showerror,e),'\n'))
end
println("HPF status: ", hpf_status)

# 2. chance-constraint RHS sign audit (fundamental results now in hdata["nw"]["1"])
function audit_v(fund)
    println("\n── voltage RMS limits  vmaxrms² − v_fund²  (negative ⇒ infeasible) ──")
    neg = 0
    for (i,bus) in sort(collect(fund["bus"]), by=x->parse(Int,x[1]))
        vm = get(bus,"vm",NaN); vmaxrms = get(bus,"vmaxrms", nothing)
        vmaxrms === nothing && continue
        lim = vmaxrms^2 - vm^2
        if lim < 0
            neg += 1
            @printf("  bus %-3s vm=%.4f vmaxrms=%.4f  lim=%.4e  *** NEG ***\n", i, vm, vmaxrms, lim)
        end
    end
    println("  negative V-RMS limits: ", neg)
end

function audit_i(fund)
    println("\n── current RMS limits  c_rating² − cm_fr²  (negative ⇒ infeasible) ──")
    neg = 0
    for (b,br) in sort(collect(fund["branch"]), by=x->parse(Int,x[1]))
        cr = get(br,"c_rating", nothing); cmfr = get(br,"cm_fr", nothing)
        (cr === nothing || cmfr === nothing) && continue
        lim = cr^2 - cmfr^2
        if lim < 0
            neg += 1
            @printf("  branch %-3s cm_fr=%.4f c_rating=%.4f  lim=%.4e  *** NEG ***\n", b, cmfr, cr, lim)
        end
    end
    println("  negative I-RMS limits: ", neg, " of ", length(fund["branch"]), " branches")
end

audit_v(hdata["nw"]["1"])
audit_i(hdata["nw"]["1"])

# 3. bisect: disable chance constraints one family at a time via typed-signature override
println("\n── bisect: which chance-constraint family causes infeasibility? ──")
noop_ihd() = @eval HPM constraint_voltage_ihd_chance(pm::sHHC_QCQP, i::Int; nw::Int=fundamental(pm)) = nothing
noop_thd() = @eval HPM constraint_voltage_thd_chance(pm::sHHC_QCQP, i::Int; harmonics, K=0:2) = nothing
noop_vrms() = @eval HPM constraint_voltage_rms_chance(pm::sHHC_QCQP, i::Int; nw_ids, K=0:2) = nothing
noop_irms() = @eval HPM constraint_current_rms_chance(pm::sHHC_QCQP, b::Int; nw_ids, K=0:2) = nothing

function solve_status()
    hd = HPM.replicate(data, H=H)
    try
        sol = HPM.solve_shhc(hd, sHHC_QCQP, SILENT; pce_params=PCE_PARAMS)
        return string(sol["termination_status"])
    catch e
        return "ERROR: " * first(split(sprint(showerror,e),'\n'))
    end
end

println("  baseline (all on)         : ", solve_status())
noop_irms();                 println("  − current RMS             : ", solve_status())
noop_vrms();                 println("  − current & voltage RMS   : ", solve_status())
noop_thd();                  println("  − RMS & THD               : ", solve_status())
noop_ihd();                  println("  − RMS & THD & IHD (none)  : ", solve_status())
println("\nDiagnostic complete.")
