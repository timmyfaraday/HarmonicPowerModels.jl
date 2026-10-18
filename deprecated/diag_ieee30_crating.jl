# diag_ieee30_crating.jl — does removing the /√3 in c_rating fix the QCQP? (throwaway)
using HarmonicPowerModels, PowerModels
using Ipopt
import JuMP
using Printf
const PMs = PowerModels
const HPM = HarmonicPowerModels
include(joinpath(HPM.BASE_DIR, "examples", "shhc", "pce_params.jl"))
const SILENT = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")

datapath = joinpath(HPM.BASE_DIR, "examples", "shhc", "data", "30bus_network_shhc.m")
H = [1, 5, 7, 11, 13]; harm = [5,7,11,13]
data = PMs.parse_file(datapath)

# replicate sets c_rating = rate_a / sqrt(3) / vmin.  Override to rate_a / vmin
# (drop the questionable √3) and check feasibility + hosting capacity.
function solve_with_crating(drop_sqrt3::Bool)
    hd = HPM.replicate(data, H=H)
    if drop_sqrt3
        for (_, ntw) in hd["nw"], (_, br) in ntw["branch"]
            haskey(br, "c_rating") && (br["c_rating"] *= sqrt(3))   # undo the /√3
        end
    end
    try
        sol = HPM.solve_shhc(hd, sHHC_QCQP, SILENT; pce_params=PCE_PARAMS)
        return string(sol["termination_status"]), sol
    catch e
        return "ERROR: " * first(split(sprint(showerror,e),'\n')), nothing
    end
end

st0, _   = solve_with_crating(false)
st1, sol = solve_with_crating(true)
println("c_rating = rate_a/√3/vmin (as-is) : ", st0)
println("c_rating = rate_a/vmin   (no √3)  : ", st1)
if sol !== nothing && haskey(sol["solution"], "nw")
    @printf("\n  hosting capacity with corrected rating:\n")
    for h in harm
        @printf("    h=%-2d  I = %.6e pu\n", h, sol["solution"]["nw"]["$h"]["load"]["1"]["cmd"])
    end
end
println("\nDone.")
