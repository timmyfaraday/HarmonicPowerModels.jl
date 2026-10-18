# diag_ieee30_slack.jl — is the branch overload caused by the 6 slack buses? (throwaway)
using HarmonicPowerModels, PowerModels
using Ipopt
import JuMP
using Printf
const PMs = PowerModels
const HPM = HarmonicPowerModels
const SILENT = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")

datapath = joinpath(HPM.BASE_DIR, "examples", "shhc", "data", "30bus_network_shhc.m")
H = [1, 5, 7, 11, 13]
watch = ["18", "27"]

function run_case(label, flip_to_pv)
    data = PMs.parse_file(datapath)
    if flip_to_pv
        for (k,b) in data["bus"]
            if b["bus_type"] == 3 && parse(Int,k) != 1   # keep bus 1 as the only slack
                b["bus_type"] = 2
            end
        end
    end
    n3 = count(b->b[2]["bus_type"]==3, data["bus"])
    hdata = HPM.replicate(data, H=H)
    HPM.update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, SILENT)
    fund = hdata["nw"]["1"]
    neg = 0
    for (b,br) in fund["branch"]
        cr = get(br,"c_rating",nothing); cm = get(br,"cm_fr",nothing)
        (cr===nothing||cm===nothing) && continue
        cr^2 - cm^2 < 0 && (neg += 1)
    end
    println("\n[$label]  #slack=$n3   negative I-RMS branches: $neg")
    for b in watch
        br = fund["branch"][b]
        @printf("    branch %-3s: rate_a=%.4f  c_rating=%.4f  cm_fr=%.4f  headroom=%.4e %s\n",
                b, get(br,"rate_a",NaN), get(br,"c_rating",NaN), get(br,"cm_fr",NaN),
                get(br,"c_rating",0.0)^2-get(br,"cm_fr",0.0)^2,
                get(br,"c_rating",0.0)^2-get(br,"cm_fr",0.0)^2 < 0 ? "NEG" : "ok")
    end
end

run_case("as-is: 6 slack buses", false)
run_case("forced: single slack (bus 1), rest PV", true)
println("\nDone.")
