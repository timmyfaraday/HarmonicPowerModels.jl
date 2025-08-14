################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Hakan Ergun                                                         #
################################################################################
# Changelog:                                                                   #
#                                                                              #
################################################################################
function update_hdata_with_fundamental_hpf_results!(hdata, model_type::Type, optimizer)
    # remove all but the fundamental network
    hpf_data = deepcopy(hdata)
    for n in keys(hpf_data["nw"])
        if n ≠ "1"
            delete!(hpf_data, n)
        end
    end

    # solve hpf problem for the fundamental harmonic
    fpf = solve_hpf(hpf_data, HarmonicPowerModel, optimizer)["solution"]["nw"]["1"]

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["v_fund_magn"] = fpf["bus"][i]["vbm"]
    end
    for (b, branch) in hdata["nw"]["1"]["branch"]
        branch["i_fund_magn"] = [fpf["branch"][b]["cbm_fr"], fpf["branch"][b]["cbm_to"]]
    end
    for (x, xfmr) in hdata["nw"]["1"]["xfmr"]
        xfmr["i_fund_magn"] = [fpf["xfmr"][x]["cxm_$nb"] for nb in xfmr["bus"]]
    end
    for (g, gen) in hdata["nw"]["1"]["gen"]
        gen["i_fund_magn"] = fpf["gen"][g]["cgm"]
    end
end

function calculate_maximum_harmonic_source_current_injection!(data, Zh)    
    for (nw, network) in data["nw"]
      if network !== "nw"
          h = parse(Int, nw)
          for (r, hsrc) in network["hsrc"]
              bus = data["hsrc"]["nw"]["1"]["bus"]
              v_ihd_max = network["bus"]["$bus"]["v_ihd_max"] 
              hb_idx = network["bus"]["$bus"]["index"]
              load["cmdmax"] = v_ihd_max / sqrt(Zh["$hb_idx"][h]["re"]^2 + Zh["$hb_idx"][h]["im"]^2)
          end
      end
    end
end