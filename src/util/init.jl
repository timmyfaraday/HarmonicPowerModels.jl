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

    # solve hpf problem for the fundamental harmonic only
    hpf_results = solve_hopf(hpf_data, HarmonicPowerModel, optimizer)

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["v_fund_magn"] = sqrt(hpf_results["solution"]["nw"]["1"]["bus"][i]["vbr"]^2 + hpf_results["solution"]["nw"]["1"]["bus"][i]["vbi"]^2)
        bus["va"] = atan(hpf_results["solution"]["nw"]["1"]["bus"][i]["vbi"] / hpf_results["solution"]["nw"]["1"]["bus"][i]["vbr"])
    end
    for (b, branch) in hdata["nw"]["1"]["branch"]
        branch["cm_fr"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cbr_fr"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["cbi_fr"]^2)
        branch["cm_to"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cbr_to"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["cbi_to"]^2)
    end
    for (x, xfmr) in hdata["nw"]["1"]["xfmr"]
        w1_idx = xfmr["bus"][1]
        w2_idx = xfmr["bus"][2] 
        xfmr["ctm_fr"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cxr_"*"$w1_idx"]^2 
                                + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cxi_"*"$w1_idx"]^2)
        xfmr["ctm_to"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cxr_"*"$w2_idx"]^2 
                                + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cxi_"*"$w2_idx"]^2)
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