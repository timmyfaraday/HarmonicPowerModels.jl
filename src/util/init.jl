################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Hakan Ergun                                                         #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

""
function update_hdata_with_fairness_principle_data!(hdata, model_type::Type, optimizer)
    if hdata["principle"] == "Kalai-Smorodinsky bargaining"
        for (l,load) in hdata["nw"]["1"]["load"]
            # make a deepcopy of the hdata
            hdata_temp = deepcopy(hdata)

            # set the principle to maximum efficiency
            hdata_temp["principle"] = "maximum efficiency"

            # remove all the loads except for l, for all harmonics
            for (nw,ntw) in hdata_temp["nw"], (nl,load) in ntw["load"] if nl ≠ l
                delete!(ntw["load"], nl)
            end end

            # solve the harmonic hosting capacity problem for the single load
            results_hhc = solve_hhc(hdata_temp, model_type, optimizer)

            # write away the solution for each network
            for (nw,ntw) in results_hhc["solution"]["nw"] if nw ≠ "1"
                hdata["nw"]["$nw"]["load"]["$l"]["cmdmax"] = ntw["load"]["$l"]["cmd"]
            end end
        end
    end
end

""
function update_hdata_with_fundamental_hpf_results!(hdata, model_type::Type, optimizer)
    # remove all but the fundamental network
    hpf_data = deepcopy(hdata)
    for n in keys(hpf_data["nw"])
        if n ≠ "1"
            delete!(hpf_data, n)
        end
    end

    # solve hpf problem for the fundamental harmonic only
    hpf_results = solve_hpf(hpf_data, model_type, optimizer)

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["vm"] = floor(hpf_results["solution"]["nw"]["1"]["bus"][i]["vm"], digits=10)
        bus["va"] = hpf_results["solution"]["nw"]["1"]["bus"][i]["va"]
    end
    for (b, branch) in hdata["nw"]["1"]["branch"]
        branch["cm_fr"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cr_fr"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["ci_fr"]^2)
        branch["cm_to"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cr_to"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["ci_to"]^2)
    end
    for (x, xfmr) in hdata["nw"]["1"]["xfmr"]
        xfmr["ctm_fr"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["crx_fr"]^2 
                                + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cix_fr"]^2)
        xfmr["ctm_to"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["crx_to"]^2 
                                + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cix_to"]^2)
    end
end