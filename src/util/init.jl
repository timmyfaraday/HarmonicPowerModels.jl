################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Hakan Ergun                                                         #
################################################################################
# Changelog:                                                                   #
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
    hpf_results = solve_hpf(hpf_data, model_type, optimizer)

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["vm"] = hpf_results["solution"]["nw"]["1"]["bus"][i]["vm"]
        bus["va"] = hpf_results["solution"]["nw"]["1"]["bus"][i]["va"]
    end
    for (b, branch) in hdata["nw"]["1"]["branch"]
        branch["cm_fr"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cr_fr"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["ci_fr"]^2)
        branch["cm_to"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cr_to"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["ci_to"]^2)
    end
    for (x, xfmr) in hdata["nw"]["1"]["xfmr"]
        xfmr["ctm_fr"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["crt_fr"]^2 
                                + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cit_fr"]^2)
        xfmr["ctm_to"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["crt_to"]^2 
                                + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cit_to"]^2)
    end
end