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

function create_pf_data_model(data)
    pf_data = deepcopy(data)

    for (n, network) in pf_data["nw"]
        if n ≠ "1"
            delete!(pf_data, n)
        end
    end

    return pf_data
end

function write_pf_results!(data, pf_result)
    for (b, bus) in data["nw"]["1"]["bus"]
        bus["vm"] = pf_result["solution"]["nw"]["1"]["bus"][b]["vm"]
        bus["va"] = pf_result["solution"]["nw"]["1"]["bus"][b]["va"]
    end

    for (br, branch) in data["nw"]["1"]["branch"]
        branch["cm_fr"] = sqrt(pf_result["solution"]["nw"]["1"]["branch"][br]["cr_fr"]^2 + pf_result["solution"]["nw"]["1"]["branch"][br]["ci_fr"]^2)
        branch["cm_to"] = sqrt(pf_result["solution"]["nw"]["1"]["branch"][br]["cr_to"]^2 + pf_result["solution"]["nw"]["1"]["branch"][br]["ci_to"]^2)
    end

    for (t, xfmr) in data["nw"]["1"]["xfmr"]
        xfmr["ctm_fr"] = sqrt(pf_result["solution"]["nw"]["1"]["xfmr"][t]["crt_fr"]^2 + pf_result["solution"]["nw"]["1"]["xfmr"][t]["cit_fr"]^2)
        xfmr["ctm_to"] = sqrt(pf_result["solution"]["nw"]["1"]["xfmr"][t]["crt_to"]^2 + pf_result["solution"]["nw"]["1"]["xfmr"][t]["cit_to"]^2)
    end
end
