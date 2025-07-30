################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

""
function _IMs.solution_preprocessor(pm::AbstractHarmonicModel, solution::Dict)
    per_unit = _IMs.get_data(x -> x["per_unit"], pm.data, _PMs.pm_it_name; apply_to_subnetworks = false)
    solution["it"][_PMs.pm_it_name]["per_unit"]     = per_unit
    solution["it"][_PMs.pm_it_name]["s_base_mva"]   = pm.data["s_base_mva"]
end

""
function sol_data_model!(pm::AbstractHarmonicModel, solution::Dict)
    _PMs.apply_pm!(_sol_data_model_ivr!, solution)
end

""
function _sol_data_model_ivr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vbr") && haskey(bus, "vbi")
                bus["vbm"] = hypot(bus["vbr"], bus["vbi"])
                bus["vba"] = atan(bus["vbi"], bus["vbr"]) * 180 / pi
            end
        end
    end

    if haskey(solution, "branch")
        for (i, branch) in solution["branch"]
            if haskey(branch, "pb_fr") && haskey(branch, "pb_to")
                branch["pb_loss"] = branch["pb_fr"] + branch["pb_to"]
            end
            if haskey(branch, "qb_fr") && haskey(branch, "qb_to")
                branch["qb_loss"] = branch["qb_fr"] + branch["qb_to"]
            end
        end
    end

    if haskey(solution, "xfmr")
        for (x, xfmr) in solution["xfmr"]
            # power loss
            xfmr["px_loss"] = sum(xfmr[nk] for nk in keys(xfmr) if startswith(nk, "px_"); init=0.0)
            xfmr["qx_loss"] = sum(xfmr[nk] for nk in keys(xfmr) if startswith(nk, "qx_"); init=0.0)
            # excitation voltage
            if haskey(xfmr, "exr") && haskey(xfmr, "exi")
                xfmr["exm"] = hypot(xfmr["exr"], xfmr["exi"])
                xfmr["exa"] = atan(xfmr["exi"], xfmr["exr"]) * 180 / pi
            end
            # winding voltage
            for nk in keys(xfmr) if startswith(nk, "vxr_")
                nb = nk[5:end]
                xfmr["vxm_$nb"] = hypot(xfmr["vxr_$nb"], xfmr["vxi_$nb"])
                xfmr["vxa_$nb"] = atan(xfmr["vxi_$nb"], xfmr["vxr_$nb"]) * 180 / pi
    end end end end

    if haskey(solution, "gen")
        for (i, gen) in solution["gen"]
            if haskey(gen, "cgr") && haskey(gen, "cgi")
                gen["cgm"] = hypot(gen["cgr"], gen["cgi"])
    end end end
    
    if haskey(solution, "load")
        for (i, load) in solution["load"]
            if haskey(load, "clr") && haskey(load, "cli")
                load["clm"] = hypot(load["clr"], load["cli"])
end end end end