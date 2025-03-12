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