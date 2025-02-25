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

# voltage ######################################################################
""
const voltage_thd_limits = Dict(
    "Clean Bus" =>                  0.00000,
    "IEC61000-2-4:2002, Cl. 2" =>   0.08000,
    "IEC61000-3-6:2008" =>          0.08000,
    "AS/NZS61000-3-6" =>            0.08000,
    "IEEE519-2022-1/69kV" =>        0.05000,
    "IEEE519-2022-69/161kV" =>      0.02500) 