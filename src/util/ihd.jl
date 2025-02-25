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
const voltage_ihd_limits = Dict(
    "Clean Bus" =>                 [Inf, zeros(49)...],
    "IEC61000-2-4:2002, Cl. 2" =>  [Inf, 0.02000, 0.05000, 0.01000, 0.06000, 
                                    0.00500, 0.05000, 0.00500, 0.01500, 0.00500, 
                                    0.03500, 0.00458, 0.03000, 0.00429, 0.00400, 
                                    0.00406, 0.02000, 0.00389, 0.01761, 0.00375, 
                                    0.00200, 0.00364, 0.01408, 0.00354, 0.01275, 
                                    0.00346, 0.00200, 0.00339, 0.01061, 0.00333, 
                                    0.00975, 0.00328, 0.00200, 0.00324, 0.00833, 
                                    0.00319, 0.00773, 0.00316, 0.00200, 0.00313, 
                                    0.00671, 0.00310, 0.00627, 0.00307, 0.00200, 
                                    0.00304, 0.00551, 0.00302, 0.00518, 0.00300],
    "IEC61000-3-6:2008" =>         [Inf, 0.01400, 0.02000, 0.00800, 0.02000,
                                    0.00400, 0.02000, 0.00400, 0.01000, 0.00350,
                                    0.01500, 0.00318, 0.01500, 0.00296, 0.00300,
                                    0.00279, 0.01200, 0.00266, 0.01074, 0.00255,
                                    0.00200, 0.00246, 0.00887, 0.00239, 0.00816,
                                    0.00233, 0.00200, 0.00228, 0.00703, 0.00223,
                                    0.00658, 0.00219, 0.00200, 0.00216, 0.00583,
                                    0.00213, 0.00551, 0.00210, 0.00200, 0.00208,
                                    0.00498, 0.00205, 0.00474, 0.00203, 0.00200,
                                    0.00201, 0.00434, 0.00200, 0.00416, 0.00198],
    "AS/NZS61000-3-6" =>           [Inf, 0.02000, 0.05000, 0.01000, 0.05000,
                                    0.00500, 0.05000, 0.00500, 0.01500, 0.00500,
                                    0.03500, 0.00200, 0.03000, 0.00200, 0.00300,
                                    0.00200, 0.02000, 0.00200, 0.01500, 0.00200,
                                    0.00200, 0.00200, 0.01500, 0.00200, 0.01500,
                                    0.00200, 0.00200, 0.00200, 0.01483, 0.00200,
                                    0.01087, 0.00200, 0.00200, 0.00200, 0.00986,
                                    0.00200, 0.00943, 0.00200, 0.00200, 0.00200,
                                    0.00871, 0.00200, 0.00840, 0.00200, 0.00200,
                                    0.00200, 0.00785, 0.00200, 0.00761, 0.00200],
    "IEEE519-2022-1/69kV" =>       [Inf, 0.030 .* ones(49)...],
    "IEEE519-2022-69/161kV" =>     [Inf, 0.015 .* ones(49)...])

# current ######################################################################
""
function current_ihd_limits_ieee519(voltage, i_ratio, harmonic)
    if voltage < 69.0
        if i_ratio == Inf
            return 0.0
        elseif i_ratio < 20
            if 2 <= harmonic < 11
                return 0.04
            elseif 11 <= harmonic < 17
                return 0.02
            elseif 17 <= harmonic < 23
                return 0.015
            elseif 23 <= harmonic < 35
                return 0.006
            elseif 35 <= harmonic <= 50
                return 0.003
            else
                return 0.0
            end
        elseif 20.0 <= i_ratio < 50.0
            if 2 <= harmonic < 11
                return 0.07
            elseif 11 <= harmonic < 17
                return 0.035
            elseif 17 <= harmonic < 23
                return 0.025
            elseif 23 <= harmonic < 35
                return 0.010
            elseif 35 <= harmonic <= 50
                return 0.005
            else
                return 0.0
            end
        elseif 50.0 <= i_ratio < 100.0
            if 2 <= harmonic < 11
                return 0.10
            elseif 11 <= harmonic < 17
                return 0.045
            elseif 17 <= harmonic < 23
                return 0.040
            elseif 23 <= harmonic < 35
                return 0.015
            elseif 35 <= harmonic <= 50
                return 0.007
            else
                return 0.0
            end
        elseif 100.0 <= i_ratio < 1000.0
            if 2 <= harmonic < 11
                return 0.12
            elseif 11 <= harmonic < 17
                return 0.055
            elseif 17 <= harmonic < 23
                return 0.050
            elseif 23 <= harmonic < 35
                return 0.020
            elseif 35 <= harmonic <= 50
                return 0.010
            else
                return 0.0
            end
        else
            if 2 <= harmonic < 11
                return 0.15
            elseif 11 <= harmonic < 17
                return 0.07
            elseif 17 <= harmonic < 23
                return 0.06
            elseif 23 <= harmonic < 35
                return 0.025
            elseif 35 <= harmonic <= 50
                return 0.014
            else
                return 0.0
            end
        end
    else
        if i_ratio == Inf
            return 0.0
        elseif i_ratio < 20       
            if 2 <= harmonic < 11
                return 0.02
            elseif 11 <= harmonic < 17
                return 0.01
            elseif 17 <= harmonic < 23
                return 0.0075
            elseif 23 <= harmonic < 35
                return 0.003
            elseif 35 <= harmonic <= 50
                return 0.0015
            else
                return 0.0
            end
        elseif 20.0 <= i_ratio < 50.0
            if 2 <= harmonic < 11
                return 0.035
            elseif 11 <= harmonic < 17
                return 0.0175
            elseif 17 <= harmonic < 23
                return 0.0125
            elseif 23 <= harmonic < 35
                return 0.005
            elseif 35 <= harmonic <= 50
                return 0.0025
            else
                return 0.0
            end
        elseif 50.0 <= i_ratio < 100.0
            if 2 <= harmonic < 11
                return 0.05
            elseif 11 <= harmonic < 17
                return 0.0225
            elseif 17 <= harmonic < 23
                return 0.02
            elseif 23 <= harmonic < 35
                return 0.0075
            elseif 35 <= harmonic <= 50
                return 0.0035
            else
                return 0.0
            end
        elseif 100.0 <= i_ratio < 1000.0
            if 2 <= harmonic < 11
                return 0.06
            elseif 11 <= harmonic < 17
                return 0.0275
            elseif 17 <= harmonic < 23
                return 0.025
            elseif 23 <= harmonic < 35
                return 0.01
            elseif 35 <= harmonic <= 50
                return 0.005
            else
                return 0.0
            end
        else
            if 2 <= harmonic < 11
                return 0.075
            elseif 11 <= harmonic < 17
                return 0.035
            elseif 17 <= harmonic < 23
                return 0.03
            elseif 23 <= harmonic < 35
                return 0.0125
            elseif 35 <= harmonic <= 50
                return 0.007
            else
                return 0.0
end end end end
