################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

"""
Deterministic Harmonic Hosting Capacity (NLP)
"""
mutable struct dHHC_NLP <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
Deterministic Harmonic Hosting Capacity (SOC)
"""
mutable struct dHHC_SOC <: _PMs.AbstractIVRModel _PMs.@pm_fields end