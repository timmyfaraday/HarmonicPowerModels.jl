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
mutable struct NLP_DHHC <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
Deterministic Harmonic Hosting Capacity (QC)
"""
mutable struct QC_DHHC <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
Deterministic Harmonic Hosting Capacity (SOC)
"""
mutable struct SOC_DHHC <: _PMs.AbstractIVRModel _PMs.@pm_fields end