################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
################################################################################

"""
Deterministic Harmonic Hosting Capacity (NLP)
"""
mutable struct dHHC_NLP <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
Deterministic Harmonic Hosting Capacity (SOC)
"""
mutable struct dHHC_SOC <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
Abstract Stochastic Harmonic Hosting Capacity model
"""
abstract type AbstractSHHCModel <: _PMs.AbstractIVRModel end

"""
Stochastic Harmonic Hosting Capacity (SOC, PCE)
"""
mutable struct sHHC_SOC <: AbstractSHHCModel _PMs.@pm_fields end