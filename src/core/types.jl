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
# v0.3.0 - redefine structs                                                    #
################################################################################

""
abstract type AbstractHarmonicModel <: _PMs.AbstractIVRModel end
abstract type AbstractHHCModel      <: AbstractHarmonicModel end

""
mutable struct HarmonicPowerModel   <: AbstractHarmonicModel    _PMs.@pm_fields end
mutable struct dHHCPowerModel       <: AbstractHHCModel         _PMs.@pm_fields end