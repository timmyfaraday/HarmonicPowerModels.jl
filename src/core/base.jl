################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
################################################################################

""
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int) = haskey(aim.ref[:it][it][:nw], nw)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol) = haskey(aim.ref[:it][it][:nw][nw],key)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx) = haskey(aim.ref[:it][it][:nw][nw][key],idx)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx, param::String) = haskey(aim.ref[:it][it][:nw][nw][key][idx],param)