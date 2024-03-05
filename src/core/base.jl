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
################################################################################

""
is_pos_sequence(nh::Int) = nh % 3 == 1
is_neg_sequence(nh::Int) = nh % 3 == 2
is_zero_sequence(nh::Int) = nh % 3 == 0

""
ids(pm::_PMs.AbstractPowerModel, key::Symbol) = _PMs.ids(pm, key, nw=fundamental(pm))

""
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int) = haskey(aim.ref[:it][it][:nw], nw)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol) = haskey(aim.ref[:it][it][:nw][nw],key)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx) = haskey(aim.ref[:it][it][:nw][nw][key],idx)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx, param::String) = haskey(aim.ref[:it][it][:nw][nw][key][idx],param)