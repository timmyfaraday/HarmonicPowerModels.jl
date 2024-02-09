################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for power harmonics     #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

""
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int) = haskey(aim.ref[:it][it][:nw], nw)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol) = haskey(aim.ref[:it][it][:nw][nw],key)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx) = haskey(aim.ref[:it][it][:nw][nw][key],idx)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx, param::String) = haskey(aim.ref[:it][it][:nw][nw][key][idx],param)