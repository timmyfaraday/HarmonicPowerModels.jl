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
# v0.3.0 - reviewed TVA                                                        #
################################################################################

""
is_pos_sequence(nh::Int)  = nh % 3 == 1
is_neg_sequence(nh::Int)  = nh % 3 == 2
is_zero_sequence(nh::Int) = nh % 3 == 0

""
ids(pm::_PMs.AbstractPowerModel, key::Symbol) = _PMs.ids(pm, key, nw=fundamental(pm))

""
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int) = haskey(aim.ref[:it][it][:nw], nw)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol) = haskey(aim.ref[:it][it][:nw][nw],key)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx) = haskey(aim.ref[:it][it][:nw][nw][key],idx)
hasref(aim::_IMs.AbstractInfrastructureModel, it::Symbol, nw::Int, key::Symbol, idx, param::String) = haskey(aim.ref[:it][it][:nw][nw][key][idx],param)

""
function solve_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method;
                        ref_extensions=[], solution_processors=[], relax_integrality=false,
                        multinetwork=false, kwargs...)

    if multinetwork != _IMs.ismultinetwork(data)
        model_requirement   = multinetwork ? "multi-network" : "single-network"
        data_type           = _IMs.ismultinetwork(data) ? "multi-network" : "single-network"
        _MEM.error(_PMs._LOGGER, "attempted to build a $(model_requirement) model with $(data_type) data")
    end

    start_time = time()
    pm = _IMs.instantiate_model(data, model_type, build_method, _HPM.ref_add_core!, _PMs._pm_global_keys, _PMs.pm_it_sym; kwargs...)
    _MEM.debug(_PMs._LOGGER, "pm model build time: $(time() - start_time)")

    start_time = time()
    result = _PMs.optimize_model!(pm, relax_integrality=relax_integrality, optimizer=optimizer, solution_processors=solution_processors)
    _MEM.debug(_PMs._LOGGER, "pm model solve and solution time: $(time() - start_time)")

    return result
end


