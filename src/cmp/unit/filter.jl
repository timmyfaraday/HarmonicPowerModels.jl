################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# parameters ###################################################################
""
function add_filter_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (f, filter) in ntw["filter"]
        h       = parse(Int, nw)
        idata   = fdata["filter"][f] 
        if nw == "1"
            filter["id"]        = idata["index"]
            filter["bus"]       = idata["bus"]
            filter["type"]      = idata["type"]
end end end

# variables ####################################################################
""
function variable_filter_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_filter_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_filter_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_filter_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cfr = _PMs.var(pm, nw)[:cfr] = 
            JuMP.@variable( pm.model,
                            [f in _PMs.ids(pm, nw, :filter)], 
                            base_name="$(nw)_cfr",
                            start=0.0)

    report && _PMs.sol_component_value(pm, nw, :filter, :cfr, _PMs.ids(pm, nw, :filter), cfr)
end
""
function variable_filter_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cfi = _PMs.var(pm, nw)[:cfi] = 
            JuMP.@variable( pm.model,
                            [f in _PMs.ids(pm, nw, :filter)], 
                            base_name="$(nw)_cfi",
                            start=0.0)

    report && _PMs.sol_component_value(pm, nw, :filter, :cfi, _PMs.ids(pm, nw, :filter), cfi)
end

# constraints ##################################################################
## active filter current constraint ############################################
""
function constraint_active_filter_current(pm::_PMs.AbstractPowerModel, f::Int)
    i       = _PMs.ref(pm, fundamental(pm), :filter, f, "bus")

    type    = _PMs.ref(pm, fundamental(pm), :filter, f, "type")

    if type == "a"
        constraint_active_filter_current(pm, f, i)
end end
""
function constraint_active_filter_current(pm::_PMs.AbstractIVRModel, f, i)
    vbr = [_PMs.var(pm, nw, :vbr, i) for nw in sorted_nw_ids(pm)]
    vbi = [_PMs.var(pm, nw, :vbi, i) for nw in sorted_nw_ids(pm)]
    
    cfr = [_PMs.var(pm, nw, :cfr, f) for nw in sorted_nw_ids(pm)]
    cfi = [_PMs.var(pm, nw, :cfi, f) for nw in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model,    vbr[fundamental(pm)] * cfr[fundamental(pm)] 
                                + vbi[fundamental(pm)] * cfi[fundamental(pm)] 
                                    == 
                                0.0)
    JuMP.@constraint(pm.model,  sum(vbr[n] * cfr[n] + vbi[n] * cfi[n] 
                                        for n in 2:lastindex(vbr)) 
                                    == 
                                0.0)
end