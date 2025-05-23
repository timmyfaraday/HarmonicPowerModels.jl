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

# util #########################################################################
""
# i_base_ka = s_base_mva / v_base_kv, see Power System Analysis, pg. 26
calc_filter_current_base(hdata::Dict{String,Any}, idata::Dict{String,Any}) =
    hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(idata["bus"])]["v_base_kv"] 
""
# i_rms_max = S_nom / s_base_mva / sqrt(3) / v_rms_max(bus)
function calc_filter_current_rms_max(hdata::Dict{String,Any}, idata::Dict{String,Any})
    S_nom       = idata["rate_a"] 
    s_base_mva  = hdata["s_base_mva"]
    v_rms_max   = hdata["nw"]["1"]["bus"][string(idata["bus"])]["v_rms_max"]
     
    return S_nom / s_base_mva / sqrt(3) ./ v_rms_max
end
""
collect_filter_current_magnitude_limits(pm::AbstractHarmonicModel, nw::Int) = 
    Dict(f => _PMs.ref(pm, fundamental(pm), :filter, f, "i_rms_max")
            for f in _PMs.ids(pm, nw, :filter))

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
            #-----------------------------------#
            filter["i_base_ka"] = calc_filter_current_base(hdata, idata)
            filter["i_rms_max"] = calc_filter_current_rms_max(hdata, idata)
end end end

# variables ####################################################################
""
function variable_filter_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_filter_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_filter_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_filter_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_filter_current_magnitude_limits(pm, nw)

    cfr = _PMs.var(pm, nw)[:cfr] = 
            JuMP.@variable( pm.model,
                            [f in _PMs.ids(pm, nw, :filter)], 
                            base_name="$(nw)_cfr",
                            start=0.0)

    if bounded
        for f in _PMs.ids(pm, nw, :filter)
            JuMP.set_lower_bound(cfr[f], -c_lim[f])
            JuMP.set_upper_bound(cfr[f],  c_lim[f])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :filter, :cfr, _PMs.ids(pm, nw, :filter), cfr)
end
""
function variable_filter_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_filter_current_magnitude_limits(pm, nw)

    cfi = _PMs.var(pm, nw)[:cfi] = 
            JuMP.@variable( pm.model,
                            [f in _PMs.ids(pm, nw, :filter)], 
                            base_name="$(nw)_cfi",
                            start=0.0)

    if bounded
        for f in _PMs.ids(pm, nw, :filter)
            JuMP.set_lower_bound(cfi[f], -c_lim[f])
            JuMP.set_upper_bound(cfi[f],  c_lim[f])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :filter, :cfi, _PMs.ids(pm, nw, :filter), cfi)
end

# constraints ##################################################################
## active filter current constraint ############################################
""
function constraint_filter_current(pm::_PMs.AbstractPowerModel, f::Int)
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

## filter root-mean-square current limit ####################################
""
function constraint_filter_current_rms_limit(pm::AbstractHarmonicModel, f::Int)
    i_rms_max   = _PMs.ref(pm, fundamental(pm), :filter, f, "i_rms_max")

    constraint_filter_current_rms_limit(pm, f, i_rms_max)
end
function constraint_filter_current_rms_limit(pm::HarmonicPowerModel, f, i_rms_max)
    cfr =  [_PMs.var(pm, n, :cfr, f) for n in sorted_nw_ids(pm)]
    cfi =  [_PMs.var(pm, n, :cfi, f) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(cfr.^2 + cfi.^2) <= i_rms_max^2)
end
""
function constraint_filter_current_rms_limit(pm::dHHCPowerModel, f, i_rms_max)
    cfr =  [_PMs.var(pm, n, :cfr, f) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cfi =  [_PMs.var(pm, n, :cfi, f) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [i_rms_max; vcat(cfr, cfi)] in JuMP.SecondOrderCone())
end