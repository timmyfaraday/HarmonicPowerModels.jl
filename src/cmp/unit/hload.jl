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

# init #########################################################################
# i_base_ka = s_base_mva / v_base_kv, see Power System Analysis, pg. 26
calc_hload_current_base(hdata::Dict{String,Any}, ldata::Dict{String,Any}) =
    hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(ldata["load_bus"])]["v_base_kv"]
# i_rms_max = S_nom / s_base_mva / sqrt(3) / v_rms_max(bus)
function calc_hload_current_rms_max(hdata::Dict{String,Any}, ldata::Dict{String,Any})
    S_nom       = sqrt(ldata["pd"]^2 + ldata["qd"]^2) 
    v_rms_max   = hdata["nw"]["1"]["bus"][string(ldata["load_bus"])]["v_rms_max"]

    v_rms_min   = hdata["nw"]["1"]["bus"][string(ldata["load_bus"])]["v_rms_min"]
     
    # return S_nom / (sqrt(3) * v_rms_max)
    return S_nom / (sqrt(3) * v_rms_min)
end
""
hload_current_magnitude_limit(pm::AbstractHarmonicModel, nw::Int, l) = 
    nw == fundamental(pm) ? _PMs.ref(pm, fundamental(pm), :hload, l, "i_rms_max") :
                            _PMs.ref(pm, fundamental(pm), :hload, l, "i_rms_max") *
                            _PMs.ref(pm, nw, :hload, l, "hcm")
""
collect_hload_current_magnitude_limits(pm::AbstractHarmonicModel, nw::Int) = 
    Dict(l => hload_current_magnitude_limit(pm, nw, l) for l in _PMs.ids(pm, nw, :hload))

# parameters ###################################################################
""
function add_hload_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, nwt) in hdata["nw"], (l, hload) in nwt["hload"]
        h       = parse(Int, nw)
        ldata   = fdata["load"][l]
        if nw == "1"
            hload["id"]         = ldata["index"]
            hload["bus"]        = ldata["load_bus"]
            #-----------------------------------#
            hload["i_base_ka"]  = calc_hload_current_base(hdata, ldata)
            hload["i_rms_max"]  = calc_hload_current_rms_max(hdata, ldata)
            hload["p_fund"]     = ldata["pd"]
            hload["q_fund"]     = ldata["qd"]
        else
            hload["hcm"]        = fdata["bus"][string(ldata["load_bus"])]["nh_$nw"]
end end end

# variables ####################################################################
""
function variable_hload_current(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_hload_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_hload_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_hload_current_real(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_hload_current_magnitude_limits(pm, nw)
    
    clr = _PMs.var(pm, nw)[:clr] = 
        JuMP.@variable( pm.model, 
                        [l in _PMs.ids(pm, nw, :hload)], 
                        base_name="$(nw)_clr",
                        start=0.0)

    for l in _PMs.ids(pm, nw, :hload)
        if bounded
            JuMP.set_lower_bound(clr[l], -c_lim[l])
            JuMP.set_upper_bound(clr[l],  c_lim[l])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :hload, :clr, _PMs.ids(pm, nw, :hload), clr)
end
""
function variable_hload_current_imaginary(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_hload_current_magnitude_limits(pm, nw)
    
    cli = _PMs.var(pm, nw)[:cli] = 
            JuMP.@variable( pm.model,
                            [l in _PMs.ids(pm, nw, :hload)], 
                            base_name="$(nw)_cli",
                            start=0.0)

    for l in _PMs.ids(pm, nw, :hload)
        if bounded
            JuMP.set_lower_bound(cli[l], -c_lim[l])
            JuMP.set_upper_bound(cli[l],  c_lim[l])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :hload, :cli, _PMs.ids(pm, nw, :hload), cli)
end

# constraints ##################################################################
""
function constraint_hload_power(pm::_PMs.AbstractPowerModel, l::Int; nw::Int=fundamental(pm))
    if nw == fundamental(pm)
        i       = _PMs.ref(pm, nw, :hload, l, "bus")
        
        p_fund  = _PMs.ref(pm, nw, :hload, l, "p_fund")
        q_fund  = _PMs.ref(pm, nw, :hload, l, "q_fund")

        constraint_hload_constant_power(pm, nw, l, i, p_fund, q_fund)
    else
        hcm = _PMs.ref(pm, nw, :hload, l, "hcm")
        
        constraint_hload_constant_current(pm, nw, l, hcm)
    end
end
""
function constraint_hload_constant_power(pm::AbstractHarmonicModel, n::Int, l, i, p_fund, q_fund)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    clr = _PMs.var(pm, n, :clr, l)
    cli = _PMs.var(pm, n, :cli, l)

    JuMP.@constraint(pm.model, p_fund == vbr * clr  + vbi * cli)
    JuMP.@constraint(pm.model, q_fund == vbi * clr  - vbr * cli)
end
""
function constraint_hload_constant_current(pm::_PMs.AbstractIVRModel, n::Int, l, hcm)
    clr     = _PMs.var(pm, n, :clr, l)
    cli     = _PMs.var(pm, n, :cli, l)

    clfr    = _PMs.var(pm, fundamental(pm), :clr, l)
    clfi    = _PMs.var(pm, fundamental(pm), :cli, l)

    JuMP.@constraint(pm.model, clr == hcm * clfr)
    JuMP.@constraint(pm.model, cli == hcm * clfi)
end