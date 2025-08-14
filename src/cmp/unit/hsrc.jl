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

# util #########################################################################
# i_base_ka = s_base_mva / v_base_kv, see Power System Analysis, pg. 26
calc_hsrc_current_base(hdata::Dict{String,Any}, rdata::Dict{String,Any}) =
    hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(rdata["load_bus"])]["v_base_kv"]
""
calc_hsrc_current_angle_ref(hdata::Dict{String,Any}, rdata::Dict{String,Any}, h) =
    if haskey(hdata["nw"]["1"]["bus"][string(rdata["load_bus"])], "ref_angle")
        if is_pos_sequence(h)
            return  bus["ref_angle"]
        elseif is_neg_sequence(h)
            return -bus["ref_angle"]
        elseif is_zero_sequence(h)
            return  0.0
        end
    else 
        return 0.0
    end

# parameters ###################################################################
""
function add_hsrc_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (r, hsrc) in ntw["hsrc"]
        h       = parse(Int, nw)
        rdata   = fdata["load"][r]
        if nw == "1"
            hsrc["id"]          = rdata["index"]
            hsrc["bus"]         = rdata["load_bus"]
            #-----------------------------------#
            hsrc["p_fund"]      = rdata["pd"]
            hsrc["q_fund"]      = rdata["qd"]
            #-----------------------------------#
            hsrc["i_base_ka"]   = calc_hsrc_current_base(hdata, rdata)
            hsrc["crar"]        = calc_hsrc_current_angle_ref(hdata, rdata, h)
        else
            hsrc["crar"]        = calc_hsrc_current_angle_ref(hdata, rdata, h)
end end end

# variables ####################################################################
""
function variable_hsrc_current(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_hsrc_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_hsrc_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    if nw ≠ fundamental(pm)
        variable_hsrc_current_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    end
end
""
function variable_hsrc_current_real(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    crr = _PMs.var(pm, nw)[:crr] = 
        JuMP.@variable( pm.model, 
                        [r in _PMs.ids(pm, nw, :hsrc)], 
                        base_name="$(nw)_crr",
                        start=0.0)

    report && _PMs.sol_component_value(pm, nw, :hsrc, :crr, _PMs.ids(pm, nw, :hsrc), crr)
end
""
function variable_hsrc_current_imaginary(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cri = _PMs.var(pm, nw)[:cri] = 
            JuMP.@variable( pm.model,
                            [r in _PMs.ids(pm, nw, :hsrc)], 
                            base_name="$(nw)_cri",
                            start=0.0)

    report && _PMs.sol_component_value(pm, nw, :hsrc, :cri, _PMs.ids(pm, nw, :hsrc), cri)
end
""
function variable_hsrc_current_magnitude(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    crm = _PMs.var(pm, nw)[:crm] = 
        JuMP.@variable( pm.model,
                        [r in _PMs.ids(pm, nw, :hsrc)], 
                        base_name="$(nw)_crm",
                        start=0.0)

    for r in _PMs.ids(pm, nw, :hsrc)
        if bounded
            JuMP.set_lower_bound(crm[r], 0.0)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :hsrc, :crm, _PMs.ids(pm, nw, :hsrc), crm)
end

# constraints ##################################################################
## hsrc current constraint #####################################################
""
function constraint_hsrc_current(pm::AbstractHarmonicModel, r::Int; nw::Int=fundamental(pm))
    if nw == fundamental(pm)
        i       = _PMs.ref(pm, nw, :hsrc, r, "bus")
        
        p_fund  = _PMs.ref(pm, nw, :hsrc, r, "p_fund")
        q_fund  = _PMs.ref(pm, nw, :hsrc, r, "q_fund")

        constraint_hsrc_constant_power(pm, nw, r, i, p_fund, q_fund)
    else
        crar = _PMs.ref(pm, nw, :hsrc, r, "crar")

        constraint_hsrc_current_angle_ref(pm, nw, r, crar)
end end
""
function constraint_hsrc_constant_power(pm::AbstractHarmonicModel, n::Int, r, i, p_fund, q_fund)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    crr = _PMs.var(pm, n, :crr, r)
    cri = _PMs.var(pm, n, :cri, r)

    JuMP.@constraint(pm.model, p_fund == vbr * crr + vbi * cri)
    JuMP.@constraint(pm.model, q_fund == vbi * crr - vbr * cri)
end
""
function constraint_hsrc_current_angle_ref(pm::AbstractHarmonicModel, n::Int, r, crar)
    crr = _PMs.var(pm, n, :crr, r)
    cri = _PMs.var(pm, n, :cri, r)
    crm = _PMs.var(pm, n, :crm, r)

    JuMP.@constraint(pm.model, crm * sind(crar) == cri)
    JuMP.@constraint(pm.model, crm * cosd(crar) == crr)
end