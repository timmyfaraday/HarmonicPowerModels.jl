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
calc_hscr_current_base(hdata::Dict{String,Any}, sdata::Dict{String,Any}) =
    hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(sdata["load_bus"])]["v_base_kv"]
""
calc_hsrc_current_angle_ref(hdata::Dict{String,Any}, sdata::Dict{String,Any}, h) =
    if haskey(hdata["nw"]["1"]["bus"][string(sdata["load_bus"])], "ref_angle")
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
function add_hscr_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (s, hscr) in ntw["hsrc"]
        h       = parse(Int, nw)
        sdata   = fdata["load"][s]
        if nw == "1"
            hscr["id"]          = sdata["index"]
            hscr["bus"]         = sdata["load_bus"]
            #-----------------------------------#
            hscr["p_fund"]      = sdata["pd"]
            hscr["q_fund"]      = sdata["qd"]
            #-----------------------------------#
            hscr["i_base_ka"]   = calc_hscr_current_base(hdata, ldata)
            hscr["csar"]        = calc_hsrc_current_angle_ref(hdata, sdata, h)
        else
            hscr["csar"]        = calc_hsrc_current_angle_ref(hdata, sdata, h)
end end end

# variables ####################################################################
""
function variable_hsrc_current(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_hsrc_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_hsrc_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    if nw ≠ fundamental(pm)
        variable_hsrc_current_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    end
end
""
function variable_hsrc_current_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csr = _PMs.var(pm, nw)[:csr] = 
        JuMP.@variable( pm.model, 
                        [s in _PMs.ids(pm, nw, :hsrc)], 
                        base_name="$(nw)_csr",
                        start=0.0)

    report && _PMs.sol_component_value(pm, nw, :hsrc, :csr, _PMs.ids(pm, nw, :hsrc), csr)
end
""
function variable_hsrc_current_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csi = _PMs.var(pm, nw)[:csi] = 
            JuMP.@variable( pm.model,
                            [s in _PMs.ids(pm, nw, :hsrc)], 
                            base_name="$(nw)_csi",
                            start=0.0)

    report && _PMs.sol_component_value(pm, nw, :hscr, :csi, _PMs.ids(pm, nw, :hscr), csi)
end
""
function variable_hsrc_current_magnitude(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csm = _PMs.var(pm, nw)[:csm] = 
        JuMP.@variable( pm.model,
                        [s in _PMs.ids(pm, nw, :hsrc)], 
                        base_name="$(nw)_csm",
                        start=0.0)

    for s in _PMs.ids(pm, nw, :hsrc)
        if bounded
            JuMP.set_lower_bound(csm[s], 0.0)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :hsrc, :csm, _PMs.ids(pm, nw, :hsrc), csm)
end

# constraints ##################################################################
## hsrc current constraint #####################################################
""
function constraint_hsrc_current(pm::HarmonicPowerModel, s::Int; nw::Int=fundamental(pm))
    if nw == fundamental(pm)
        i       = _PMs.ref(pm, nw, :hsrc, s, "bus")
        
        p_fund  = _PMs.ref(pm, nw, :hsrc, s, "p_fund")
        q_fund  = _PMs.ref(pm, nw, :hsrc, s, "q_fund")

        constraint_hsrc_constant_power(pm, nw, s, i, p_fund, q_fund)
    else
        csar = _PMs.ref(pm, nw, :hscr, s, "csar")

        constraint_hsrc_current_angle_ref(pm, nw, s, csar)
end end
""
function constraint_hsrc_constant_power(pm::HarmonicPowerModel, n::Int, s, i, p_fund, q_fund)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    csr = _PMs.var(pm, n, :csr, s)
    csi = _PMs.var(pm, n, :csi, s)

    JuMP.@constraint(pm.model, p_fund == vbr * csr  + vbi * csi)
    JuMP.@constraint(pm.model, q_fund == vbi * csr  - vbr * csi)
end
""
function constraint_hsrc_current_angle_ref(pm::HarmonicPowerModel, n::Int, s, csar)
    csr = _PMs.var(pm, n, :csr, s)
    csi = _PMs.var(pm, n, :csi, s)
    csm = _PMs.var(pm, n, :csm, s)

    JuMP.@constraint(pm.model, csm * sind(csar) == csi)
    JuMP.@constraint(pm.model, csm * cosd(csar) == csr)
end