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

# Open topics
# 1) the bounds on the harmonic currents are not very tight, look into standard
#    to find tightning bounds. 

# util #########################################################################
""
# i_base_ka = s_base_mva / v_base_kv, see Power System Analysis, pg. 26
calc_branch_current_base(hdata::Dict{String,Any}, bdata::Dict{String,Any}) =
    hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(bdata["f_bus"])]["v_base_kv"]
""
# i_rms_max = S_nom / s_base_mva / sqrt(3) / min(v_rms_max(f_bus), v_rms_max(t_bus))
function calc_branch_current_rms_max(hdata::Dict{String,Any}, bdata::Dict{String,Any})
    S_nom       = bdata["rate_a"] 
    s_base_mva  = hdata["s_base_mva"]
    v_rms_max   = min(hdata["nw"]["1"]["bus"][string(bdata["f_bus"])]["v_rms_max"],
                      hdata["nw"]["1"]["bus"][string(bdata["t_bus"])]["v_rms_max"])
     
    return S_nom / s_base_mva / sqrt(3) / v_rms_max
end
""
collect_branch_current_magnitude_limits(pm::_PMs.AbstractPowerModel, nw::Int) = 
    Dict(b => _PMs.ref(pm, fundamental(pm), :branch, b, "i_rms_max") 
            for b in _PMs.ids(pm, nw, :branch))

# parameters ###################################################################
""
function add_branch_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (nb, branch) in ntw["branch"]
        h       = parse(Int, nw)
        bdata   = fdata["branch"][nb]
        if nw == "1"
            branch = Dict(  "id"            => bdata["index"],
                            "f_bus"         => bdata["f_bus"],
                            "t_bus"         => bdata["t_bus"],
                            #-----------------------------------#
                            "r"             => bdata["br_r"],
                            "x"             => bdata["br_x"],
                            "g_fr"          => bdata["g_fr"],
                            "b_fr"          => bdata["b_fr"],
                            "g_to"          => bdata["g_fr"],
                            "b_to"          => bdata["b_to"],
                            #-----------------------------------#
                            "i_base_ka"     => calc_branch_current_base(hdata, bdata),
                            "i_fund_magn"   => [0.0, 0.0],
                            "i_rms_max"     => calc_branch_current_rms_max(hdata, bdata))
        else
            branch = Dict(  "r"             => bdata["br_r"] * sqrt(h),
                            "x"             => bdata["br_x"] * h,
                            "g_fr"          => bdata["g_fr"] / sqrt(h),
                            "b_fr"          => bdata["b_fr"] / h,
                            "g_to"          => bdata["g_fr"] / sqrt(h),
                            "b_to"          => bdata["b_to"] / h)
end end end

# variables ####################################################################
""
function variable_branch_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_branch_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_branch_current_magnitude_limits(pm, nw)

    cbr     = _PMs.var(pm, nw)[:cbr] = 
                JuMP.@variable( pm.model,
                                [(b,i,j) in _PMs.ref(pm, nw, :arcs)], 
                                base_name="$(nw)_cbr",
                                start=0.0)

    if bounded
        for (b,i,j) in _PMs.ref(pm, nw, :arcs)
            JuMP.set_lower_bound(cr[(b,i,j)], -c_lim[b])
            JuMP.set_upper_bound(cr[(b,i,j)],  c_lim[b])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :cbr_fr, :cbr_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), cbr)
end
""
function variable_branch_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_branch_current_magnitude_limits(pm, nw)

    cbi     = _PMs.var(pm, nw)[:cbi] = 
                JuMP.@variable( pm.model,
                                [(b,i,j) in _PMs.ref(pm, nw, :arcs)], 
                                base_name="$(nw)_ci",
                                start=0.0)

    if bounded
        for (b,i,j) in _PMs.ref(pm, nw, :arcs)
            JuMP.set_lower_bound(ci[(b,i,j)], -c_lim[b])
            JuMP.set_upper_bound(ci[(b,i,j)],  c_lim[b])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :cbi_fr, :cbi_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), cbi)
end
""
function variable_branch_series_current_real(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_branch_current_magnitude_limits(pm, nw)

    cbsr    = _PMs.var(pm, nw)[:cbsr] = 
                JuMP.@variable( pm.model,
                                [b in _PMs.ids(pm, nw, :branch)], 
                                base_name="$(nw)_cbsr",
                                start=0.0)

    if bounded
        for b in _PMs.ids(pm, nw, :branch)
            JuMP.set_lower_bound(cbsr[b], -c_lim[b])
            JuMP.set_upper_bound(cbsr[b],  c_lim[b])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :branch, :cbsr_fr, _PMs.ids(pm, nw, :branch), cbsr)
end
""
function variable_branch_series_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_branch_current_magnitude_limits(pm, nw)

    cbsi    = _PMs.var(pm, nw)[:cbsi] = 
                JuMP.@variable( pm.model,
                                [b in _PMs.ids(pm, nw, :branch)], 
                                base_name="$(nw)_cbsi",
                                start=0.0)

    if bounded
        for b in _PMs.ids(pm, nw, :branch)
            JuMP.set_lower_bound(cbsi[b], -c_lim[b])
            JuMP.set_upper_bound(cbsi[b],  c_lim[b])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :branch, :cbsi_fr, _PMs.ids(pm, nw, :branch), cbsi)
end

# constraints ##################################################################
## branch current constraint ###################################################
""
function constraint_branch_current_from(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    branch  = _PMs.ref(pm, nw, :branch, i)

    idx     = (i, branch["f_bus"], branch["t_bus"])

    g       = branch["g_fr"]
    b       = branch["b_fr"]

    constraint_branch_current(pm, nw, idx, g, b, 1)
end
""
function constraint_branch_current_to(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    branch  = _PMs.ref(pm, nw, :branch, i)
    
    idx     = (i, branch["t_bus"], branch["f_bus"])

    g       = branch["g_to"]
    b       = branch["b_to"]

    constraint_branch_current(pm, nw, idx, g, b, -1)
end
""
function constraint_branch_current(pm::HarmonicPowerModel, n::Int, idx, g, b, sign)
    vr      = _PMs.var(pm, n, :vr, idx[2])
    vi      = _PMs.var(pm, n, :vi, idx[2])

    cbsr    = _PMs.var(pm, n, :cbsr, idx[1]) * sign
    cbsi    = _PMs.var(pm, n, :cbsi, idx[1]) * sign

    cbr     = _PMs.var(pm, n, :cbr, idx)
    cbi     = _PMs.var(pm, n, :cbi, idx)

    JuMP.@constraint(pm.model, cbr == cbsr + g * vr - b * vi)
    JuMP.@constraint(pm.model, cbi == cbsi + g * vi + b * vr)
end
## voltage drop constraint #####################################################
""
function constraint_branch_voltage_drop(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    branch  = _PMs.ref(pm, nw, :branch, i)

    idx     = (i, branch["f_bus"], branch["t_bus"])

    r       = branch["r"]
    x       = branch["x"]

    constraint_voltage_drop(pm, nw, i, idx, r, x)
end
""
function constraint_branch_voltage_drop(pm::HarmonicPowerModel, n::Int, idx, r, x)
    cbsr    = _PMs.var(pm, n, :cbsr, idx[1])
    cbsi    = _PMs.var(pm, n, :cbsi, idx[1])
    
    vr_fr   = _PMs.var(pm, n, :vr, idx[2])
    vi_fr   = _PMs.var(pm, n, :vi, idx[2])

    vr_to   = _PMs.var(pm, n, :vr, idx[3])
    vi_to   = _PMs.var(pm, n, :vi, idx[3])

    JuMP.@constraint(pm.model, vr_to == vr_fr - r * cbsr + x * cbsi)
    JuMP.@constraint(pm.model, vi_to == vi_fr - r * cbsi - x * cbsr)
end
## root-mean-square current limit ##############################################
""
function constraint_branch_current_rms_limit(pm::HarmonicPowerModel, i::Int)
    branch      = _PMs.ref(pm, fundamental(pm), :branch, i)
    idx_fr      = (i, branch["f_bus"], branch["t_bus"])
    idx_to      = (i, branch["t_bus"], branch["f_bus"])

    i_rms_max   = branch["i_rms_max"]
    i_fund_magn = branch["i_fund_magn"]

    constraint_branch_current_rms_limit(pm, idx_fr, idx_to, i_rms_max, i_fund_magn)
end
# branch
""
function constraint_branch_current_rms_limit(pm::HarmonicPowerModel, idx_fr, idx_to, i_rms_max, i_fund_magn)
    cbr_fr  = [_PMs.var(pm, n, :cbr, idx_fr) for n in sorted_nw_ids(pm)]
    cbi_fr  = [_PMs.var(pm, n, :cbi, idx_fr) for n in sorted_nw_ids(pm)]

    cbr_to  = [_PMs.var(pm, n, :cbr, idx_to) for n in sorted_nw_ids(pm)]
    cbi_to  = [_PMs.var(pm, n, :cbi, idx_to) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(cbr_fr.^2 + cbi_fr.^2) <= i_rms_max^2)
    JuMP.@constraint(pm.model, sum(cbr_to.^2 + cbi_to.^2) <= i_rms_max^2)
end
""
function constraint_branch_current_rms_limit(pm::dHHCPowerModel, idx_fr, idx_to, i_rms_max, i_fund_magn)
    cbr_fr  = [_PMs.var(pm, n, :cbr, idx_fr) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cbi_fr  = [_PMs.var(pm, n, :cbi, idx_fr) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    cbr_to  = [_PMs.var(pm, n, :cbr, idx_to) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cbi_to  = [_PMs.var(pm, n, :cbi, idx_to) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(i_rms_max^2 - i_fund_magn[1]^2); vcat(cbr_to, cbi_fr)] in JuMP.SecondOrderCone())
    JuMP.@constraint(pm.model, [sqrt(i_rms_max^2 - i_fund_magn[2]^2); vcat(cbr_to, cbi_to)] in JuMP.SecondOrderCone())
end