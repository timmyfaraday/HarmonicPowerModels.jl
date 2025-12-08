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
# i_rms_max = S_nom / s_base_mva / sqrt(3) / min(v_rms_max(bus_fr), v_rms_max(bus_to))
function calc_branch_current_rms_max(hdata::Dict{String,Any}, bdata::Dict{String,Any})
    S_nom       = bdata["rate_a"] 
    v_rms_max   = min(hdata["nw"]["1"]["bus"][string(bdata["f_bus"])]["v_rms_max"],
                      hdata["nw"]["1"]["bus"][string(bdata["t_bus"])]["v_rms_max"])
    
    v_rms_min   = min(hdata["nw"]["1"]["bus"][string(bdata["f_bus"])]["v_rms_min"],
                      hdata["nw"]["1"]["bus"][string(bdata["t_bus"])]["v_rms_min"])
     
    #return S_nom / (sqrt(3) * v_rms_max)
    #return S_nom / (sqrt(3) * v_rms_min)

    return S_nom / (v_rms_min)
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
            branch["id"]            = bdata["index"]
            branch["bus"]           = [bdata["f_bus"], bdata["t_bus"]]
            #-----------------------------------#
            branch["r"]             = bdata["br_r"]
            branch["x"]             = bdata["br_x"]
            branch["g_fr"]          = bdata["g_fr"]
            branch["b_fr"]          = bdata["b_fr"]
            branch["g_to"]          = bdata["g_fr"]
            branch["b_to"]          = bdata["b_to"]
            #-----------------------------------#
            branch["i_base_ka"]     = calc_branch_current_base(hdata, bdata)
            branch["i_fund_magn"]   = [0.0, 0.0]
            branch["i_rms_max"]     = calc_branch_current_rms_max(hdata, bdata)
        else
            branch["r"]             = bdata["br_r"] * sqrt(h)
            branch["x"]             = bdata["br_x"] * h
            branch["g_fr"]          = bdata["g_fr"] / sqrt(h)
            branch["b_fr"]          = bdata["b_fr"] * h^(sign(bdata["b_fr"]))
            branch["g_to"]          = bdata["g_fr"] / sqrt(h)
            branch["b_to"]          = bdata["b_to"] * h^(sign(bdata["b_to"]))
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
                                [(b,i,j) in _PMs.ref(pm, nw, :arcs_branch)], 
                                base_name="$(nw)_cbr",
                                start=0.0)

    if bounded
        for (b,i,j) in _PMs.ref(pm, nw, :arcs_branch)
            JuMP.set_lower_bound(cbr[(b,i,j)], -c_lim[b])
            JuMP.set_upper_bound(cbr[(b,i,j)],  c_lim[b])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :cbr_fr, :cbr_to, _PMs.ref(pm, nw, :arcs_branch_from), _PMs.ref(pm, nw, :arcs_branch_to), cbr)
end
""
function variable_branch_current_imaginary(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_branch_current_magnitude_limits(pm, nw)

    cbi     = _PMs.var(pm, nw)[:cbi] = 
                JuMP.@variable( pm.model,
                                [(b,i,j) in _PMs.ref(pm, nw, :arcs_branch)], 
                                base_name="$(nw)_cbi",
                                start=0.0)

    if bounded
        for (b,i,j) in _PMs.ref(pm, nw, :arcs_branch)
            JuMP.set_lower_bound(cbi[(b,i,j)], -c_lim[b])
            JuMP.set_upper_bound(cbi[(b,i,j)],  c_lim[b])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :cbi_fr, :cbi_to, _PMs.ref(pm, nw, :arcs_branch_from), _PMs.ref(pm, nw, :arcs_branch_to), cbi)
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

""
function variable_branch_power(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_power_active(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_power_reactive(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_branch_power_active(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    pb = Dict()
    for (b,i,j) in _PMs.ref(pm, nw, :arcs_branch)
        vbr_fr = _PMs.var(pm, nw, :vbr, i)
        vbi_fr = _PMs.var(pm, nw, :vbi, i)
        cbr_fr = _PMs.var(pm, nw, :cbr, (b,i,j))
        cbi_fr = _PMs.var(pm, nw, :cbi, (b,i,j))

        vbr_to = _PMs.var(pm, nw, :vbr, j)
        vbi_to = _PMs.var(pm, nw, :vbi, j)
        cbr_to = _PMs.var(pm, nw, :cbr, (b,j,i))
        cbi_to = _PMs.var(pm, nw, :cbi, (b,j,i))
        
        pb[(b,i,j)] = vbr_fr * cbr_fr  + vbi_fr * cbi_fr
        pb[(b,j,i)] = vbr_to * cbr_to  + vbi_to * cbi_to
    end
    _PMs.var(pm, nw)[:pb] = pb
    
    report && _PMs.sol_component_value_edge(pm, nw, :branch, :pb_fr, :pb_to, _PMs.ref(pm, nw, :arcs_branch_from), _PMs.ref(pm, nw, :arcs_branch_to), pb)
end
""
function variable_branch_power_reactive(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    qb = Dict()
    for (b,i,j) in _PMs.ref(pm, nw, :arcs_branch)
        vbr_fr = _PMs.var(pm, nw, :vbr, i)
        vbi_fr = _PMs.var(pm, nw, :vbi, i)
        cbr_fr = _PMs.var(pm, nw, :cbr, (b,i,j))
        cbi_fr = _PMs.var(pm, nw, :cbi, (b,i,j))

        vbr_to = _PMs.var(pm, nw, :vbr, j)
        vbi_to = _PMs.var(pm, nw, :vbi, j)
        cbr_to = _PMs.var(pm, nw, :cbr, (b,j,i))
        cbi_to = _PMs.var(pm, nw, :cbi, (b,j,i))
        
        qb[(b,i,j)] = vbi_fr * cbr_fr  - vbr_fr * cbi_fr
        qb[(b,j,i)] = vbi_to * cbr_to  - vbr_to * cbi_to
    end
    _PMs.var(pm, nw)[:qb] = qb
    
    report && _PMs.sol_component_value_edge(pm, nw, :branch, :qb_fr, :qb_to, _PMs.ref(pm, nw, :arcs_branch_from), _PMs.ref(pm, nw, :arcs_branch_to), qb)
end

# constraints ##################################################################
## branch current constraint ###################################################
""
function constraint_branch_current_from(pm::AbstractHarmonicModel, i::Int; nw::Int=fundamental(pm))
    idx     = (i, _PMs.ref(pm, fundamental(pm), :branch, i, "bus")...)

    g       = _PMs.ref(pm, nw, :branch, i, "g_fr")
    b       = _PMs.ref(pm, nw, :branch, i, "b_fr")

    constraint_branch_current(pm, nw, idx, g, b, 1)
end
""
function constraint_branch_current_to(pm::AbstractHarmonicModel, i::Int; nw::Int=fundamental(pm))
    idx     = (i, reverse(_PMs.ref(pm, fundamental(pm), :branch, i, "bus"))...)

    g       = _PMs.ref(pm, nw, :branch, i, "g_to")
    b       = _PMs.ref(pm, nw, :branch, i, "b_to")

    constraint_branch_current(pm, nw, idx, g, b, -1)
end
""
function constraint_branch_current(pm::AbstractHarmonicModel, n::Int, idx, g, b, sign)
    vbr     = _PMs.var(pm, n, :vbr, idx[2])
    vbi     = _PMs.var(pm, n, :vbi, idx[2])

    cbsr    = _PMs.var(pm, n, :cbsr, idx[1]) * sign
    cbsi    = _PMs.var(pm, n, :cbsi, idx[1]) * sign

    cbr     = _PMs.var(pm, n, :cbr, idx)
    cbi     = _PMs.var(pm, n, :cbi, idx)

    JuMP.@constraint(pm.model, cbr == cbsr + g * vbr - b * vbi)
    JuMP.@constraint(pm.model, cbi == cbsi + g * vbi + b * vbr)
end
## voltage drop constraint #####################################################
""
function constraint_branch_voltage_drop(pm::AbstractHarmonicModel, i::Int; nw::Int=fundamental(pm))
    idx     = (i, _PMs.ref(pm, fundamental(pm), :branch, i, "bus")...)

    r       = _PMs.ref(pm, nw, :branch, i, "r")
    x       = _PMs.ref(pm, nw, :branch, i, "x")

    constraint_branch_voltage_drop(pm, nw, idx, r, x)
end
""
function constraint_branch_voltage_drop(pm::AbstractHarmonicModel, n::Int, idx, r, x)
    cbsr    = _PMs.var(pm, n, :cbsr, idx[1])
    cbsi    = _PMs.var(pm, n, :cbsi, idx[1])
    
    vbr_fr  = _PMs.var(pm, n, :vbr, idx[2])
    vbi_fr  = _PMs.var(pm, n, :vbi, idx[2])

    vbr_to  = _PMs.var(pm, n, :vbr, idx[3])
    vbi_to  = _PMs.var(pm, n, :vbi, idx[3])

    JuMP.@constraint(pm.model, vbr_to == vbr_fr - r * cbsr + x * cbsi)
    JuMP.@constraint(pm.model, vbi_to == vbi_fr - r * cbsi - x * cbsr)
end
## root-mean-square current limit ##############################################
""
function constraint_branch_current_rms_limit(pm::AbstractHarmonicModel, i::Int)
    branch      = _PMs.ref(pm, fundamental(pm), :branch, i)
    idx_fr      = (i, _PMs.ref(pm, fundamental(pm), :branch, i, "bus")...)
    idx_to      = (i, reverse(_PMs.ref(pm, fundamental(pm), :branch, i, "bus"))...)

    i_rms_max   = _PMs.ref(pm, fundamental(pm), :branch, i, "i_rms_max")
    i_fund_magn = _PMs.ref(pm, fundamental(pm), :branch, i, "i_fund_magn")

    constraint_branch_current_rms_limit(pm, idx_fr, idx_to, i_rms_max, i_fund_magn)
end
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

    JuMP.@constraint(pm.model, [sqrt(i_rms_max^2 - i_fund_magn[1]^2); vcat(cbr_fr, cbi_fr)] in JuMP.SecondOrderCone())
    JuMP.@constraint(pm.model, [sqrt(i_rms_max^2 - i_fund_magn[2]^2); vcat(cbr_to, cbi_to)] in JuMP.SecondOrderCone())
end