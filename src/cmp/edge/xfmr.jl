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
# 1) the bounds on the harmonic currents are not very tight, look into standards
#    to find tightning bounds.

# util #########################################################################
""
calc_xfmr_configuration(hdata::Dict{String,Any}, xdata::Dict{String,Any}) = 
    haskey(xdata, "vg") ? [uppercase(xdata["vg"][n]) for n in 1:2] : ['Y', 'Y'] ;
""
calc_xfmr_grounding(hdata::Dict{String,Any}, xdata::Dict{String,Any}) = 
    haskey(xdata, "gnd1") && haskey(xdata, "gnd2") ? [xdata["gnd1"], xdata["gnd2"]] : Bool[0, 0] ;
""
function calc_xfmr_shift_real(hdata::Dict{String,Any}, xdata::Dict{String,Any}, h::Real)
    shift = haskey(xdata, "vg") ? parse(Int, xdata["vg"][3]) : 0 ;
    if is_pos_sequence(h)
        return cosd(-30.0 * shift)
    elseif is_neg_sequence(h)
        return cosd(30.0 * shift)
    elseif is_zero_sequence(h)
        return 1.0
end end
""
function calc_xfmr_shift_imaginary(hdata::Dict{String,Any}, xdata::Dict{String,Any}, h::Real)
    shift = haskey(xdata, "vg") ? parse(Int, xdata["vg"][3]) : 0 ;
    if is_pos_sequence(h)
        return sind(-30.0 * shift)
    elseif is_neg_sequence(h)
        return sind(30.0 * shift)
    elseif is_zero_sequence(h)
        return 0.0
end end
""
# i_base_ka = s_base_mva / v_base_kv, see Power System Analysis, pg. 26
calc_xfmr_current_base(hdata::Dict{String,Any}, xdata::Dict{String,Any}) =
    [hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(nb)]["v_base_kv"] 
        for nb in [xdata["f_bus"], xdata["t_bus"]]]
""
# i_rms_max = S_nom / s_base_mva / sqrt(3) / min(v_rms_max(f_bus), v_rms_max(t_bus))
function calc_xfmr_current_rms_max(hdata::Dict{String,Any}, xdata::Dict{String,Any})
    S_nom       = xdata["rate_a"] 
    s_base_mva  = hdata["s_base_mva"]
    v_rms_max   = [ hdata["nw"]["1"]["bus"][string(xdata["f_bus"])]["v_rms_max"],
                    hdata["nw"]["1"]["bus"][string(xdata["t_bus"])]["v_rms_max"]]
     
    return S_nom / s_base_mva / sqrt(3) ./ v_rms_max
end

""
collect_xfmr_voltage_magnitude_limits(pm::HarmonicPowerModel, nw::Int) =
    Dict(x => Dict(i => ifelse( nw == fundamental(pm),
                                _PMs.ref(pm, nw, :bus, i, "v_rms_max"), 
                                _PMs.ref(pm, nw, :bus, i, "v_ihd_max") * 
                                fundamental_bus_voltage_multiplier(pm, i))
                    for i in _PMs.ref(pm, fundamental(pm), :xfmr, x, "bus"))
            for x in _PMs.ids(pm, nw, :xfmr))
""
collect_xfmr_current_magnitude_limits(pm::HarmonicPowerModel, nw::Int) = 
    Dict(x => Dict(i => _PMs.ref(pm, fundamental(pm), :xfmr, x, "i_rms_max")[ni]
            for (ni,i) in enumerate(_PMs.ref(pm, fundamental(pm), :xfmr, x, "bus")))
            for x in _PMs.ids(pm, nw, :xfmr))

# parameters ###################################################################
""
function add_xfmr_hdata!(hdata::Dict{String,Any}, 
                         fdata::Dict{String,Any}, 
                         xfmr_magn::Dict{String,Any})

    for (nw, ntw) in hdata["nw"], (nx, xfmr) in ntw["xfmr"]
        h       = parse(Int, nw)
        xdata   = fdata["xfmr"][nx]
        if nw == "1"
            xfmr["id"]          = bdata["index"]
            xfmr["Nw"]          = 2
            xfmr["bus"]         = [bdata["f_bus"], bdata["t_bus"]]
            xfmr["linear_magn"] = isempty(xfmr_magn)
            xfmr["cnf"]         = calc_xfmr_configuration(hdata, xdata)
            xfmr["gnd"]         = calc_xfmr_grounding(hdata, xdata)
            xfmr["Hᴵ"]          = Int[]
            xfmr["Hᴱ"]          = Int[]
            #-----------------------------------#
            xfmr["x_core"]      = [xdata["xsc"]]
            xfmr["b_core"]      = 0.0
            xfmr["g_core"]      = xdata["gsh"]
            #-----------------------------------#
            xfmr["r_wnd"]       = [xdata["r1"], xdata["r2"]]
            xfmr["b_wnd"]       = 0.0
            xfmr["g_wnd"]       = 0.0
            xfmr["r_gnd"]       = [xdata["re1"], xdata["re2"]]
            xfmr["x_gnd"]       = [xdata["xe1"], xdata["xe2"]]
            #-----------------------------------#
            xfmr["tr"]          = calc_xfmr_shift_real(hdata, xdata, h)
            xfmr["ti"]          = calc_xfmr_shift_imaginary(hdata, xdata, h)
            #-----------------------------------#
            xfmr["cxmfr"]       = nothing
            xfmr["cxmfi"]       = nothing
            #-----------------------------------#
            xfmr["i_base_ka"]   = calc_xfmr_current_base(hdata, xdata)
            xfmr["i_fund_magn"] = [0.0, 0.0]
            xfmr["i_rms_max"]   = calc_xfmr_current_rms_max(hdata, xdata)
        else
            xfmr["x_core"]      = [xdata["xsc"]] .* h
            xfmr["b_core"]      = 0.0
            xfmr["g_core"]      = xdata["gsh"] / sqrt(h)
            #-----------------------------------#
            xfmr["r_wnd"]       = [xdata["r1"], xdata["r2"]] .* h
            xfmr["b_wnd"]       = 0.0
            xfmr["g_wnd"]       = 0.0
            xfmr["r_gnd"]       = [xdata["re1"], xdata["re2"]] .* sqrt(h)
            xfmr["x_gnd"]       = [xdata["xe1"], xdata["xe2"]] .* h
            #-----------------------------------#
            xfmr["tr"]          = calc_xfmr_shift_real(hdata, xdata, h)
            xfmr["ti"]          = calc_xfmr_shift_imaginary(hdata, xdata, h)
            #-----------------------------------#
            xfmr["cxmfr"]       = nothing
            xfmr["cxmfi"]       = nothing
    end end
    
    !isempty(xfmr_magn) ? sample_magnetizing_current(hdata, xfmr_magn) : ~ ;
end

# variables ####################################################################
""
function variable_xfmr_voltage(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_xfmr_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_xfmr_voltage_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_voltage_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_xfmr_voltage_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_xfmr_voltage_magnitude_limits(pm, nw)

    vxr     = _PMs.var(pm, nw)[:vxr] =
                JuMP.@variable( pm.model, 
                                [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], 
                                base_name="$(nw)_vxr",
                                start=v_lim[x][i])

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            JuMP.set_lower_bound(vxr[(x,i,j)], -v_lim[x][i])
            JuMP.set_upper_bound(vxr[(x,i,j)],  v_lim[x][i])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :vxr_fr, :vxr_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), vxr) # this report does not make sense
end
""
function variable_xfmr_voltage_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_xfmr_voltage_magnitude_limits(pm, nw)

    vxi     = _PMs.var(pm, nw)[:vxi] = 
                JuMP.@variable( pm.model,
                                [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], 
                                base_name="$(nw)_vxi",
                                start=0.0)

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            JuMP.set_lower_bound(vxi[(x,i,j)], -v_lim[x][i])
            JuMP.set_upper_bound(vxi[(x,i,j)],  v_lim[x][i])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :vxi_fr, :vxi_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), vxi)
end
""
function variable_xfmr_voltage_excitation_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_xfmr_voltage_magnitude_limits(pm, nw)

    exr     = _PMs.var(pm, nw)[:exr] = 
                JuMP.@variable( pm.model,
                                [x in _PMs.ids(pm, nw, :xfmr)], 
                                base_name="$(nw)_exr",
                                start = maximum(values(v_lim[x])))  

    if bounded
        for x in _PMs.ids(pm, nw, :xfmr)
            JuMP.set_lower_bound(exr[x], -maximum(values(v_lim[x])))
            JuMP.set_upper_bound(exr[x],  maximum(values(v_lim[x])))
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :exr, _PMs.ids(pm, nw, :xfmr), exr)
end
function variable_xfmr_voltage_excitation_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_xfmr_voltage_magnitude_limits(pm, nw)

    exi     = _PMs.var(pm, nw)[:exi] = 
                JuMP.@variable( pm.model,
                                [x in _PMs.ids(pm, nw, :xfmr)], 
                                base_name="$(nw)_exi",
                                start=0.0)

    if bounded
        for x in _PMs.ids(pm, nw, :xfmr)
            JuMP.set_lower_bound(exi[x], -maximum(values(v_lim[x])))
            JuMP.set_upper_bound(exi[x],  maximum(values(v_lim[x])))
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :exi, _PMs.ids(pm, nw, :xfmr), exi)
end

""
function variable_xfmr_current(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_xfmr_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_xfmr_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_xfmr_current_magnetizing_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_current_magnetizing_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_xfmr_current_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_xfmr_current_magnitude_limits(pm, nw)
    
    cxr     = _PMs.var(pm, nw)[:cxr] = 
                JuMP.@variable( pm.model,
                                [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], 
                                base_name="$(nw)_cxr",
                                start=0.0)

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            JuMP.set_lower_bound(cxr[(x,i,j)], -c_lim[x][i])
            JuMP.set_upper_bound(cxr[(x,i,j)],  c_lim[x][i])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :cxr_fr, :cxr_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), cxr)
end
""
function variable_xfmr_current_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_xfmr_current_magnitude_limits(pm, nw)
    
    cxi     = _PMs.var(pm, nw)[:cxi] = 
                JuMP.@variable( pm.model,
                                [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], 
                                base_name="$(nw)_cxi",
                                start=0.0)

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            JuMP.set_lower_bound(cxi[(x,i,j)], -c_lim[x][i])
            JuMP.set_upper_bound(cxi[(x,i,j)],  c_lim[x][i])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :cxi_fr, :cxi_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), cxi)
end
""
function variable_xfmr_current_series_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_xfmr_current_magnitude_limits(pm, nw)
    
    cxsr    = _PMs.var(pm, nw)[:cxsr] = 
                JuMP.@variable( pm.model,
                                [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], 
                                base_name="$(nw)_cxsr",
                                start=0.0)

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            JuMP.set_lower_bound(cxsr[(x,i,j)], -c_lim[x][i])
            JuMP.set_upper_bound(cxsr[(x,i,j)],  c_lim[x][i])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :cxsr_fr, :cxsr_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), cxsr)
end
""
function variable_xfmr_current_series_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_xfmr_current_magnitude_limits(pm, nw)
    
    cxsi    = _PMs.var(pm, nw)[:cxsi] = 
                JuMP.@variable( pm.model,
                                [(x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], 
                                base_name="$(nw)_cxsi",
                                start=0.0)

    if bounded
        for (x,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)
            JuMP.set_lower_bound(cxsi[(x,i,j)], -c_lim[x][i])
            JuMP.set_upper_bound(cxsi[(x,i,j)],  c_lim[x][i])
        end
    end

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :cxsi_fr, :cxsi_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), cxsi)
end
""
function variable_xfmr_current_magnetizing_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_xfmr_current_magnitude_limits(pm, nw)
    
    cxmr    = _PMs.var(pm, nw)[:cxmr] = 
                JuMP.@variable( pm.model,
                                [x in _PMs.ids(pm, nw, :xfmr)], 
                                base_name="$(nw)_cxmr",
                                start=0.0)

    if bounded
        for x in _PMs.ids(pm, nw, :xfmr)
            JuMP.set_lower_bound(cxmr[x], -maximum(values(c_lim[x])))
            JuMP.set_upper_bound(cxmr[x],  maximum(values(c_lim[x])))
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :cxmr, _PMs.ids(pm, nw, :xfmr), cxmr)
end
""
function variable_xfmr_current_magnetizing_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_xfmr_current_magnitude_limits(pm, nw)
    
    cxmi    = _PMs.var(pm, nw)[:cxmi] = 
                JuMP.@variable( pm.model,
                                [x in _PMs.ids(pm, nw, :xfmr)], 
                                base_name="$(nw)_cxmi",
                                start=0.0)

    if bounded
        for x in _PMs.ids(pm, nw, :xfmr)
            JuMP.set_lower_bound(cxmi[x], -maximum(values(c_lim[x])))
            JuMP.set_upper_bound(cxmi[x],  maximum(values(c_lim[x])))
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :cxmi, _PMs.ids(pm, nw, :xfmr), cxmi)
end

# constraints ##################################################################
## core constraints ############################################################
### core voltage drop constraint ###############################################
""
function constraint_xfmr_core_voltage_drop(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx     = (x, _PMs.ref(pm, nw, :xfmr, x, "bus")...)
    
    x_core  = _PMs.ref(pm, nw, :xfmr, x, "x_core")
    
    constraint_xfmr_core_voltage_drop(pm, nw, idx, x_core)
end
""
function constraint_xfmr_core_voltage_drop(pm::HarmonicPowerModel, n::Int, idx, x_core)
    exr     = _PMs.var(pm, n, :exr, idx[1])
    exi     = _PMs.var(pm, n, :exi, idx[1])

    vxr     = _PMs.var(pm, n, :vxr, idx)
    vxi     = _PMs.var(pm, n, :vxi, idx)

    cxsr    = _PMs.var(pm, n, :cxsr, idx)
    cxsi    = _PMs.var(pm, n, :cxsi, idx)
    
    JuMP.@constraint(pm.model, vxr == exr - x_core * cxsi)
    JuMP.@constraint(pm.model, vxi == exi + x_core * cxsr)
end

### core current balance constraint ############################################
""
function constraint_xfmr_core_current_balance(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx_fr  = (x, _PMs.ref(pm, nw, :xfmr, x, "bus")...)
    idx_to  = (x, reverse(_PMs.ref(pm, nw, :xfmr, x, "bus"))...)
    
    tr      = _PMs.ref(pm, nw, :xfmr, x, "tr")
    ti      = _PMs.ref(pm, nw, :xfmr, x, "ti")

    g_core  = _PMs.ref(pm, nw, :xfmr, x, "g_core")

    constraint_xfmr_core_current_balance(pm, nw, x, idx_fr, idx_to, tr, ti, g_core)
end
"""
first principles: conj(tₓᵢⱼ) * iₓᵢⱼ + iₓⱼᵢ = 0
conj(tₓᵢⱼₕ) * (iˢₓᵢⱼₕ - iᵐₓₕ - eₓₕ) + iₓⱼᵢₕ = 0
(tʳₓᵢⱼₕ - j tⁱₓᵢⱼₕ) * (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ + j (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ)) + iˢ⁻ʳₓⱼᵢₕ + j iˢ⁻ⁱₓⱼᵢₕ = 0
tʳₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ) + j tʳₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ) - j tⁱₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ) - j² tⁱₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ) + iˢ⁻ʳₓⱼᵢₕ + j iˢ⁻ⁱₓⱼᵢₕ = 0
tʳₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ) + j tʳₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ) - j tⁱₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ) + tⁱₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ) + iˢ⁻ʳₓⱼᵢₕ + j iˢ⁻ⁱₓⱼᵢₕ = 0

Re: tʳₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ) + tⁱₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ) + iˢ⁻ʳₓⱼᵢₕ = 0
Im: tʳₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ -  gˢʰₓₕ * eⁱₓₕ) - tⁱₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ -  gˢʰₓₕ * eʳₓₕ) + iˢ⁻ⁱₓⱼᵢₕ = 0
"""
function constraint_xfmr_core_current_balance(pm::HarmonicPowerModel, n::Int, x, idx_fr, idx_to, tr, ti, g_core)
    cxsr_fr = _PMs.var(pm, n, :cxsr, idx_fr)
    cxsi_fr = _PMs.var(pm, n, :cxsi, idx_fr)

    cxsr_to = _PMs.var(pm, n, :cxsr, idx_to)
    cxsi_to = _PMs.var(pm, n, :cxsi, idx_to)

    cxmr    = _PMs.var(pm, n, :cxmr, x)
    cxmi    = _PMs.var(pm, n, :cxmi, x)

    exr     = _PMs.var(pm, n, :exr, x)
    exi     = _PMs.var(pm, n, :exi, x)

    JuMP.@constraint(pm.model,    cxsr_to
                                + tr * (cxsr_fr - cxmr - g_core * exr)
                                + ti * (cxsi_fr - cxmi - g_core * exi)
                                    == 
                                0.0)
    JuMP.@constraint(pm.model,    cxsi_to
                                + tr * (cxsi_fr - cxmi - g_core * exi)
                                - ti * (cxsr_fr - cxmr - g_core * exr)
                                    == 
                                0.0)
end

### core voltage phase shift constraint ########################################
"" 
function constraint_xfmr_core_voltage_phase_shift(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx = (x, reverse(_PMs.ref(pm, nw, :xfmr, x, "bus"))...)
    
    tr  = _PMs.ref(pm, nw, :xfmr, x, "tr")
    ti  = _PMs.ref(pm, nw, :xfmr, x, "ti")
    
    constraint_xfmr_core_voltage_phase_shift(pm, nw, idx, tr, ti)
end 
"""
first principles: uᵢ = tₓᵢⱼ * uⱼ
eₓₕ = tₓᵢⱼₕ * vₓⱼᵢₕ
eʳₓₕ + j eⁱₓₕ = (tʳₓᵢⱼₕ + j tⁱₓᵢⱼₕ) * (vʳₓⱼᵢₕ + j vⁱₓⱼᵢₕ)
eʳₓₕ + j eⁱₓₕ = tʳₓᵢⱼₕ vʳₓⱼᵢₕ + j tʳₓᵢⱼₕ vⁱₓⱼᵢₕ + j tⁱₓᵢⱼₕ vʳₓⱼᵢₕ + j² tⁱₓᵢⱼₕ vⁱₓⱼᵢₕ
eʳₓₕ + j eⁱₓₕ = tʳₓᵢⱼₕ vʳₓⱼᵢₕ + j tʳₓᵢⱼₕ vⁱₓⱼᵢₕ + j tⁱₓᵢⱼₕ vʳₓⱼᵢₕ - tⁱₓᵢⱼₕ vⁱₓⱼᵢₕ

Re: eʳₓₕ = tʳₓᵢⱼₕ vʳₓⱼᵢₕ - tⁱₓᵢⱼₕ vⁱₓⱼᵢₕ
Im: eⁱₓₕ = tʳₓᵢⱼₕ vⁱₓⱼᵢₕ + tⁱₓᵢⱼₕ vʳₓⱼᵢₕ
"""
function constraint_xfmr_core_voltage_phase_shift(pm::HarmonicPowerModel, n::Int, idx, tr, ti)
    exr = _PMs.var(pm, n, :exr, idx[1])
    exi = _PMs.var(pm, n, :exi, idx[1])

    vxr = _PMs.var(pm, n, :vxr, idx)
    vxi = _PMs.var(pm, n, :vxi, idx)

    JuMP.@constraint(pm.model, exr == tr * vxr - ti * vxi)
    JuMP.@constraint(pm.model, exi == tr * vxi + ti * vxr)
end

### core magnitization constraint ############################################## 
""
function constraint_xfmr_core_magnetization(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    if _PMs.ref(pm, fundamental(pm), :xfmr, x, "linear_magn")
        b_core  = _PMs.ref(pm, nw, :xfmr, x, "b_core")

        constraint_xfmr_core_magnetization(pm, nw, x, b_core)
    else
        if nw ∉ _PMs.ref(pm, nw, :xfmr, x, "Hᴵ")
            b_core  = _PMs.ref(pm, nw, :xfmr, x, "b_core")

            constraint_xfmr_core_magnetization(pm, nw, x, b_core)
        else
            cxmfr = _PMs.ref(pm, nw, :xfmr, x, "cxmfr")
            cxmfi = _PMs.ref(pm, nw, :xfmr, x, "cxmfi")

            constraint_xfmr_core_magnetization(pm, nw, x, cxmfr, cxmfi)
end end end
""
function constraint_xfmr_core_magnetization(pm::HarmonicPowerModel, n::Int, x, b_core)
    exr = _PMs.var(pm, n, :exr, x)
    exi = _PMs.var(pm, n, :exi, x)

    cmrx = _PMs.var(pm, n, :cmrx, x)
    cmix = _PMs.var(pm, n, :cmix, x)

    JuMP.@constraint(pm.model, cmrx == -b_core * exi)
    JuMP.@constraint(pm.model, cmix ==  b_core * exr)
end
""
function constraint_xfmr_core_magnetization(pm::HarmonicPowerModel, n::Int, x, cxmfr, cxmfi)
    cxmr    = _PMs.var(pm, n, :cmrx, x)
    cxmi    = _PMs.var(pm, n, :cmix, x)

    ex      = reduce(vcat,[[_PMs.var(pm, nw, :exr, x), _PMs.var(pm, nw, :exi, x)] 
                            for nw in _PMs.ref(pm, n, :xfmr, x, "Hᴱ")])

    sym_exr = Symbol("exc_re_", n, "_", x)
    sym_exi = Symbol("exc_im_", n, "_", x)

    JuMP.register(pm.model, sym_exr, length(ex), cxmfr; autodiff=true)
    JuMP.register(pm.model, sym_exi, length(ex), cxmfi; autodiff=true)

    JuMP.add_nonlinear_constraint(pm.model, :($(cxmr) == $(sym_exr)($(ex...))))
    JuMP.add_nonlinear_constraint(pm.model, :($(cxmi) == $(sym_exi)($(ex...))))
end

## winding constraints #########################################################
### winding current balance constraint #########################################
""
function constraint_xfmr_winding_current_balance(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx = [arc for arc in _PMs.ref(pm, fundamental(pm), :bus_arcs_xfmr) if arc[1] == x]

    cnf = _PMs.ref(pm, nw, :xfmr, x, "cnf")

    r   = _PMs.ref(pm, nw, :xfmr, x, "r_wnd")
    b   = _PMs.ref(pm, nw, :xfmr, x, "b_wnd")
    g   = _PMs.ref(pm, nw, :xfmr, x, "g_wnd")

    for wnd in 1:_PMs.ref(pm, fundamental(pm), :xfmr, x, "Nw")
        constraint_xfmr_winding_current_balance(pm, nw, idx[wnd], cnf[wnd], r[wnd], b[wnd], g[wnd])
end end
""
function constraint_xfmr_winding_current_balance(pm::HarmonicPowerModel, n::Int, idx, cnf, r_wnd, b_wnd, g_wnd)
    vxr     = _PMs.var(pm, n, :vxr, idx)
    vxi     = _PMs.var(pm, n, :vxi, idx)
    
    cxr     = _PMs.var(pm, n, :cxr, idx)
    cxi     = _PMs.var(pm, n, :cxi, idx)

    cxsr    = _PMs.var(pm, n, :cxsr, idx)
    cxsi    = _PMs.var(pm, n, :cxsi, idx)

    JuMP.@constraint(pm.model, cxr == cxsr - g_wnd * vxr + b_wnd * vxi) # to be checked, is all the same equation with different parameter values, so do this appriory
    JuMP.@constraint(pm.model, cxi == cxsi - g_wnd * vxi - b_wnd * vxr)
end

### winding voltage drop constraint ############################################
""
function constraint_xfmr_winding_voltage_drop(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx     = [arc for arc in _PMs.ref(pm, fundamental(pm), :bus_arcs_xfmr) if arc[1] == x]

    gnd     = _PMs.ref(pm, nw, :xfmr, x, "gnd")

    r_wnd   = _PMs.ref(pm, nw, :xfmr, x, "r_wnd")
    r_gnd   = _PMs.ref(pm, nw, :xfmr, x, "r_gnd")
    x_gnd   = _PMs.ref(pm, nw, :xfmr, x, "x_gnd")

    for wnd in 1:_PMs.ref(pm, fundamental(pm), :xmfr, nw, "Nw")
        constraint_xfmr_winding_voltage_drop(pm, nw, idx[wnd], gnd[w], r_wnd[wnd], r_gnd[wnd], x_gnd[wnd])
end end
""
function constraint_xfmr_winding_voltage_drop(pm::HarmonicPowerModel, n::Int, idx, gnd, r_wnd, r_gnd, x_gnd)
    vbr = _PMs.var(pm, n, :vbr, idx[2])
    vbi = _PMs.var(pm, n, :vbi, idx[2])

    vxr = _PMs.var(pm, n, :vxr, idx)
    vxi = _PMs.var(pm, n, :vxi, idx)

    cxr = _PMs.var(pm, n, :cxr, idx)
    cxi = _PMs.var(pm, n, :cxi, idx)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(n)
        JuMP.@constraint(pm.model, vxr == vbr - r_wnd * cxr)
        JuMP.@constraint(pm.model, vxi == vbi - r_wnd * cxi)
    end

    # h ∈ 𝓗⁰, gnd == true -> cnf ∈ {Ye, Ze}
    if is_zero_sequence(n) && gnd
        JuMP.@constraint(pm.model, vxr == vbr - (r_wnd + 3r_gnd) * cxr + 3x_gnd * cxi)
        JuMP.@constraint(pm.model, vxi == vbi - (r_wnd + 3r_gnd) * cxi - 3x_gnd * cxr)
end end

### winding zero sequence current blocking constraint ##########################
""
function constraint_xfmr_winding_zero_seq_current_blocking(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx = [arc for arc in _PMs.ref(pm, fundamental(pm), :bus_arcs_xfmr) if arc[1] == x]

    gnd = _PMs.ref(pm, nw, :xfmr, x, "gnd")

    for wnd in 1:_PMs.ref(pm, fundamental(pm), :xmfr, nw, "Nw")
        constraint_xfmr_winding_zero_seq_current_blocking(pm, nw, idx[wnd], gnd[w])
    end
end
""
function constraint_xfmr_winding_zero_seq_current_blocking(pm::HarmonicPowerModel, n::Int, idx, gnd)
    cxr = _PMs.var(pm, n, :cxr, idx)
    cxi = _PMs.var(pm, n, :cxi, idx)

    # h ∈ 𝓗⁰, gnd == false -> cnf ∈ {D, Y, Z}
    if is_zero_sequence(n) && !gnd
        JuMP.@constraint(pm.model, cxr == 0)
        JuMP.@constraint(pm.model, cxi == 0)
end end

### winding root-mean-square current limit #####################################
""
function constraint_xfmr_winding_current_rms_limit(pm::HarmonicPowerModel, x::Int)
    idx         = [arc for arc in _PMs.ref(pm, fundamental(pm), :bus_arcs_xfmr) if arc[1] == x]

    i_rms_max   = _PMs.ref(pm, fundamental(pm), :xfmr, x, "i_rms_max")
    i_fund_magn = _PMs.ref(pm, fundamental(pm), :xfmr, x, "i_fund_magn")

    for wnd in 1:_PMs.ref(pm, fundamental(pm), :xfmr, x, "Nw")
        constraint_xfmr_winding_current_rms_limit(pm, idx[wnd], i_rms_max[wnd], i_fund_magn[wnd])
end end
""
function constraint_xfmr_winding_current_rms_limit(pm::HarmonicPowerModel, idx, i_rms_max, i_fund_magn)
    cxr =  [_PMs.var(pm, n, :cxr, idx) for n in sorted_nw_ids(pm)]
    cxi =  [_PMs.var(pm, n, :cxi, idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(cxr.^2 + cxi.^2) <= i_rms_max^2)
end
""
function constraint_xfmr_winding_current_rms_limit(pm::dHHCPowerModel, idx, i_rms_max, i_fund_magn)
    cxr =  [_PMs.var(pm, n, :cxr, idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cxi =  [_PMs.var(pm, n, :cxi, idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(i_rms_max^2 - i_fund_magn^2); vcat(cxr, cxi)] in JuMP.SecondOrderCone())
end