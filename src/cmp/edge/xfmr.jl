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
calc_xfmr_current_base(hdata::Dict{String,Any}, bdata::Dict{String,Any}) =
    [hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(nb)]["v_base_kv"] 
        for nb in [bdata["f_bus"], bdata["t_bus"]]]
""
# i_rms_max = S_nom / s_base_mva / sqrt(3) / min(v_rms_max(f_bus), v_rms_max(t_bus))
function calc_xfmr_current_rms_max(hdata::Dict{String,Any}, bdata::Dict{String,Any})
    S_nom       = bdata["rate_a"] 
    s_base_mva  = hdata["s_base_mva"]
    v_rms_max   = [ hdata["nw"]["1"]["bus"][string(bdata["f_bus"])]["v_rms_max"],
                    hdata["nw"]["1"]["bus"][string(bdata["t_bus"])]["v_rms_max"]]
     
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
function add_xfmr_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (nb, xfmr) in ntw["xfmr"]
        h       = parse(Int, nw)
        bdata   = fdata["xfmr"][nb]
        if nw == "1"
            xfmr = Dict("id"            => bdata["index"],
                        "nw"            => 2,
                        "bus"           => [bdata["f_bus"], bdata["t_bus"]],
                        "cnf"           => [],
                        "gnd"           => [],
                        #-----------------------------------#
                        "r"             => bdata["br_r"],
                        "x"             => bdata["br_x"],
                        "g_fr"          => bdata["g_fr"],
                        "b_fr"          => bdata["b_fr"],
                        "g_to"          => bdata["g_fr"],
                        "b_to"          => bdata["b_to"],
                        #-----------------------------------#
                        "i_base_ka"     => calc_xfmr_current_base(hdata, bdata),
                        "i_fund_magn"   => [0.0, 0.0],
                        "i_rms_max"     => calc_xfmr_current_rms_max(hdata, bdata))
        else
            xfmr = Dict("r"             => bdata["br_r"] * sqrt(h),
                        "x"             => bdata["br_x"] * h,
                        "g_fr"          => bdata["g_fr"] / sqrt(h),
                        "b_fr"          => bdata["b_fr"] / h,
                        "g_to"          => bdata["g_fr"] / sqrt(h),
                        "b_to"          => bdata["b_to"] / h)
end end end

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
## core constraint #############################################################
""
function constraint_xfmr_core_voltage_drop(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx = (x, _PMs.ref(pm, nw, :xfmr, x, "bus")...)
    
    xsc = _PMs.ref(pm, nw, :xfmr, x, "x")
    
    constraint_xfmr_core_voltage_drop(pm, nw, idx, xsc)
end
""
function constraint_xfmr_core_voltage_drop(pm::HarmonicPowerModel, n::Int, idx, xsc)
    exr     = _PMs.var(pm, n, :exr, idx[1])
    exi     = _PMs.var(pm, n, :exi, idx[1])

    vxr     = _PMs.var(pm, n, :vxr, idx)
    vxi     = _PMs.var(pm, n, :vxi, idx)

    cxsr    = _PMs.var(pm, n, :cxsr, idx)
    cxsi    = _PMs.var(pm, n, :cxsi, idx)
    
    JuMP.@constraint(pm.model, vxr == exr - x * cxsi)
    JuMP.@constraint(pm.model, vxi == exi + x * cxsr)
end
""
function constraint_xfmr_core_current_balance(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    idx_fr  = (x, _PMs.ref(pm, nw, :xfmr, x, "bus")...)
    idx_to  = (x, reverse(_PMs.ref(pm, nw, :xfmr, x, "bus"))...)
    
    tr      = _PMs.ref(pm, nw, :xfmr, x, "tr")
    ti      = _PMs.ref(pm, nw, :xfmr, x, "ti")

    g       = _PMs.ref(pm, nw, :xfmr, x, "g")

    constraint_xfmr_core_current_balance(pm, nw, x, idx_fr, idx_to, tr, ti, g)
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
function constraint_xfmr_core_current_balance(pm::HarmonicPowerModel, n::Int, x, idx_fr, idx_to, tr, ti, g)
    cxsr_fr = _PMs.var(pm, n, :cxsr, idx_fr)
    cxsi_fr = _PMs.var(pm, n, :cxsi, idx_fr)

    cxsr_to = _PMs.var(pm, n, :cxsr, idx_to)
    cxsi_to = _PMs.var(pm, n, :cxsi, idx_to)

    cxmr    = _PMs.var(pm, n, :cxmr, x)
    cxmi    = _PMs.var(pm, n, :cxmi, x)

    exr     = _PMs.var(pm, n, :exr, x)
    exi     = _PMs.var(pm, n, :exi, x)

    JuMP.@constraint(pm.model,    cxsr_to
                                + tr * (cxsr_fr - cxmr - g * exr)
                                + ti * (cxsi_fr - cxmi - g * exi)
                                    == 
                                0.0)
    JuMP.@constraint(pm.model,    cxsi_to
                                + tr * (cxsi_fr - cxmi - g * exi)
                                - ti * (cxsr_fr - cxmr - g * exr)
                                    == 
                                0.0)
end
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
""
function constraint_xfmr_core_magnetization(pm::HarmonicPowerModel, x::Int; nw::Int=fundamental(pm))
    xfmr    = _PMs.ref(pm, nw, :xfmr, x)
    Hᴵ      = haskey(xfmr, "Hᴵ") ? xfmr["Hᴵ"] : Int[] ;

    if nw in Hᴵ
        int_a = _PMs.ref(pm, nw, :xfmr, x, "Im_A")
        int_b = _PMs.ref(pm, nw, :xfmr, x, "Im_B")

        constraint_xfmr_core_magnetization(pm, nw, x, int_a, int_b)
    else 
        constraint_xfmr_core_magnetization(pm, nw, x)
    end
end
""
function constraint_xfmr_core_magnetization(pm::_PMs.AbstractIVRModel, n::Int, x)
    cmrx = _PMs.var(pm, n, :cmrx, x)
    cmix = _PMs.var(pm, n, :cmix, x)

    JuMP.@constraint(pm.model, cmrx == 0.0)
    JuMP.@constraint(pm.model, cmix == 0.0)
end
""
function constraint_xfmr_core_magnetization(pm::_PMs.AbstractIVRModel, n::Int, x, int_a, int_b)
    cmrx = _PMs.var(pm, n, :cmrx, x)
    cmix = _PMs.var(pm, n, :cmix, x)

    et = reduce(vcat,[[_PMs.var(pm, nw, :erx, x), _PMs.var(pm, nw, :eix, x)] 
                                for nw in _PMs.ref(pm, n, :xfmr, x, "Hᴱ")])

    sym_exc_a = Symbol("exc_a_", n, "_", x)
    sym_exc_b = Symbol("exc_b_", n, "_", x)

    JuMP.register(pm.model, sym_exc_a, length(et), int_a; autodiff=true)
    JuMP.register(pm.model, sym_exc_b, length(et), int_b; autodiff=true)

    JuMP.add_nonlinear_constraint(pm.model, :($(cmrx) == $(sym_exc_a)($(et...))))
    JuMP.add_nonlinear_constraint(pm.model, :($(cmix) == $(sym_exc_b)($(et...))))
end