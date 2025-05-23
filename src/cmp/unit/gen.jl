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
""
calc_gen_admittance_real(hdata::Dict{String,Any}, gdata::Dict{String,Any}, h) = # to be reviewed by Hakan
    sqrt( (1 + gdata["xr_ratio"]^2) / (gdata["pmax"]^2 + gdata["qmax"]^2)) / sqrt(h) # gdata["rx_ratio"] * 
""
calc_gen_admittance_imaginary(hdata::Dict{String,Any}, gdata::Dict{String,Any}, h) = # to be reviewed by Hakan
    sqrt( (1 + 1/gdata["xr_ratio"]^2) / (gdata["pmax"]^2 + gdata["qmax"]^2))/ h # gdata["xr_ratio"] * 
""
# i_base_ka = s_base_mva / v_base_kv, see Power System Analysis, pg. 26
calc_gen_current_base(hdata::Dict{String,Any}, gdata::Dict{String,Any}) =
    hdata["s_base_mva"] / hdata["nw"]["1"]["bus"][string(gdata["gen_bus"])]["v_base_kv"] 
""
# i_rms_max = S_nom / s_base_mva / sqrt(3) / v_rms_max(bus)
function calc_gen_current_rms_max(hdata::Dict{String,Any}, gdata::Dict{String,Any})
    S_nom       = sqrt(gdata["pmax"]^2 + gdata["qmax"]^2) 
    v_rms_max   = hdata["nw"]["1"]["bus"][string(gdata["gen_bus"])]["v_rms_max"]
     
    return S_nom / (sqrt(3) * v_rms_max)
end
""
collect_gen_current_magnitude_limits(pm::AbstractHarmonicModel, nw::Int) = 
    Dict(g => _PMs.ref(pm, fundamental(pm), :gen, g, "i_rms_max")
            for g in _PMs.ids(pm, nw, :gen))

# parameters ###################################################################
""
function add_gen_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (g, gen) in ntw["gen"]
        h       = parse(Int, nw)
        gdata   = fdata["gen"][g] 
        if nw == "1"
            gen["id"]           = gdata["index"]
            gen["bus"]          = gdata["gen_bus"]
            #-----------------------------------#
            gen["bsc"]          = calc_gen_admittance_imaginary(hdata, gdata, h)
            gen["gsc"]          = calc_gen_admittance_real(hdata, gdata, h)
            #-----------------------------------#
            gen["i_base_ka"]    = calc_gen_current_base(hdata, gdata)
            gen["i_fund_magn"]  = 0.0
            gen["i_rms_max"]    = calc_gen_current_rms_max(hdata, gdata)
            gen["p_fund_min"]   = gdata["pmin"]
            gen["p_fund_max"]   = gdata["pmax"]
            gen["q_fund_min"]   = gdata["qmin"]
            gen["q_fund_max"]   = gdata["qmax"]
        else
            gen["gsc"]          = calc_gen_admittance_real(hdata, gdata, h)
            gen["bsc"]          = calc_gen_admittance_imaginary(hdata, gdata, h)
end end end

# variables ####################################################################
""
function variable_gen_current(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_gen_current_real(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_gen_current_magnitude_limits(pm, nw)

    cgr = _PMs.var(pm, nw)[:cgr] = 
            JuMP.@variable( pm.model,
                            [g in _PMs.ids(pm, nw, :gen)], 
                            base_name="$(nw)_cgr",
                            start=0.0)

    if bounded
        for g in _PMs.ids(pm, nw, :gen)
            JuMP.set_lower_bound(cgr[g], -c_lim[g])
            JuMP.set_upper_bound(cgr[g],  c_lim[g])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :cgr, _PMs.ids(pm, nw, :gen), cgr)
end
""
function variable_gen_current_imaginary(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    c_lim   = collect_gen_current_magnitude_limits(pm, nw)
    
    cgi = _PMs.var(pm, nw)[:cgi] = 
            JuMP.@variable( pm.model,
                            [g in _PMs.ids(pm, nw, :gen)], 
                            base_name="$(nw)_cgi",
                            start=0.0)

    if bounded
        for g in _PMs.ids(pm, nw, :gen)
            JuMP.set_lower_bound(cgi[g], -c_lim[g])
            JuMP.set_upper_bound(cgi[g],  c_lim[g])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :cgi, _PMs.ids(pm, nw, :gen), cgi)
end

# constraints ##################################################################
## generator current constraint ################################################
""
function constraint_gen_current(pm::AbstractHarmonicModel, g::Int; nw::Int=fundamental(pm))
    i   = _PMs.ref(pm, fundamental(pm), :gen, g, "bus")

    gsc = _PMs.ref(pm, nw, :gen, g, "gsc")
    bsc = _PMs.ref(pm, nw, :gen, g, "bsc")

    if nw ≠ fundamental(pm)
        constraint_gen_current(pm, nw, g, i, gsc, bsc)
end end
""
function constraint_gen_current(pm::AbstractHarmonicModel, n::Int, g, i, gsc, bsc)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    cgr = _PMs.var(pm, n, :cgr, g)
    cgi = _PMs.var(pm, n, :cgi, g)

    JuMP.@constraint(pm.model, cgr == gsc * vbr - bsc * vbi)
    JuMP.@constraint(pm.model, cgi == gsc * vbi + bsc * vbr)
end

## generator root-mean-square current limit ####################################
""
function constraint_gen_current_rms_limit(pm::AbstractHarmonicModel, g::Int)
    i_rms_max   = _PMs.ref(pm, fundamental(pm), :gen, g, "i_rms_max")
    i_fund_magn = _PMs.ref(pm, fundamental(pm), :gen, g, "i_fund_magn")

    constraint_gen_current_rms_limit(pm, g, i_rms_max, i_fund_magn)
end
function constraint_gen_current_rms_limit(pm::HarmonicPowerModel, g, i_rms_max, i_fund_magn)
    cgr =  [_PMs.var(pm, n, :cgr, g) for n in sorted_nw_ids(pm)]
    cgi =  [_PMs.var(pm, n, :cgi, g) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(cgr.^2 + cgi.^2) <= i_rms_max^2)
end
""
function constraint_gen_current_rms_limit(pm::dHHCPowerModel, g, i_rms_max, i_fund_magn)
    cgr =  [_PMs.var(pm, n, :cgr, g) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cgi =  [_PMs.var(pm, n, :cgi, g) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(i_rms_max^2 - i_fund_magn^2); vcat(cgr, cgi)] in JuMP.SecondOrderCone())
end

## generator fundamental active power limit ####################################
""
function constraint_gen_power_active_fundamental_limit(pm::AbstractHarmonicModel, g::Int)
    i           = _PMs.ref(pm, fundamental(pm), :gen, g, "bus")

    p_fund_min  = _PMs.ref(pm, fundamental(pm), :gen, g, "p_fund_min")
    p_fund_max  = _PMs.ref(pm, fundamental(pm), :gen, g, "p_fund_max")

    constraint_gen_power_active_fundamental_limit(pm, g, i, p_fund_min, p_fund_max)
end
""
function constraint_gen_power_active_fundamental_limit(pm::AbstractHarmonicModel, g, i, p_fund_min, p_fund_max)
    vbr = _PMs.var(pm, fundamental(pm), :vbr, i)
    vbi = _PMs.var(pm, fundamental(pm), :vbi, i)

    cgr = _PMs.var(pm, fundamental(pm), :cgr, g)
    cgi = _PMs.var(pm, fundamental(pm), :cgi, g)

    JuMP.@constraint(pm.model, p_fund_min <= vbr * cgr  + vbi * cgi)
    JuMP.@constraint(pm.model,               vbr * cgr  + vbi * cgi <= p_fund_max)
end

## generator fundamental reactive power limit ##################################
""
function constraint_gen_power_reactive_fundamental_limit(pm::AbstractHarmonicModel, g::Int)
    i           = _PMs.ref(pm, fundamental(pm), :gen, g, "bus")

    q_fund_min  = _PMs.ref(pm, fundamental(pm), :gen, g, "q_fund_min")
    q_fund_max  = _PMs.ref(pm, fundamental(pm), :gen, g, "q_fund_max")

    constraint_gen_power_reactive_fundamental_limit(pm, g, i, q_fund_min, q_fund_max)
end
""
function constraint_gen_power_reactive_fundamental_limit(pm::AbstractHarmonicModel, g, i, q_fund_min, q_fund_max)
    vbr = _PMs.var(pm, fundamental(pm), :vbr, i)
    vbi = _PMs.var(pm, fundamental(pm), :vbi, i)

    cgr = _PMs.var(pm, fundamental(pm), :cgr, g)
    cgi = _PMs.var(pm, fundamental(pm), :cgi, g)

    JuMP.@constraint(pm.model, q_fund_min <= vbi * cgr  - vbr * cgi)
    JuMP.@constraint(pm.model,               vbi * cgr  - vbr * cgi <= q_fund_max)
end





