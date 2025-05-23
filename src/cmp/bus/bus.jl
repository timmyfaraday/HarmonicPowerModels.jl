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
fundamental_bus_voltage_multiplier(pm::HarmonicPowerModel, i) = 
    _PMs.ref(pm, fundamental(pm), :bus, i, "v_rms_max")
""
fundamental_bus_voltage_multiplier(pm::dHHCPowerModel, i) = 
    _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")
""
bus_voltage_magnitude_limit(pm::AbstractHarmonicModel, nw, i) = 
    nw == fundamental(pm) ? _PMs.ref(pm, nw, :bus, i, "v_rms_max") :
                            _PMs.ref(pm, nw, :bus, i, "v_ihd_max") *
                            fundamental_bus_voltage_multiplier(pm, i)
""
collect_bus_voltage_magnitude_limits(pm::AbstractHarmonicModel, nw::Int) =
    Dict(i => bus_voltage_magnitude_limit(pm, nw, i) for i in _PMs.ids(pm, nw, :bus))

# parameters ###################################################################
""
function add_bus_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (nb, bus) in ntw["bus"]
        h       = parse(Int, nw)
        bdata   = fdata["bus"][nb]
        if nw == "1"
            bus["id"]           = bdata["index"]
            bus["std"]          = bdata["std"]
            bus["type"]         = bdata["bus_type"]
            #-----------------------------------#
            bus["v_base_kv"]    = bdata["base_kv"]
            bus["v_fund_magn"]  = 1.0
            bus["v_rms_min"]    = bdata["vmin"]
            bus["v_rms_max"]    = bdata["vmax"]
            bus["v_thd_max"]    = voltage_thd_limits[bdata["std"]]
        else
            bus["v_ihd_max"]    = voltage_ihd_limits[bdata["std"]][h]
end end end

# variables ####################################################################
""
function variable_bus_voltage(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_bus_voltage_real(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_bus_voltage_magnitude_limits(pm, nw)

    vbr     = _PMs.var(pm, nw)[:vbr] = 
                JuMP.@variable( pm.model,
                                [i in _PMs.ids(pm, nw, :bus)], 
                                base_name="$(nw)_vbr",
                                start=v_lim[i])

    if bounded
        for i in _PMs.ids(pm, nw, :bus)
            JuMP.set_lower_bound(vbr[i], -v_lim[i])
            JuMP.set_upper_bound(vbr[i],  v_lim[i])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vbr, _PMs.ids(pm, nw, :bus), vbr)
end

""
function variable_bus_voltage_imaginary(pm::AbstractHarmonicModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_bus_voltage_magnitude_limits(pm, nw)

    vbi     = _PMs.var(pm, nw)[:vbi] = 
                JuMP.@variable( pm.model,
                                [i in _PMs.ids(pm, nw, :bus)], 
                                base_name="$(nw)_vbi",
                                start=0.0)

    if bounded
        for i in _PMs.ids(pm, nw, :bus)
            JuMP.set_lower_bound(vbi[i], -v_lim[i])
            JuMP.set_upper_bound(vbi[i],  v_lim[i])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vbi, _PMs.ids(pm, nw, :bus), vbi)
end

# constraints ##################################################################
## Kirchhoff's current law #####################################################
""
function constraint_bus_current_balance(pm::AbstractHarmonicModel, i::Int; nw::Int=fundamental(pm))
    bus_arcs_branch = _PMs.ref(pm, nw, :bus_arcs_branch, i)

    bus_wnds_xfmr   = _PMs.ref(pm, nw, :bus_wnds_xfmr, i)

    bus_filter      = _PMs.ref(pm, nw, :bus_filter, i)
    bus_gen         = _PMs.ref(pm, nw, :bus_gen, i)
    bus_hload       = _PMs.ref(pm, nw, :bus_hload, i)
    bus_hsrc        = _PMs.ref(pm, nw, :bus_hsrc, i)
    bus_shunt       = _PMs.ref(pm, nw, :bus_shunt, i)

    constraint_bus_current_balance(pm, nw, bus_arcs_branch, bus_wnds_xfmr, 
                                           bus_filter, bus_gen, bus_hload, bus_hsrc, bus_shunt)
end
""
function constraint_bus_current_balance(pm::AbstractHarmonicModel, n::Int, bus_arcs_branch, bus_wnds_xfmr, bus_filter, bus_gen, bus_hload, bus_hsrc, bus_shunt)
    cbr = _PMs.var(pm, n, :cbr)
    cbi = _PMs.var(pm, n, :cbi)
    cxr = _PMs.var(pm, n, :cxr)
    cxi = _PMs.var(pm, n, :cxi)

    cfr = _PMs.var(pm, n, :cfr)
    cfi = _PMs.var(pm, n, :cfi)
    cgr = _PMs.var(pm, n, :cgr)
    cgi = _PMs.var(pm, n, :cgi)
    clr = _PMs.var(pm, n, :clr)
    cli = _PMs.var(pm, n, :cli)
    crr = _PMs.var(pm, n, :crr)
    cri = _PMs.var(pm, n, :cri)
    csr = _PMs.var(pm, n, :csr)
    csi = _PMs.var(pm, n, :csi)

    JuMP.@constraint(pm.model,    sum(cbr[b] for b in bus_arcs_branch)
                                + sum(cxr[x] for x in bus_wnds_xfmr)
                                ==
                                  sum(cfr[f] for f in bus_filter)
                                + sum(cgr[g] for g in bus_gen) 
                                - sum(clr[l] for l in bus_hload)
                                - sum(crr[r] for r in bus_hsrc)
                                - sum(csr[s] for s in bus_shunt))

    JuMP.@constraint(pm.model,    sum(cbi[b] for b in bus_arcs_branch)
                                + sum(cxi[x] for x in bus_wnds_xfmr)
                                ==
                                  sum(cfi[f] for f in bus_filter)
                                + sum(cgi[g] for g in bus_gen)
                                - sum(cli[l] for l in bus_hload)
                                - sum(cri[r] for r in bus_hsrc)
                                - sum(csi[s] for s in bus_shunt))
end

## individual harmonic voltage distortion limit ################################
""
function constraint_bus_voltage_ihd_limit(pm::AbstractHarmonicModel, i::Int; nw::Int=fundamental(pm))
    if nw ≠ fundamental(pm)
        v_ihd_max   = _PMs.ref(pm, nw, :bus, i, "v_ihd_max")
        v_fund_magn = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")

        println(v_fund_magn)

        constraint_bus_voltage_ihd_limit(pm, nw, i, v_ihd_max, v_fund_magn)
    end
end
""
function constraint_bus_voltage_ihd_limit(pm::HarmonicPowerModel, n::Int, i, v_ihd_max, v_fund_magn)
    vbr = [_PMs.var(pm, 1, :vbr, i), _PMs.var(pm, n, :vbr, i)] 
    vbi = [_PMs.var(pm, 1, :vbi, i), _PMs.var(pm, n, :vbi, i)]

    JuMP.@constraint(pm.model, (vbr[2]^2 + vbi[2]^2) <= v_ihd_max^2 * (vbr[1]^2 + vbi[1]^2))
end
""
function constraint_bus_voltage_ihd_limit(pm::dHHCPowerModel, n::Int, i, v_ihd_max, v_fund_magn)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    JuMP.@constraint(pm.model, [v_ihd_max * v_fund_magn; vcat(vbr, vbi)] in JuMP.SecondOrderCone())
end

## root-mean-square voltage limit ##############################################
""
function constraint_bus_voltage_rms_limit(pm::AbstractHarmonicModel, i::Int)
    v_rms_min   = _PMs.ref(pm, fundamental(pm), :bus, i, "v_rms_min")
    v_rms_max   = _PMs.ref(pm, fundamental(pm), :bus, i, "v_rms_max")
    v_fund_magn = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")

    constraint_bus_voltage_rms_limit(pm, i, v_rms_min, v_rms_max, v_fund_magn)
end
""
function constraint_bus_voltage_rms_limit(pm::HarmonicPowerModel, i, v_rms_min, v_rms_max, v_fund_magn)
    vbr = [_PMs.var(pm, n, :vbr, i) for n in sorted_nw_ids(pm)]
    vbi = [_PMs.var(pm, n, :vbi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, v_rms_min^2 <= sum(vbr.^2 + vbi.^2)                )
    JuMP.@constraint(pm.model,                sum(vbr.^2 + vbi.^2)  <= v_rms_max^2)
end
""
function constraint_bus_voltage_rms_limit(pm::dHHCPowerModel, i, v_rms_min, v_rms_max, v_fund_magn)
    vbr = [_PMs.var(pm, n, :vbr, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    vbi = [_PMs.var(pm, n, :vbi, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    println(v_fund_magn, v_rms_max^2 - v_fund_magn^2)

    JuMP.@constraint(pm.model, [sqrt(v_rms_max^2 - v_fund_magn^2); vcat(vbr, vbi)] in JuMP.SecondOrderCone())
end

## total harmonic voltage distortion limit #####################################
""
function constraint_bus_voltage_thd_limit(pm::AbstractHarmonicModel, i::Int)
    v_thd_max   = _PMs.ref(pm, fundamental(pm), :bus, i, "v_thd_max")
    v_fund_magn = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")
   
    constraint_bus_voltage_thd_limit(pm, i, v_thd_max, v_fund_magn)
end
""
function constraint_bus_voltage_thd_limit(pm::HarmonicPowerModel, i, v_thd_max, v_fund_magn)
    vbr = [_PMs.var(pm, n, :vbr, i) for n in sorted_nw_ids(pm)]
    vbi = [_PMs.var(pm, n, :vbi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(vbr[2:end].^2 + vbi[2:end].^2) <= v_thd_max^2 * (vbr[1]^2 + vbi[1]^2))
end
""
function constraint_bus_voltage_thd_limit(pm::dHHCPowerModel, i, v_thd_max, v_fund_magn)
    vbr = [_PMs.var(pm, n, :vbr, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    vbi = [_PMs.var(pm, n, :vbi, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [v_thd_max * v_fund_magn; vcat(vbr, vbi)] in JuMP.SecondOrderCone())
end