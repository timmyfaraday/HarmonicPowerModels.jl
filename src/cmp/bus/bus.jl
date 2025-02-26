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
fundamental_voltage_multiplier(pm::HarmonicPowerModel, i) = 
    _PMs.ref(pm, fundamental(pm), :bus, i, "v_rms_max")
fundamental_voltage_multiplier(pm::dHHCPowerModel, i) = 
    _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")

""
collect_bus_voltage_magnitude_limits(pm::HarmonicPowerModel, nw::Int) =
    Dict(i => ifelse(   nw == fundamental(pm),
                        _PMs.ref(pm, nw, :bus, i, "v_rms_max"), 
                        _PMs.ref(pm, nw, :bus, i, "v_ihd_max") * 
                        fundamental_voltage_multiplier(pm, i)
                    ) 
            for i in _PMs.ids(pm, nw, :bus))

# parameters ###################################################################
""
function add_bus_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (nb, bus) in ntw["bus"]
        bdata = fdata["bus"][nb]
        if nw == "1"
            bus = Dict( "id"            => bdata["index"],
                        "std"           => bdata["std"],
                        "type"          => bdata["type"],
                        #-----------------------------------#
                        "v_base_kv"     => bdata["base_kv"],
                        "v_fund_magn"   => 1.0,
                        "v_rms_min"     => bdata["vmin"],
                        "v_rms_max"     => bdata["vmax"],
                        "v_thd_max"     => voltage_thd_limits[bdata["std"]])
        else
            bus = Dict( "v_ihd_max"     => voltage_ihd_limits[bdata["std"]])
end end end

# variables ####################################################################
""
function variable_bus_voltage(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_bus_voltage_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_bus_voltage_magnitude_limits(pm, nw)

    vr      = _PMs.var(pm, nw)[:vr] = 
                JuMP.@variable( pm.model,
                                [i in _PMs.ids(pm, nw, :bus)], 
                                base_name="$(nw)_vr",
                                start=v_lim)

    if bounded
        for i in _PMs.ids(pm, nw, :bus)
            JuMP.set_lower_bound(vr[i], -v_lim[i])
            JuMP.set_upper_bound(vr[i],  v_lim[i])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vr, _PMs.ids(pm, nw, :bus), vr)
end

""
function variable_bus_voltage_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    v_lim   = collect_bus_voltage_magnitude_limits(pm, nw)

    vi      = _PMs.var(pm, nw)[:vi] = 
                JuMP.@variable( pm.model,
                                [i in _PMs.ids(pm, nw, :bus)], 
                                base_name="$(nw)_vi",
                                start=0.0)

    if bounded
        for i in _PMs.ids(pm, nw, :bus)
            JuMP.set_lower_bound(vi[i], -v_lim[i])
            JuMP.set_upper_bound(vi[i],  v_lim[i])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :bus, :vi, _PMs.ids(pm, nw, :bus), vi)
end

# constraints ##################################################################
## Kirchhoff's current law #####################################################
""
function constraint_bus_current_balance(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    bus_arcs      = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_xfmr = _PMs.ref(pm, nw, :bus_arcs_xfmr, i)

    bus_filters   = _PMs.ref(pm, nw, :bus_filters, i)
    bus_gens      = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads     = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts    = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_gs  = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs  = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
    gen_bg  = Dict(k => _PMs.ref(pm, nw, :gen, k, "bg") for k in bus_gens)

    constraint_bus_current_balance(pm, nw, i,   bus_arcs, bus_arcs_xfmr, 
                                            bus_filters, bus_gens, bus_loads, 
                                            bus_gs, bus_bs, gen_bg)
end
""
function constraint_bus_current_balance(pm::HarmonicPowerModel, n::Int, i, bus_arcs, bus_arcs_xfmr, bus_filters, bus_gens, bus_loads, bus_gs, bus_bs, gen_bg)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    cr = _PMs.var(pm, n, :cr)
    ci = _PMs.var(pm, n, :ci)
    crx = _PMs.var(pm, n, :crx)
    cix = _PMs.var(pm, n, :cix)

    crf = _PMs.var(pm, n, :crf)
    cif = _PMs.var(pm, n, :cif)
    # crg = _PMs.var(pm, n, :crg)
    # cig = _PMs.var(pm, n, :cig)
    crd = _PMs.var(pm, n, :crd)
    cid = _PMs.var(pm, n, :cid)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crx[t] for t in bus_arcs_xfmr)
                                ==
                                sum(crf[f] for f in bus_filters)
                                #+ sum(crg[g] for g in bus_gens) # remove later for generators + add sum(gs for gs in values(gen_gs))*vr .....
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr 
                                + sum(bs for bs in values(bus_bs))*vi
                                + sum(bg for bg in values(gen_bg))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cix[t] for t in bus_arcs_xfmr)
                                ==
                                sum(cif[f] for f in bus_filters)
                                # + sum(cig[g] for g in bus_gens) # check later for generators
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi 
                                - sum(bs for bs in values(bus_bs))*vr
                                - sum(bg for bg in values(gen_bg))*vr
                                )
end

## individual harmonic voltage distortion limit ################################
""
function constraint_bus_voltage_ihd_limit(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    if nw ≠ fundamental(pm)
        v_ihd_max   = _PMs.ref(pm, nw, :bus, i, "v_ihd_max")
        v_fund_magn = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")

        constraint_bus_voltage_ihd_limit(pm, nw, i, v_ihd_max, v_fund_magn)
    end
end
""
function constraint_bus_voltage_ihd_limit(pm::HarmonicPowerModel, n::Int, i, v_ihd_max, v_fund_magn)
    vr = [_PMs.var(pm, 1, :vr, i), _PMs.var(pm, n, :vr, i)] 
    vi = [_PMs.var(pm, 1, :vi, i), _PMs.var(pm, n, :vi, i)]

    JuMP.@constraint(pm.model, (vr[2]^2 + vi[2]^2) <= v_ihd_max^2 * (vr[1]^2 + vi[1]^2))
end
""
function constraint_bus_voltage_ihd_limit(pm::dHHCPowerModel, n::Int, i, v_ihd_max, v_fund_magn)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, [v_ihd_max * v_fund_magn; vcat(vr, vi)] in JuMP.SecondOrderCone())
end

## root-mean-square voltage limit ##############################################
""
function constraint_bus_voltage_rms_limit(pm::HarmonicPowerModel, i::Int)
    v_min_rms   = _PMs.ref(pm, fundamental(pm), :bus, i, "v_rms_min")
    v_max_rms   = _PMs.ref(pm, fundamental(pm), :bus, i, "v_rms_max")
    v_fund_magn = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")

    constraint_bus_voltage_rms_limit(pm, i, v_rms_min, v_rms_max, v_fund_magn)
end
""
function constraint_bus_voltage_rms_limit(pm::HarmonicPowerModel, i, v_rms_min, v_rms_max, v_fund_magn)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, v_rms_min^2 <= sum(vr.^2 + vi.^2)                )
    JuMP.@constraint(pm.model,                sum(vr.^2 + vi.^2)  <= v_rms_max^2)
end
""
function constraint_bus_voltage_rms_limit(pm::dHHCPowerModel, i, v_rms_max, v_fund_magn)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(v_rms_max^2 - v_fund_magn^2); vcat(vr, vi)] in JuMP.SecondOrderCone())
end

## total harmonic voltage distortion limit #####################################
""
function constraint_bus_voltage_thd_limit(pm::HarmonicPowerModel, i::Int)
    v_thd_max   = _PMs.ref(pm, fundamental(pm), :bus, i, "v_thd_max")
    v_fund_magn = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_magn")
   
    constraint_bus_voltage_thd_limit(pm, i, v_thd_max, v_fund_magn)
end
""
function constraint_bus_voltage_thd_limit(pm::HarmonicPowerModel, i, v_thd_max, v_fund_magn)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(vr[2:end].^2 + vi[2:end].^2) <= v_thd_max^2 * (vr[1]^2 + vi[1]^2))
end
""
function constraint_bus_voltage_thd_limit(pm::HarmonicPowerModel, i, v_thd_max, v_fund_magn)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [v_thd_max * v_fund_magn; vcat(vr, vi)] in JuMP.SecondOrderCone())
end