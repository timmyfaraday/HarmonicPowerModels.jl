################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
################################################################################

# ref bus
""
function constraint_voltage_ref_bus(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=fundamental(pm))
    if nw == 1
        vref = 1.0
    else
        if hasref(pm, _PMs.pm_it_sym, nw, :bus, i, "ihdmax")
            vref = _PMs.ref(pm, nw, :bus, i, "ihdmax")
        else
            vref = 0.0 
    end end

    constraint_voltage_ref_bus(pm, nw, i, vref)
end

# bus
""
function constraint_voltage_rms_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=fundamental(pm))
    vminrms = _PMs.ref(pm, nw, :bus, i, "vminrms")
    vmaxrms = _PMs.ref(pm, nw, :bus, i, "vmaxrms")

    constraint_voltage_rms_limit(pm, i, vminrms, vmaxrms)
end
""
function constraint_voltage_rms_limit(pm::dHHC_SOC, i::Int; nw::Int=fundamental(pm))
    vmaxrms = _PMs.ref(pm, nw, :bus, i, "vmaxrms")
    vmfund  = _PMs.ref(pm, fundamental(pm), :bus, i, "vm")

    constraint_voltage_rms_limit(pm, i, vmaxrms, vmfund)
end
""
function constraint_voltage_thd_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=fundamental(pm))
    thdmax = _PMs.ref(pm, nw, :bus, i, "thdmax")
   
    constraint_voltage_thd_limit(pm, i, thdmax)
end
""
function constraint_voltage_thd_limit(pm::dHHC_SOC, i::Int; nw::Int=fundamental(pm))
    thdmax = _PMs.ref(pm, nw, :bus, i, "thdmax")
    vmfund = _PMs.ref(pm, fundamental(pm), :bus, i, "vm")
   
    constraint_voltage_thd_limit(pm, i, thdmax, vmfund)
end
""
function constraint_voltage_ihd_limit(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=fundamental(pm))
    ihdmax = _PMs.ref(pm, nw, :bus, i, "ihdmax")

    if nw ≠ fundamental(pm)
        constraint_voltage_ihd_limit(pm, nw, i, ihdmax)
    end
end
""
function constraint_voltage_ihd_limit(pm::dHHC_SOC, i::Int; nw::Int=fundamental(pm))
    ihdmax = _PMs.ref(pm, nw, :bus, i, "ihdmax")
    vmfund = _PMs.ref(pm, fundamental(pm), :bus, i, "vm")

    if nw ≠ fundamental(pm)
        constraint_voltage_ihd_limit(pm, nw, i, ihdmax, vmfund)
    end
end
""
function constraint_current_balance(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=fundamental(pm))
    bus_arcs      = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_xfmr = _PMs.ref(pm, nw, :bus_arcs_xfmr, i)

    bus_filters   = _PMs.ref(pm, nw, :bus_filters, i)
    bus_gens      = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads     = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts    = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_gs  = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs  = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i,   bus_arcs, bus_arcs_xfmr, 
                                            bus_filters, bus_gens, bus_loads, 
                                            bus_gs, bus_bs)
end

# branch
""
function constraint_current_rms_limit(pm::_PMs.AbstractPowerModel, b::Int; nw::Int=fundamental(pm))
    branch = _PMs.ref(pm, nw, :branch, b)
    f_bus, t_bus = branch["f_bus"], branch["t_bus"]
    f_idx, t_idx = (b, f_bus, t_bus), (b, t_bus, f_bus)

    c_rating = branch["c_rating"]

    constraint_current_rms_limit(pm, f_idx, t_idx, c_rating)
end
""
function constraint_current_rms_limit(pm::dHHC_SOC, b::Int; nw::Int=fundamental(pm))
    branch = _PMs.ref(pm, nw, :branch, b)
    f_bus, t_bus = branch["f_bus"], branch["t_bus"]
    f_idx, t_idx = (b, f_bus, t_bus), (b, t_bus, f_bus)

    c_rating = branch["c_rating"]
    cm_fund_fr = branch["cm_fr"]
    cm_fund_to = branch["cm_to"]

    constraint_current_rms_limit(pm, f_idx, t_idx, c_rating, cm_fund_fr, cm_fund_to)
end

# fairness principle
""
function constraint_fairness_principle(pm::_PMs.AbstractPowerModel)
    principle   = pm.data["principle"]

    load_ids    = sort(collect(_PMs.ids(pm, :load, nw=fundamental(pm))))

    constraint_fairness_principle(pm, principle, load_ids)
end

# filter
""
function constraint_active_filter_current(pm::_PMs.AbstractPowerModel, f::Int; nw::Int=fundamental(pm))
    filter  = _PMs.ref(pm, nw, :filter, f)
    bus     = filter["bus"]

    if filter["a/p"] == "a"
        constraint_active_filter_current(pm, nw, f, bus)
    end
end

# load
""
function constraint_load_current(pm::_PMs.AbstractPowerModel, l::Int; nw::Int=fundamental(pm))
    load = _PMs.ref(pm, nw, :load, l)

    i      = load["load_bus"]
    pd, qd = load["pd"], load["qd"]

    angmin = load["reference_harmonic_angle"] - load["harmonic_angle_range"] / 2
    angmax = load["reference_harmonic_angle"] + load["harmonic_angle_range"] / 2

    if nw == 1
        constraint_load_constant_power(pm, nw, l, i, pd, qd)
    else
        constraint_load_current_angle(pm, nw, l, angmin, angmax)
    end  
end
""
function constraint_load_power(pm::_PMs.AbstractPowerModel, l::Int; nw::Int=fundamental(pm))
    load = _PMs.ref(pm, nw, :load, l)

    i       = load["load_bus"]
    pd, qd  = load["pd"], load["qd"]
    mult    = load["multiplier"]

    if nw == 1
        constraint_load_constant_power(pm, nw, l, i, pd, qd)
    else
        constraint_load_constant_current(pm, nw, l, mult)
    end
end

# xfmr
""
function constraint_transformer_core_magnetization(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    xfmr = _PMs.ref(pm, nw, :xfmr, t)
    Hᴵ = haskey(xfmr, "Hᴵ") ? xfmr["Hᴵ"] : Int[] ;

    if nw in Hᴵ
        int_a = _PMs.ref(pm, nw, :xfmr, t, "Im_A")
        int_b = _PMs.ref(pm, nw, :xfmr, t, "Im_B")

        constraint_transformer_core_magnetization(pm, nw, t, int_a, int_b)
    else 
        constraint_transformer_core_magnetization(pm, nw, t)
    end
end
""
function constraint_transformer_core_voltage_drop(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    f_idx = (t,f_bus,t_bus)
    
    xsc = _PMs.ref(pm, nw, :xfmr, t, "xsc")
    
    constraint_transformer_core_voltage_drop(pm, nw, t, f_idx, xsc)
end
"" 
function constraint_transformer_core_voltage_phase_shift(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    t_idx = (t,t_bus,f_bus)
    
    tr = _PMs.ref(pm, nw, :xfmr, t, "tr")
    ti = _PMs.ref(pm, nw, :xfmr, t, "ti")

    gnd1 = _PMs.ref(pm, nw, :xfmr, t, "gnd1")
    gnd2 = _PMs.ref(pm, nw, :xfmr, t, "gnd2")
    
    constraint_transformer_core_voltage_phase_shift(pm, nw, t, t_idx, tr, ti, gnd1, gnd2)
end
""
function constraint_transformer_core_current_balance(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    f_idx = (t,f_bus,t_bus)
    t_idx = (t,t_bus,f_bus)
    
    tr = _PMs.ref(pm, nw, :xfmr, t, "tr")
    ti = _PMs.ref(pm, nw, :xfmr, t, "ti")

    rsh = _PMs.ref(pm, nw, :xfmr, t, "rsh")

    constraint_transformer_core_current_balance(pm, nw, t, f_idx, t_idx, tr, ti, rsh)
end
""
function constraint_transformer_winding_config(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    w_bus = [f_bus, t_bus]
    w_idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]

    r  = [_PMs.ref(pm, nw, :xfmr, t, nk) for nk in ["r1","r2"]]
    re = [_PMs.ref(pm, nw, :xfmr, t, nk) for nk in ["re1","re2"]]
    xe = [_PMs.ref(pm, nw, :xfmr, t, nk) for nk in ["xe1","xe2"]]

    gnd = [_PMs.ref(pm, nw, :xfmr, t, nk) for nk in ["gnd1","gnd2"]]

    for w in 1:2
        constraint_transformer_winding_config(pm, nw, w_bus[w], w_idx[w], r[w], re[w], xe[w], gnd[w])
    end
end
""
function constraint_transformer_winding_current_balance(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    w_idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]
    
    r = [_PMs.ref(pm, nw, :xfmr, t, nk) for nk in ["r1","r2"]]
    b_sh = [0.0,0.0]
    g_sh = [0.0,0.0]

    cnf = [_PMs.ref(pm, nw, :xfmr, t, nk) for nk in ["cnf1","cnf2"]]

    for w in 1:2
        constraint_transformer_winding_current_balance(pm, nw, w_idx[w], r[w], b_sh[w], g_sh[w], cnf[w])
    end
end
""
function constraint_transformer_winding_current_rms_limit(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=fundamental(pm))
    xfmr = _PMs.ref(pm, nw, :xfmr, t)

    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    w_idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]

    c_rating = xfmr["c_rating"]

    for w in 1:2
        constraint_transformer_winding_current_rms_limit(pm, w_idx[w], c_rating)
    end
end
""
function constraint_transformer_winding_current_rms_limit(pm::dHHC_SOC, t::Int; nw::Int=fundamental(pm))
    xfmr = _PMs.ref(pm, nw, :xfmr, t)

    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    w_idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]

    c_rating = xfmr["c_rating"]

    cm_fund = [xfmr["ctm_fr"], xfmr["ctm_to"]]

    for w in 1:2
        constraint_transformer_winding_current_rms_limit(pm, w_idx[w], c_rating, cm_fund[w])
    end
end