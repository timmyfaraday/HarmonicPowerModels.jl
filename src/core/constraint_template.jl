################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# active filter
""
function constraint_active_filter(pm::_PMs.AbstractPowerModel, g::Int; nw::Int=nw_id_default(pm))
    gen = _PMs.ref(pm, nw, :gen, g)

    if haskey(gen, "isfilter") && gen["isfilter"] == 1
        constraint_active_filter(pm, nw, g)
    end
end

# ref bus
""
function constraint_ref_bus(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=nw_id_default(pm))
    vref = ref(pm, nw, :bus, i, "ihdmax")

    constraint_ref_bus(pm, nw, i, vref)
end

# bus
""
function constraint_voltage_rms_limit(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default(pm))
    vminrms = ref(pm, nw, :bus, i, "vminrms")
    vmaxrms = ref(pm, nw, :bus, i, "vmaxrms")

    constraint_voltage_rms_limit(pm, i, vminrms, vmaxrms)
end
""
function constraint_voltage_rms_limit(pm::QC_DHHC, i::Int; nw::Int=nw_id_default(pm))
    vmaxrms = ref(pm, nw, :bus, i, "vmaxrms")

    constraint_voltage_rms_limit(pm, i, vmaxrms)
end
""
function constraint_voltage_rms_limit(pm::SOC_DHHC, i::Int; nw::Int=nw_id_default(pm))
    vmaxrms = ref(pm, nw, :bus, i, "vmaxrms")

    constraint_voltage_rms_limit(pm, i, vmaxrms)
end
""
function constraint_voltage_thd_limit(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default(pm))
    thdmax = ref(pm, nw, :bus, i, "thdmax")
   
    constraint_voltage_thd_limit(pm, i, thdmax)
end
""
function constraint_voltage_ihd_limit(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default(pm))
    ihdmax = ref(pm, nw, :bus, i, "ihdmax")

    constraint_voltage_ihd_limit(pm, nw, i, ihdmax)
end
""
function constraint_voltage_magnitude_sqr(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default(pm))
    constraint_voltage_magnitude_sqr(pm, nw, i)
end
""
function constraint_current_balance(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=nw_id_default(pm))
    bus_arcs      = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_xfmr = _PMs.ref(pm, nw, :bus_arcs_xfmr, i)
    bus_gens      = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads     = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts    = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_arcs_xfmr, bus_gens, bus_loads, bus_gs, bus_bs)
end

# branch
function constraint_current_rms_limit(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default(pm))
    branch = ref(pm, nw, :branch, b)
    
    f_bus, t_bus = branch["f_bus"], branch["t_bus"]
    f_idx, t_idx = (b, f_bus, t_bus), (b, t_bus, f_bus)

    if haskey(branch, "c_rating_a")
        c_rating = branch["c_rating_a"]

        constraint_current_rms_limit(pm, f_idx, t_idx, c_rating)
    end
end

# load
""
function constraint_load_current(pm::_PMs.AbstractPowerModel, l::Int; nw::Int=nw_id_default(pm))
    load = _PMs.ref(pm, nw, :load, l)

    i       = load["load_bus"]
    pd, qd  = load["pd"], load["qd"]

    angmin = load["reference_harmonic_angle"] - load["harmonic_angle_range"] / 2
    angmax = load["reference_harmonic_angle"] + load["harmonic_angle_range"] / 2

    if nw == 1
        constraint_load_constant_power(pm, nw, l, i, pd, qd)
    else
        # constraint_load_current_fixed_angle(pm, nw, l)
        constraint_load_current_variable_angle(pm, nw, l, angmin, angmax)
        # constraint_load_current_variable_angle_relative(pm, nw, l)
    end  
end
""
function constraint_load_power(pm::_PMs.AbstractPowerModel, l::Int; nw::Int=nw_id_default(pm))
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
function constraint_transformer_core_excitation(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=nw_id_default(pm))
    xfmr = _PMs.ref(pm, nw, :xfmr, t)
    Hᴵ = haskey(xfmr, "Hᴵ") ? xfmr["Hᴵ"] : Int[] ;

    if nw in Hᴵ
        int_a = _PMs.ref(pm, nw, :xfmr, t, "Im_A")
        int_b = _PMs.ref(pm, nw, :xfmr, t, "Im_B")

        constraint_transformer_core_excitation(pm, nw, t, int_a, int_b)
    else 
        constraint_transformer_core_excitation(pm, nw, t)
    end
end
""
function constraint_transformer_core_voltage_drop(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=nw_id_default(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    f_idx = (t,f_bus,t_bus)
    
    xsc = _PMs.ref(pm, nw, :xfmr, t, "xsc")
    
    constraint_transformer_core_voltage_drop(pm, nw, t, f_idx, xsc)
end
"" 
function constraint_transformer_core_voltage_balance(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=nw_id_default(pm))
    f_bus = _PMs.ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = _PMs.ref(pm, nw, :xfmr, t, "t_bus")
    t_idx = (t,t_bus,f_bus)
    
    tr = _PMs.ref(pm, nw, :xfmr, t, "tr")
    ti = _PMs.ref(pm, nw, :xfmr, t, "ti")

    constraint_transformer_core_voltage_balance(pm, nw, t, t_idx, tr, ti)
end
""
function constraint_transformer_core_current_balance(pm::_PMs.AbstractPowerModel, t::Int; nw::Int=nw_id_default(pm))
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    f_idx = (t,f_bus,t_bus)
    t_idx = (t,t_bus,f_bus)
    
    tr = ref(pm, nw, :xfmr, t, "tr")
    ti = ref(pm, nw, :xfmr, t, "ti")

    rsh = ref(pm, nw, :xfmr, t, "rsh")

    constraint_transformer_core_current_balance(pm, nw, t, f_idx, t_idx, tr, ti, rsh)
end
""
function constraint_transformer_winding_config(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default(pm))
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    w_bus = [f_bus, t_bus]
    w_idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]

    r  = [ref(pm, nw, :xfmr, t, nk) for nk in ["r1","r2"]]
    re = [ref(pm, nw, :xfmr, t, nk) for nk in ["re1","re2"]]
    xe = [ref(pm, nw, :xfmr, t, nk) for nk in ["xe1","xe2"]]

    gnd = [ref(pm, nw, :xfmr, t, nk) for nk in ["gnd1","gnd2"]]

    for w in 1:2
        constraint_transformer_winding_config(pm, nw, w_bus[w], w_idx[w], r[w], re[w], xe[w], gnd[w])
    end
end
""
function constraint_transformer_winding_current_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default(pm))
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    w_idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]
    
    r = [ref(pm, nw, :xfmr, t, nk) for nk in ["r1","r2"]]
    b_sh = [0.0,0.0]
    g_sh = [0.0,0.0]

    cnf = [ref(pm, nw, :xfmr, t, nk) for nk in ["cnf1","cnf2"]]

    for w in 1:2
        constraint_transformer_winding_current_balance(pm, nw, w_idx[w], r[w], b_sh[w], g_sh[w], cnf[w])
    end
end

function constraint_load_current_variable_angle_relative(pm::AbstractPowerModel, nw::Int, l::Int)
    load = _PMs.ref(pm, nw, :load, l)
    bus_idx = load["load_bus"]
    c1 = cos(load["reference_harmonic_angle"])
    c2 = sin(load["reference_harmonic_angle"])

    constraint_load_current_variable_angle_relative(pm, nw, l, bus_idx, c1, c2)
end