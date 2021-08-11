
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus_arcs      = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc   = _PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_xfmr = _PMs.ref(pm, nw, :bus_arcs_xfmr, i)
    bus_gens      = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads     = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts    = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_xfmr, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end

""
function constraint_transformer_core_excitation(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    int_a  = ref(pm, nw, :xfmr, t, "INT_A")
    int_b  = ref(pm, nw, :xfmr, t, "INT_B")
    grad_a = ref(pm, nw, :xfmr, t, "GRAD_A")
    grad_b = ref(pm, nw, :xfmr, t, "GRAD_B")

    current_harmonics_ntws = [1]

    if nw in current_harmonics_ntws
        constraint_transformer_core_excitation(pm, nw, t, int_a, int_b, grad_a, grad_b)
    else 
        constraint_transformer_core_excitation(pm, nw, t)
    end
end

""
function constraint_transformer_core_voltage_drop(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    f_idx = (t,f_bus,t_bus)
    
    xsc = ref(pm, nw, :xfmr, t, "xsc")
    
    constraint_transformer_core_voltage_drop(pm, nw, t, f_idx, xsc)
end

"" 
function constraint_transformer_core_voltage_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    t_idx = (t,t_bus,f_bus)
    
    tr = ref(pm, nw, :xfmr, t, "tr")
    ti = ref(pm, nw, :xfmr, t, "ti")

    constraint_transformer_core_voltage_balance(pm, nw, t, t_idx, tr, ti)
end

""
function constraint_transformer_core_current_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    f_idx = (t,f_bus,t_bus)
    t_idx = (t,t_bus,f_bus)
    
    tr = ref(pm, nw, :xfmr, t, "tr")
    ti = ref(pm, nw, :xfmr, t, "ti")

    constraint_transformer_core_current_balance(pm, nw, t, f_idx, t_idx, tr, ti)
end

""
function constraint_transformer_winding_config(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    nh = pm.data["harmonics"]["$nw"]
    
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    bus = [f_bus, t_bus]
    idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]

    r  = [ref(pm, nw, :xfmr, t, nk) for nk in ["r1","r2"]]
    re = [ref(pm, nw, :xfmr, t, nk) for nk in ["re1","re2"]]
    xe = [ref(pm, nw, :xfmr, t, nk) for nk in ["xe1","xe2"]]

    gnd = [ref(pm, nw, :xfmr, t, nk) for nk in ["gnd1","gnd2"]]

    for w in 1:2
        constraint_transformer_winding_config(pm, nw, nh, bus[w], idx[w], r[w], re[w], xe[w], gnd[w])
    end
end

""
function constraint_transformer_winding_current_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    nh = pm.data["harmonics"]["$nw"]
    
    f_bus = ref(pm, nw, :xfmr, t, "f_bus")
    t_bus = ref(pm, nw, :xfmr, t, "t_bus")
    idx = [(t,f_bus,t_bus), (t,t_bus,f_bus)]
    
    r = [ref(pm, nw, :xfmr, t, nk) for nk in ["r1","r2"]]
    # b_sh = [ref(pm, nw, :xfmr, t, nk) for nk in ["b_sh1","b_sh2"]] ## TODO
    # g_sh = [ref(pm, nw, :xfmr, t, nk) for nk in ["g_sh1","g_sh2"]]
    b_sh, g_sh = [0.0,0.0], [0.0,0.0]

    cnf = [ref(pm, nw, :xfmr, t, nk) for nk in ["cnf1","cnf2"]]

    for w in 1:2
        constraint_transformer_winding_current_balance(pm, nw, nh, idx[w], r[w], b_sh[w], g_sh[w], cnf[w])
    end
end

""
function constraint_voltage_magnitude_rms(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = ref(pm, nw, :bus, i, "vmin")
    vmax = ref(pm, nw, :bus, i, "vmax")

    constraint_voltage_magnitude_rms(pm, i, vmin, vmax)
end