
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(con(pm, nw), :kcl_cr)
        con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(con(pm, nw), :kcl_ci)
        con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_xfmr = ref(pm, nw, :bus_arcs_xfmr, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_arcs_xfmr, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs)
end

""
function constraint_current_transformer(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    r = ref(pm, nw, :xfmr, t, "r")
    x = ref(pm, nw, :xfmr, t, "x")
    b_sh = ref(pm, nw, :xfmr, t, "b_sh")
    g_sh = ref(pm, nw, :xfmr, t, "g_sh")
    
    tr = ref(pm, nw, :xfmr, t, "tr")
    ti = ref(pm, nw, :xfmr, t, "ti")

    buses = ref(pm, nw, :xfmr, t, "buses")
    config = ref(pm, nw, :xfmr, t, "config")
    earthed = ref(pm, nw, :xfmr, t, "earthed")
    windings = ref(pm, nw, :xfmr, t, "windings")

    for w in windings
        constraint_current_winding(pm, nw, buses[w], t, w, b_sh[w], g_sh[w], earthed[w])
        constraint_voltage_drop_winding(pm, nw, buses[w], t, w, r[w], x[w])
    end
    constraint_current_balance_transformer(pm, nw, t, tr, ti, config, windings)
end

""
function constraint_voltage_transformer(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)

    constraint_voltage_transformer(pm::AbstractIVRModel, t)
end

""
function constraint_voltage_magnitude_rms(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = ref(pm, nw, :bus, i, "vmin")
    vmax = ref(pm, nw, :bus, i, "vmax")

    constraint_voltage_magnitude_rms(pm, i, vmin, vmax)
end