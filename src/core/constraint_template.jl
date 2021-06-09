
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

    constraint_current_balance(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_xfmr, bus_gens, bus_loads, bus_gs, bus_bs)
end

""
function constraint_transformer_core_voltage_drop(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    xs = ref(pm, nw, :xfmr, t, "xs")

    constraint_transformer_core_voltage_drop(pm, nw, t, xs)
end

"" 
function constraint_transformer_core_voltage_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    tr = ref(pm, nw, :xfmr, t, "tr")
    ti = ref(pm, nw, :xfmr, t, "ti")

    constraint_transformer_core_voltage_balance(pm, nw, t, tr, ti)
end

""
function constraint_transformer_core_current_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    tr = ref(pm, nw, :xfmr, t, "tr")
    ti = ref(pm, nw, :xfmr, t, "ti")

    windings = ref(pm, nw, :xfmr, t, "windings")

    constraint_current_balance_transformer_core(pm, nw, t, tr, ti, windings)
end

""
function constraint_transformer_winding_config(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    r = ref(pm, nw, :xfmr, t, "r")

    re = ref(pm, nw, :xfmr, t, "re")
    xe = ref(pm, nw, :xfmr, t, "xe")

    buses = ref(pm, nw, :xfmr, t, "buses")
    earthed = ref(pm, nw, :xfmr, t, "earthed")

    for w in ref(pm, nw, :xfmr, t, "windings")
        constraint_transformer_winding_config(pm, nw, buses[w], t, w, r[w], re[w], xe[w], earthed[w])
    end
end

""
function constraint_transformer_winding_current_balance(pm::AbstractPowerModel, t::Int; nw::Int=nw_id_default)
    b_sh = ref(pm, nw, :xfmr, t, "b_sh")
    g_sh = ref(pm, nw, :xfmr, t, "g_sh")

    config = ref(pm, nw, :xfmr, t, "config")

    for w in ref(pm, nw, :xfmr, t, "windings")
        constraint_transformer_winding_current_balance(pm, nw, t, w, b_sh[w], g_sh[w], config[w])
    end
end

""
function constraint_voltage_magnitude_rms(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = ref(pm, nw, :bus, i, "vmin")
    vmax = ref(pm, nw, :bus, i, "vmax")

    constraint_voltage_magnitude_rms(pm, nw, i, vmin, vmax)
end