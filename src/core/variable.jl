""
function variable_transformer_voltage_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vrt = _PMs.var(pm, nw)[:vrt] = JuMP.@variable(pm.model, 
        [(t,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_vrt",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "vrt_start", 1.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :vrt_fr, :vrt_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), vrt)
end

""
function variable_transformer_voltage_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vit = _PMs.var(pm, nw)[:vit] = JuMP.@variable(pm.model,
        [(t,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_vit",
        start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "vit_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :vit_fr, :vit_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), vit)
end

""
function variable_transformer_voltage_excitation_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ert = _PMs.var(pm, nw)[:ert] = JuMP.@variable(pm.model,
            [t in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_ert",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "ert_start", 0.0)
    )

    if bounded
        for (t, xfmr) in ref(pm, nw, :xfmr)
            JuMP.set_lower_bound(ert[t], xfmr["ert_min"])
            JuMP.set_upper_bound(ert[t], xfmr["ert_max"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :ert, _PMs.ids(pm, nw, :xfmr), ert)
end

""
function variable_transformer_voltage_excitation_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    eit = _PMs.var(pm, nw)[:eit] = JuMP.@variable(pm.model,
            [t in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_eit",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "eit_start", 0.0)
    )

    if bounded
        for (t, xfmr) in ref(pm, nw, :xfmr)
            JuMP.set_lower_bound(eit[t], xfmr["eit_min"])
            JuMP.set_upper_bound(eit[t], xfmr["eit_max"])
        end
    end

    report && _PMs.sol_component_value(pm, nw, :xfmr, :eit, _PMs.ids(pm, nw, :xfmr), eit)
end

""
function variable_transformer_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crt = _PMs.var(pm, nw)[:crt] = JuMP.@variable(pm.model,
            [(t,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_crt",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "crt_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :crt_fr, :crt_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), crt)
end

""
function variable_transformer_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cit = _PMs.var(pm, nw)[:cit] = JuMP.@variable(pm.model,
            [(t,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_cit",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "cit_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :cit_fr, :cit_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), cit)
end

""
function variable_transformer_current_series_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    csrt = _PMs.var(pm, nw)[:csrt] = JuMP.@variable(pm.model,
            [(t,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_csrt",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "csrt_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :csrt_fr, :csrt_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), csrt)
end

""
function variable_transformer_current_series_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    csit = _PMs.var(pm, nw)[:csit] = JuMP.@variable(pm.model,
            [(t,i,j) in _PMs.ref(pm, nw, :xfmr_arcs)], base_name="$(nw)_csit",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "csit_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :xfmr, :csit_fr, :csit_to, _PMs.ref(pm, nw, :xfmr_arcs_from), _PMs.ref(pm, nw, :xfmr_arcs_to), csit)
end

""
function variable_transformer_current_excitation_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cert = _PMs.var(pm, nw)[:cert] = JuMP.@variable(pm.model,
            [t in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_cert",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "cert_start", 0.0)
    )

    ## bounds are needed

    report && _PMs.sol_component_value(pm, nw, :xfmr, :cert, _PMs.ids(pm, nw, :xfmr), cert)
end

""
function variable_transformer_current_excitation_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ceit = _PMs.var(pm, nw)[:ceit] = JuMP.@variable(pm.model,
            [t in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_ceit",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "ceit_start", 0.0)
    )

    ## bounds are needed

    report && _PMs.sol_component_value(pm, nw, :xfmr, :ceit, _PMs.ids(pm, nw, :xfmr), ceit)
end

""
function variable_load_current(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd = _PMs.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
            [d in _PMs.ids(pm, nw, :load)], base_name="$(nw)_crd",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, d), "crd_start", 0.0)
    )

    report && _PMs.sol_component_value(pm, nw, :load, :crd, _PMs.ids(pm, nw, :load), crd)

    cid = _PMs.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
            [d in _PMs.ids(pm, nw, :load)], base_name="$(nw)_cid",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :load, d), "cid_start", 0.0)
    )

    report && _PMs.sol_component_value(pm, nw, :load, :cid, _PMs.ids(pm, nw, :load), cid)

    if bounded
        for (d, load) in ref(pm, nw, :load)
            c_rating = load["c_rating"]
            JuMP.set_lower_bound(crd[d], -c_rating)
            JuMP.set_upper_bound(crd[d],  c_rating)
            JuMP.set_lower_bound(cid[d], -c_rating)
            JuMP.set_upper_bound(cid[d],  c_rating)
        end
    end

    # store active and reactive power expressions for use in objective + post processing
    pd = Dict()
    qd = Dict()
    for (i,load) in ref(pm, nw, :load)
        busid = load["load_bus"]
        vr = var(pm, nw, :vr, busid)
        vi = var(pm, nw, :vi, busid)
        crd = var(pm, nw, :crd, i)
        cid = var(pm, nw, :cid, i)
        pd[i] = JuMP.@NLexpression(pm.model, vr*crd  + vi*cid)
        qd[i] = JuMP.@NLexpression(pm.model, vi*crd  - vr*cid)
    end
    var(pm, nw)[:pd] = pd
    var(pm, nw)[:qd] = qd
    report && _PMs.sol_component_value(pm, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && _PMs.sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
end

""
function variable_gen_current(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crg = _PMs.var(pm, nw)[:crg] = JuMP.@variable(pm.model,
            [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_crg",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :gen, g), "crg_start", 0.0)
    )

    report && _PMs.sol_component_value(pm, nw, :gen, :crg, _PMs.ids(pm, nw, :gen), crg)

    cig = _PMs.var(pm, nw)[:cig] = JuMP.@variable(pm.model,
            [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_cig",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :gen, g), "cig_start", 0.0)
    )

    report && _PMs.sol_component_value(pm, nw, :gen, :cig, _PMs.ids(pm, nw, :gen), cig)

    if bounded
        for (g, gen) in ref(pm, nw, :gen)
            c_rating = gen["c_rating"]
            JuMP.set_lower_bound(crg[g], -c_rating)
            JuMP.set_upper_bound(crg[g],  c_rating)
            JuMP.set_lower_bound(cig[g], -c_rating)
            JuMP.set_upper_bound(cig[g],  c_rating)
        end
    end

    # store active and reactive power expressions for use in objective + post processing
    pg = Dict()
    qg = Dict()
    for (i,gen) in ref(pm, nw, :gen)
        busid = gen["gen_bus"]
        vr = var(pm, nw, :vr, busid)
        vi = var(pm, nw, :vi, busid)
        crg = var(pm, nw, :crg, i)
        cig = var(pm, nw, :cig, i)
        pg[i] = JuMP.@NLexpression(pm.model, vr*crg  + vi*cig)
        qg[i] = JuMP.@NLexpression(pm.model, vi*crg  - vr*cig)
    end
    var(pm, nw)[:pg] = pg
    var(pm, nw)[:qg] = qg
    report && _PMs.sol_component_value(pm, nw, :load, :pg, ids(pm, nw, :gen), pg)
    report && _PMs.sol_component_value(pm, nw, :load, :qg, ids(pm, nw, :gen), qg)
end

