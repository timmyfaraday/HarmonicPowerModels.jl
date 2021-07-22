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

    ## bounds are needed

    report && _PMs.sol_component_value(pm, nw, :xfmr, :ert, _PMs.ids(pm, nw, :xfmr), ert)
end

""
function variable_transformer_voltage_excitation_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    eit = _PMs.var(pm, nw)[:eit] = JuMP.@variable(pm.model,
            [t in _PMs.ids(pm, nw, :xfmr)], base_name="$(nw)_eit",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :xfmr, t), "eit_start", 0.0)
    )

    ## bounds are needed

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