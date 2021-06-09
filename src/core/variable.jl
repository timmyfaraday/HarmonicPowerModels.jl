
""
function variable_transformer_voltage_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    vrt = var(pm, nw)[:vrt] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_vrt_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "vrt_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :vrt, ids(pm, nw, :xfmr), vrt)
end

""
function variable_transformer_voltage_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    vit = var(pm, nw)[:vit] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_vit_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "vit_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :vit, ids(pm, nw, :xfmr), vit)
end

""
function variable_transformer_voltage_excitation_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ert = var(pm, nw)[:ert] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :xfmr)], base_name="$(nw)_ert",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "ert_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :ert, ids(pm, nw, :xfmr), ert)
end

""
function variable_transformer_voltage_excitation_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    eit = var(pm, nw)[:eit] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :xfmr)], base_name="$(nw)_eit",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "eit_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :eit, ids(pm, nw, :xfmr), eit)
end

""
function variable_transformer_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    crt = var(pm, nw)[:crt] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_crt_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "crt_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :crt, ids(pm, nw, :xfmr), crt)
end

""
function variable_transformer_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    cit = var(pm, nw)[:cit] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_cit_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "cit_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :cit, ids(pm, nw, :xfmr), cit)
end

""
function variable_transformer_current_series_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    csrt = var(pm, nw)[:csrt] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_csrt_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "csrt_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :csrt, ids(pm, nw, :xfmr), csrt)
end

""
function variable_transformer_current_series_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    csit = var(pm, nw)[:csit] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_csit_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "csit_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :csit, ids(pm, nw, :xfmr), csit)
end

""
function variable_transformer_current_excitation_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cert = var(pm, nw)[:cert] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :xfmr)], base_name="$(nw)_cert",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "cert_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :cert, ids(pm, nw, :xfmr), cert)
end

""
function variable_transformer_current_excitation_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ceit = var(pm, nw)[:ceit] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :xfmr)], base_name="$(nw)_ceit",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "ceit_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :ceit, ids(pm, nw, :xfmr), ceit)
end