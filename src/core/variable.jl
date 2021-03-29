
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
function variable_transformer_current_delta_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cdrt = var(pm, nw)[:cdrt] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :xfmr)], base_name="$(nw)_cdrt",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "cdrt_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :cdrt, ids(pm, nw, :xfmr), cdrt)
end

""
function variable_transformer_current_delta_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cdit = var(pm, nw)[:cdit] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :xfmr)], base_name="$(nw)_cdit",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "cdit_start", 0.0)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :cdit, ids(pm, nw, :xfmr), cdit)
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

""
function variable_transformer_voltage_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    ert = var(pm, nw)[:ert] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_ert_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "ert_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :ert, ids(pm, nw, :xfmr), ert)
end

""
function variable_transformer_voltage_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => xfmr["windings"] for (t,xfmr) in ref(pm, nw, :xfmr))

    eit = var(pm, nw)[:eit] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_eit_$(t)",
            start = _PMs.comp_start_value(ref(pm, nw, :xfmr, t), "eit_start", 0.0)
        ) for t in ids(pm, nw, :xfmr)
    )

    ## bounds are needed

    report && _IMs.sol_component_value(pm, nw, :xfmr, :eit, ids(pm, nw, :xfmr), eit)
end