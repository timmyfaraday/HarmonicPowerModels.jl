
""
function variable_transformer_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => transformer["windings"] for (t,transformer) in ref(pm, nw, :transformer))

    ctr = var(pm, nw)[:ctr] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_ctr_$(t)",
            start = comp_start_value(ref(pm, nw, :transformer, t), "ctr_start", 0.0)
        ) for t in ids(pm, nw, :transformer)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :ctr, ids(pm, nw, :transformer), ctr)
end

""
function variable_transformer_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => transformer["windings"] for (t,transformer) in ref(pm, nw, :transformer))

    cti = var(pm, nw)[:cti] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_cti_$(t)",
            start = comp_start_value(ref(pm, nw, :transformer, t), "cti_start", 0.0)
        ) for t in ids(pm, nw, :transformer)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :cti, ids(pm, nw, :transformer), cti)
end

""
function variable_transformer_series_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => transformer["windings"] for (t,transformer) in ref(pm, nw, :transformer))

    ctsr = var(pm, nw)[:ctsr] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_ctsr_$(t)",
            start = comp_start_value(ref(pm, nw, :transformer, t), "ctsr_start", 0.0)
        ) for t in ids(pm, nw, :transformer)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :ctsr, ids(pm, nw, :transformer), ctsr)
end

""
function variable_transformer_series_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => transformer["windings"] for (t,transformer) in ref(pm, nw, :transformer))

    ctsi = var(pm, nw)[:ctsi] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_ctsi_$(t)",
            start = comp_start_value(ref(pm, nw, :transformer, t), "ctsi_start", 0.0)
        ) for t in ids(pm, nw, :transformer)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :ctsi, ids(pm, nw, :transformer), ctsi)
end

""
function variable_transformer_delta_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ctdr = var(pm, nw)[:ctdr] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :transformer)], base_name="$(nw)_ctdr",
            start = comp_start_value(ref(pm, nw, :transformer, t), "ctdr_start", 0.0)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :ctdr, ids(pm, nw, :transformer), ctdr)
end

""
function variable_transformer_delta_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ctdi = var(pm, nw)[:ctdi] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :transformer)], base_name="$(nw)_ctdi",
            start = comp_start_value(ref(pm, nw, :transformer, t), "ctdi_start", 0.0)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :ctdi, ids(pm, nw, :transformer), ctdi)
end

""
function variable_transformer_excitation_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cter = var(pm, nw)[:cter] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :transformer)], base_name="$(nw)_cter",
            start = comp_start_value(ref(pm, nw, :transformer, t), "cter_start", 0.0)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :cter, ids(pm, nw, :transformer), cter)
end

""
function variable_transformer_excitation_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ctei = var(pm, nw)[:ctei] = JuMP.@variable(pm.model,
            [t in ids(pm, nw, :transformer)], base_name="$(nw)_ctei",
            start = comp_start_value(ref(pm, nw, :transformer, t), "ctei_start", 0.0)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :ctei, ids(pm, nw, :transformer), ctei)
end

""
function variable_transformer_voltage_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => transformer["windings"] for (t,transformer) in ref(pm, nw, :transformer))

    etr = var(pm, nw)[:etr] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_etr_$(t)",
            start = comp_start_value(ref(pm, nw, :transformer, t), "etr_start", 0.0)
        ) for t in ids(pm, nw, :transformer)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :etr, ids(pm, nw, :transformer), etr)
end

""
function variable_transformer_voltage_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    windings = Dict(t => transformer["windings"] for (t,transformer) in ref(pm, nw, :transformer))

    eti = var(pm, nw)[:eti] = Dict(t => JuMP.@variable(pm.model,
            [w in windings[t]], base_name="$(nw)_eti_$(t)",
            start = comp_start_value(ref(pm, nw, :transformer, t), "eti_start", 0.0)
        ) for t in ids(pm, nw, :transformer)
    )

    ## bounds are needed

    report && _IM.sol_component_value(pm, nw, :transformer, :eti, ids(pm, nw, :transformer), eti)
end