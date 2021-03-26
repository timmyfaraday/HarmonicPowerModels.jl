
# variables
""
function variable_transformer_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_delta_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_delta_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_excitation_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_excitation_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# constraints
""
function constraint_transformer_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function constraint_current_winding(pm::AbstractIVRModel, n::Int, i, t, w, g_sh, b_sh, earthed)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
    
    ctr = var(pm, n, :ctr, t)[w]
    cti = var(pm, n, :cti, t)[w]

    ctsr = var(pm, n, :ctsr, t)[w]
    ctsi = var(pm, n, :ctsi, t)[w]

    if is_zero_sequence(n) && !earthed
        JuMP.@constraint(pm.model, ctsr == 0.0)
        JuMP.@constraint(pm.model, ctsi == 0.0)
    end

    JuMP.@constraint(pm.model, ctr == ctsr + g_sh * vr - b_sh * vi)
    JuMP.@constraint(pm.model, cti == ctsi + g_sh * vi + b_sh * vr)
end

""
function constraint_voltage_drop_winding(pm::AbstractIVRModel, n::Int, i, t, w, r, x)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    er = var(pm, n, :er, t)[w]
    ei = var(pm, n, :ei, t)[w]

    ctsr = var(pm, n, :ctsr, t)[w]
    ctsi = var(pm, n, :ctsi, t)[w]

    JuMP.@constraint(pm.model, er == vr - r * ctsr + x * csti)
    JuMP.@constraint(pm.model, ei == vi - r * ctsi - x * cstr)
end

""
function constraint_current_balance_transformer(pm::AbstractIVRModel, n::Int, t, tr, ti, config, windings)
    ctdr = var(pm, n, :ctdr, t)
    ctdi = var(pm, n, :ctdi, t)
    
    cter = var(pm, n, :cter, t)
    ctei = var(pm, n, :ctei, t)

    ctsr = var(pm, n, :ctsr, t)
    ctsi = var(pm, n, :ctsi, t)

    if !is_zero_sequence(n) || !any(config, :Delta)
        JuMP.@constraint(pm.model, ctdr == 0.0)
        JuMP.@constraint(pm.model, ctdi == 0.0)
    end

    JuMP.@constraint(pm.model, sum(tr[w] * ctsr[w] - ti[w] * ctsi[w] for w in windings)
                                == cter + ctdr
                    )
    JuMP.@constraint(pm.model, sum(tr[w] * ctsi[w] - ti[w] * ctsr[w] for w in windings)
                                == ctei + ctdi
                    )
end

""
function constraint_voltage_transformer(pm::AbstractIVRModel, t)
    er = [var(pm, n, :er, t) for n in nws(pm)]
    ei = [var(pm, n, :ei, t) for n in nws(pm)]

    cter = [var(pm, n, :cter, t) for n in nws(pm)]
    ctei = [var(pm, n, :ctei, t) for n in nws(pm)]

    ### Fre's spline magix
end