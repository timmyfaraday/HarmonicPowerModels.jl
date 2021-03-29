
# variables
""
function variable_transformer_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_delta_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_delta_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_transformer_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# constraints
""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_xfmr, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    cr = var(pm, n, :cr)
    ci = var(pm, n, :ci)
    crt = var(pm, n, :crt)
    cit = var(pm, n, :cit)
    crdc = var(pm, n, :crdc)
    cidc = var(pm, n, :cidc)

    crd = var(pm, n, :crd)
    cid = var(pm, n, :cid)
    crg = var(pm, n, :crg)
    cig = var(pm, n, :cig)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crt[t][w] for (t,w) in bus_arcs_xfmr)
                                + sum(crdc[d] for d in bus_arcs_dc)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cit[t][w] for (t,w) in bus_arcs_xfmr)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

""
function constraint_current_winding(pm::AbstractIVRModel, n::Int, i, t, w, b_sh, g_sh, earthed)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
    
    crt = var(pm, n, :crt, t)[w]
    cit = var(pm, n, :cit, t)[w]

    csrt = var(pm, n, :csrt, t)[w]
    csit = var(pm, n, :csit, t)[w]

    if is_zero_sequence(n) && !earthed
        JuMP.@constraint(pm.model, csrt == 0.0)
        JuMP.@constraint(pm.model, csit == 0.0)
    end

    JuMP.@constraint(pm.model, crt == csrt + g_sh * vr - b_sh * vi)
    JuMP.@constraint(pm.model, cit == csit + g_sh * vi + b_sh * vr)
end

""
function constraint_voltage_drop_winding(pm::AbstractIVRModel, n::Int, i, t, w, r, x)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    er = var(pm, n, :er, t)[w]
    ei = var(pm, n, :ei, t)[w]

    csrt = var(pm, n, :csrt, t)[w]
    csit = var(pm, n, :csit, t)[w]

    JuMP.@constraint(pm.model, er == vr - r * csrt + x * csti)
    JuMP.@constraint(pm.model, ei == vi - r * csit - x * cstr)
end

""
function constraint_current_balance_transformer(pm::AbstractIVRModel, n::Int, t, tr, ti, config, windings)
    cdrt = var(pm, n, :cdrt, t)
    cdit = var(pm, n, :cdit, t)
    
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    csrt = var(pm, n, :csrt, t)
    csit = var(pm, n, :csit, t)

    if !is_zero_sequence(n) || !any(config, :Delta)
        JuMP.@constraint(pm.model, cdrt == 0.0)
        JuMP.@constraint(pm.model, cdit == 0.0)
    end

    JuMP.@constraint(pm.model, sum(tr[w] * csrt[w] - ti[w] * csit[w] for w in windings)
                                == cdrt + cert
                    )
    JuMP.@constraint(pm.model, sum(tr[w] * csit[w] - ti[w] * csrt[w] for w in windings)
                                == cdit + ceit
                    )
end

""
function constraint_voltage_transformer(pm::AbstractIVRModel, t)
    er = [var(pm, n, :er, t) for n in nws(pm)]
    ei = [var(pm, n, :ei, t) for n in nws(pm)]

    cert = [var(pm, n, :cert, t) for n in nws(pm)]
    ceit = [var(pm, n, :ceit, t) for n in nws(pm)]
    
    ### Fre's spline magix
end

""
function constraint_voltage_magnitude_rms(pm::AbstractIVRModel, i::Int, vmin, vmax)
    vr = [var(pm, nw, :vr, i) for nw in nw_ids(pm)]
    vi = [var(pm, nw, :vi, i) for nw in nw_ids(pm)]

    JuMP.@constraint(pm.model, vmin^2 <= sum(vr.^2 + vi.^2))
    JuMP.@constraint(pm.model, sum(vr.^2 + vi.^2) <= vmax^2)
end