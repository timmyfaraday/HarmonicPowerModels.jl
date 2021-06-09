
# variables
""
function variable_transformer_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_transformer_voltage_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_transformer_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# constraints
""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_arcs_xfmr, bus_gens, bus_gs, bus_bs)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    cr = var(pm, n, :cr)
    ci = var(pm, n, :ci)
    crt = var(pm, n, :crt)
    cit = var(pm, n, :cit)
    crdc = var(pm, n, :crdc)
    cidc = var(pm, n, :cidc)

    crg = var(pm, n, :crg)
    cig = var(pm, n, :cig)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crdc[d] for d in bus_arcs_dc)
                                + sum(crt[t] for t in bus_arcs_xfmr)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - (sum(pd for pd in values(bus_pd))*vr + sum(qd for qd in values(bus_qd))*vi)/(vr^2 + vi^2)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                + sum(cit[t] for t in bus_arcs_xfmr)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - (sum(pd for pd in values(bus_pd))*vi - sum(qd for qd in values(bus_qd))*vr)/(vr^2 + vi^2)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

""
function constraint_transformer_core_voltage_drop(pm::AbstractIVRModel, n::Int, t, xs)
    ert = var(pm, n, :ert, t)
    eit = var(pm, n, :eit, t)

    vrt = var(pm, n, :vrt, t)[1]
    vit = var(pm, n, :vit, t)[1]

    csrt = var(pm, n, :csrt, t)[1]
    csit = var(pm, n, :csit, t)[1]
    
    JuMP.@constraint(pm.model, vrt == ert - xs * csit)
    JuMP.@constraint(pm.model, vit == eit + xs * csrt)
end

""
function constraint_transformer_core_voltage_balance(pm::AbstractIVRModel, n::Int, t, tr, ti)
    ert = var(pm, n, :ert, t)
    eit = var(pm, n, :eit, t)

    vrt = var(pm, n, :vrt, t)[2]
    vit = var(pm, n, :vit, t)[2]

    JuMP.@constraint(pm.model, vrt == tr * ert + ti * eit)
    JuMP.@constraint(pm.model, vit == tr * eit + ti * ert)
end

""
function constraint_transformer_core_current_balance(pm::AbstractIVRModel, n::Int, t, tr, ti, windings)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    csrt = var(pm, n, :csrt, t)
    csit = var(pm, n, :csit, t)

    JuMP.@constraint(pm.model, sum(tr[w] * csrt[w] - ti[w] * csit[w] for w in windings)
                                == cert
                    )
    JuMP.@constraint(pm.model, sum(tr[w] * csit[w] - ti[w] * csrt[w] for w in windings)
                                == ceit
                    )
end

""
function constraint_transformer_winding_config(pm::AbstractIVRModel, n::Int, i, t, w, r, re, xe, earthed)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    vrt = var(pm, n, :vrt, t)[w]
    vit = var(pm, n, :vit, t)[w]

    crt = var(pm, n, :crt, t)[w]
    cit = var(pm, n, :cit, t)[w]

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(n)
        JuMP.@constraint(pm.model, vrt == vr - r * crt)
        JuMP.@constraint(pm.model, vit == vi - r * cit)
    end

    # h ∈ 𝓗⁰, conf ∈ {:Ye, :Ze}
    if is_zero_sequence(n) && earthed
        JuMP.@constraint(pm.model, vrt == vr - (r + 3re) * csrt + 3xe * csit)
        JuMP.@constraint(pm.model, vit == vi - (r + 3re) * csit - 3xe * cstr)
    end

    # h ∈ 𝓗⁰, conf ∈ {:D, :Y, :Z}
    if is_zero_sequence(n) && !earthed
        JuMP.@constraint(pm.model, crt == 0)
        JuMP.@constraint(pm.model, cit == 0)
    end
end

""
function constraint_transformer_winding_current_balance(pm::AbstractIVRModel, n::Int, t, w, r, b_sh, g_sh, config)
    vrt = var(pm, n, :vrt, t)[w]
    vit = var(pm, n, :vit, t)[w]
    
    crt = var(pm, n, :crt, t)[w]
    cit = var(pm, n, :cit, t)[w]

    csrt = var(pm, n, :csrt, t)[w]
    csit = var(pm, n, :csit, t)[w]

    if is_zero_sequence(n) && config == :D
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vr + b_sh * vi - vr / r)
        JuMP.@constraint(pm.model, crt == csit - g_sh * vi - b_sh * vr - vi / r)
    else 
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vr + b_sh * vi)
        JuMP.@constraint(pm.model, cit == csit - g_sh * vi - b_sh * vr)
    end
end

""
function constraint_voltage_magnitude_rms(pm::AbstractIVRModel, n::Int, i, vmin, vmax)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vmin^2 <= vr^2 + vi^2)
    JuMP.@constraint(pm.model, vr^2 + vi^2 <= vmax^2)
end