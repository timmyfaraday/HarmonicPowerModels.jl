################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

## variables
# xfmr
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

    expression_transformer_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_transformer_excitation_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

## objective
""
function objective_power_flow(pm::AbstractIVRModel)
    JuMP.@objective(pm.model, Min, 0.0)
end
""
function objective_voltage_distortion_minimization(pm::AbstractIVRModel; bus_id=6) # @F: this needs to be generalized, suggestion find the bus where active filter is connected
    vr = [var(pm, n, :vr, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]
    vi = [var(pm, n, :vi, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]

    JuMP.@objective(pm.model, Min, sum(vr.^2 + vi.^2))
end

## constraints
# active filter
""
function constraint_active_filter(pm::AbstractIVRModel, n::Int, g)
    pg = [var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

    JuMP.@NLconstraint(pm.model, pg[1] == 0)
    JuMP.@NLconstraint(pm.model, sum(pg[n] for n in 1:length(pg)) == 0)         # NB: This ugly sum construction is needed to allow sum of expressions pg
end

# ref bus
""
function constraint_ref_bus(pm::AbstractIVRModel, n::Int, i::Int)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    if n == 1 
        # for fundemental frequency: fix reference angle
        JuMP.@constraint(pm.model, vr == 1.0)
        JuMP.@constraint(pm.model, vi == 0.0)
    else 
        # for non-fundamental frequency: fix harmonic voltage to 0+j0
        JuMP.@constraint(pm.model, vr == 0.0)
        JuMP.@constraint(pm.model, vi == 0.0)
    end
end

# bus
""
function constraint_voltage_rms_limit(pm::AbstractIVRModel, i, vminrms, vmaxrms)
    w = [var(pm, n, :w, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, vminrms^2 <= sum(w)               )
    JuMP.@constraint(pm.model,              sum(w)  <= vmaxrms^2 )
end
""
function constraint_voltage_thd_limit(pm::AbstractIVRModel, i, thdmax)
    w = [var(pm, n, :w, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(w[2:end]) <= thdmax^2 * w[1])
end
""
function constraint_voltage_ihd_limit(pm::AbstractIVRModel, n::Int, i, ihd)
    v  = var(pm, 1, :w, i)
    w  = var(pm, n, :w, i)

    JuMP.@constraint(pm.model, w <= ihd * v)
end
""
function constraint_voltage_magnitude_sqr(pm::AbstractIVRModel, n::Int, i)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
    w  = var(pm, n, :w, i)
    
    JuMP.@constraint(pm.model, w == vr^2  + vi^2)                               # @F: Changed this to an equality, as the inequality only works when minimizing the distortion
end
""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_xfmr, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    cr = var(pm, n, :cr)
    ci = var(pm, n, :ci)
    crt = var(pm, n, :crt)
    cit = var(pm, n, :cit)

    crg = var(pm, n, :crg)
    cig = var(pm, n, :cig)
    crd = var(pm, n, :crd)
    cid = var(pm, n, :cid)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crt[t] for t in bus_arcs_xfmr)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cit[t] for t in bus_arcs_xfmr)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )

end

# branch
function constraint_current_rms_limit(pm::AbstractIVRModel, f_idx, t_idx, c_rating)
    crf =  [var(pm, n, :cr, f_idx) for n in sorted_nw_ids(pm)]
    cif =  [var(pm, n, :ci, f_idx) for n in sorted_nw_ids(pm)]

    crt =  [var(pm, n, :cr, t_idx) for n in sorted_nw_ids(pm)]
    cit =  [var(pm, n, :ci, t_idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crf^2 + cif^2) <= c_rating^2)
    JuMP.@constraint(pm.model, sum(crt^2 + cit^2) <= c_rating^2)
end

# load
""
function constraint_load_constant_power(pm::AbstractIVRModel, n::Int, l, i, pd, qd)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
    crd = var(pm, n, :crd, l)
    cid = var(pm, n, :cid, l)

    JuMP.@constraint(pm.model, pd == vr*crd  + vi*cid)
    JuMP.@constraint(pm.model, qd == vi*crd  - vr*cid)
end
""
function constraint_load_constant_current(pm::AbstractIVRModel, n::Int, l, mult)
    crd = var(pm, n, :crd, l)
    cid = var(pm, n, :cid, l)
    fund_crd = var(pm, 1, :crd, l)
    fund_cid = var(pm, 1, :cid, l)

    JuMP.@constraint(pm.model, crd == mult * fund_crd)
    JuMP.@constraint(pm.model, cid == mult * fund_cid)
end

# xfmr
""
function constraint_transformer_core_excitation(pm::AbstractIVRModel, n::Int, t, int_a, int_b)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    et = reduce(vcat,[[var(pm, nw, :ert, t),var(pm, nw, :eit, t)] 
                                for nw in _PMs.ref(pm, n, :xfmr, t, "Hᴱ")])

    sym_exc_a = Symbol("exc_a_",n,"_",t)
    sym_exc_b = Symbol("exc_b_",n,"_",t)

    JuMP.register(pm.model, sym_exc_a, length(et), int_a; autodiff=true)
    JuMP.register(pm.model, sym_exc_b, length(et), int_b; autodiff=true)

    JuMP.add_nonlinear_constraint(pm.model, :($(cert) == $(sym_exc_a)($(et...))))
    JuMP.add_nonlinear_constraint(pm.model, :($(ceit) == $(sym_exc_b)($(et...))))
end
""
function constraint_transformer_core_excitation(pm::AbstractIVRModel, n::Int, t)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    JuMP.@constraint(pm.model, cert == 0.0)
    JuMP.@constraint(pm.model, ceit == 0.0)
end
""
function constraint_transformer_core_voltage_drop(pm::_PMs.AbstractIVRModel, n::Int, t, f_idx, xsc)
    ert = var(pm, n, :ert, t)
    eit = var(pm, n, :eit, t)

    vrt = var(pm, n, :vrt, f_idx)
    vit = var(pm, n, :vit, f_idx)

    csrt = var(pm, n, :csrt, f_idx)
    csit = var(pm, n, :csit, f_idx)
    
    JuMP.@constraint(pm.model, vrt == ert - xsc * csit)
    JuMP.@constraint(pm.model, vit == eit + xsc * csrt)
end
""
function constraint_transformer_core_voltage_balance(pm::_PMs.AbstractIVRModel, n::Int, t, t_idx, tr, ti)
    ert = var(pm, n, :ert, t)
    eit = var(pm, n, :eit, t)

    vrt = var(pm, n, :vrt, t_idx)
    vit = var(pm, n, :vit, t_idx)

    JuMP.@constraint(pm.model, vrt == tr * ert - ti * eit)
    JuMP.@constraint(pm.model, vit == tr * eit + ti * ert)
end
""
function constraint_transformer_core_current_balance(pm::_PMs.AbstractIVRModel, n::Int, t, f_idx, t_idx, tr, ti, rsh)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    csrt_fr = var(pm, n, :csrt, f_idx)
    csit_fr = var(pm, n, :csit, f_idx)

    csrt_to = var(pm, n, :csrt, t_idx)
    csit_to = var(pm, n, :csit, t_idx)

    ert = var(pm, n, :ert, t)
    eit = var(pm, n, :eit, t)

    JuMP.@constraint(pm.model, csrt_fr + tr * csrt_to + ti * csit_to 
                                == cert + ert/rsh
                    )
    JuMP.@constraint(pm.model, csit_fr + tr * csit_to - ti * csrt_to
                                == ceit + eit/rsh
                    )
end
""
function constraint_transformer_winding_config(pm::AbstractIVRModel, n::Int, i, idx, r, re, xe, gnd)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    vrt = var(pm, n, :vrt, idx)
    vit = var(pm, n, :vit, idx)

    crt = var(pm, n, :crt, idx)
    cit = var(pm, n, :cit, idx)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(n)
        JuMP.@constraint(pm.model, vrt == vr - r * crt)
        JuMP.@constraint(pm.model, vit == vi - r * cit)
    end

    # h ∈ 𝓗⁰, gnd == true -> cnf ∈ {Ye, Ze}
    if is_zero_sequence(n) && gnd == 1
        JuMP.@constraint(pm.model, vrt == vr - (r + 3re) * crt + 3xe * cit)
        JuMP.@constraint(pm.model, vit == vi - (r + 3re) * cit - 3xe * crt)
    end

    # h ∈ 𝓗⁰, gnd == false -> cnf ∈ {D, Y, Z}
    if is_zero_sequence(n) && gnd != 1
        JuMP.@constraint(pm.model, crt == 0)
        JuMP.@constraint(pm.model, cit == 0)
    end
end
""
function constraint_transformer_winding_current_balance(pm::AbstractIVRModel, n::Int, idx, r, b_sh, g_sh, cnf)
    vrt = var(pm, n, :vrt, idx)
    vit = var(pm, n, :vit, idx)
    
    crt = var(pm, n, :crt, idx)
    cit = var(pm, n, :cit, idx)

    csrt = var(pm, n, :csrt, idx)
    csit = var(pm, n, :csit, idx)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(n)
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit)
        JuMP.@constraint(pm.model, cit == csit - g_sh * vit - b_sh * vrt)
    end

    # h ∈ 𝓗⁰, cnf ∈ {Y(e), Z(e)}
    if is_zero_sequence(n) && cnf in ['Y','Z']
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit)
        JuMP.@constraint(pm.model, cit == csit - g_sh * vit - b_sh * vrt)
    end

    # h ∈ 𝓗⁰, cnf ∈ {D}
    if is_zero_sequence(n) && cnf in ['D']
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit - vrt / r)
        JuMP.@constraint(pm.model, crt == csit - g_sh * vit - b_sh * vrt - vit / r)
    end
end