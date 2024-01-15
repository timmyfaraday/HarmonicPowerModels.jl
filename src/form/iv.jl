################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

## variables
# bus
""
function variable_bus_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# branch 
""
function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# xfmr
""
function variable_transformer_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_transformer_voltage_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_transformer_current(pm::AbstractIVRModel; nw::Int=nw_id_default(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_magnetizing_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_magnetizing_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# load 
""
function variable_gen_current(pm::AbstractIVRModel; nw::Int=nw_id_default(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# load 
""
function variable_load_current(pm::AbstractIVRModel; nw::Int=nw_id_default(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_load_current_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

## objective
""
function objective_power_flow(pm::_PMs.AbstractIVRModel)
    JuMP.@objective(pm.model, Min, 0.0)
end
""
function objective_maximum_hosting_capacity(pm::_PMs.AbstractIVRModel)
    cmd = [var(pm, n, :cmd, l)  for n in _PMs.nw_ids(pm) 
                                for l in _PMs.ids(pm, :load, nw=n) 
                                if n ≠ 1]
    
    JuMP.@objective(pm.model, Max, sum(cmd))
end
""
function objective_voltage_distortion_minimization(pm::_PMs.AbstractIVRModel; bus_id=1) # TODO @F: this needs to be generalized, suggestion find the bus where active filter is connected
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
function constraint_ref_bus(pm::AbstractIVRModel, n::Int, i::Int, vref)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vr == vref)
    JuMP.@constraint(pm.model, vi == 0.0)
end

# bus
""
function constraint_voltage_rms_limit(pm::AbstractIVRModel, i, vminrms, vmaxrms)
    w = [var(pm, n, :w, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, vminrms^2 <= sum(w)               )
    JuMP.@constraint(pm.model,              sum(w)  <= vmaxrms^2 )
end
""
function constraint_voltage_rms_limit(pm::dHHC_NLP, i, vminrms, vmaxrms)
    vr = [var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, vminrms^2 <= sum(vr^2 + vi^2)               )
    JuMP.@constraint(pm.model,              sum(vr^2 + vi^2)  <= vmaxrms^2 )
end
""
function constraint_voltage_rms_limit(pm::SOC_DHHC, i, vmaxrms)
    vr = [var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, [sqrt(vmaxrms^2 - 1.0^2); vcat(vr, vi)] in JuMP.SecondOrderCone()) # TODO: change 1.0 to input vmagfund
end
""
function constraint_voltage_thd_limit(pm::AbstractIVRModel, i, thdmax)
    w = [var(pm, n, :w, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(w[2:end]) <= thdmax^2 * w[1])
end
""
function constraint_voltage_thd_limit(pm::dHHC_NLP, i, thdmax)
    vr = [var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(vr[2:end]^2 + vi[2:end]^2) <= thdmax^2 * (vr[1]^2 + vi[1]^2))
end
""
function constraint_voltage_thd_limit(pm::SOC_DHHC, i, thdmax)
    vr = [var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, [thdmax * 1.0; vcat(vr, vi)] in JuMP.SecondOrderCone()) # TODO: change 1.0 to input vmagfund
end
""
function constraint_voltage_ihd_limit(pm::AbstractIVRModel, n::Int, i, ihdmax)
    v  = var(pm, 1, :w, i)
    w  = var(pm, n, :w, i)

    JuMP.@constraint(pm.model, w <= ihdmax^2 * v)
end
""
function constraint_voltage_ihd_limit(pm::dHHC_NLP, n::Int, i, ihdmax)
    vr = [var(pm, 1, :vr, i), var(pm, n, :vr, i)] 
    vi = [var(pm, 1, :vi, i), var(pm, n, :vr, i)]

    JuMP.@constraint(pm.model, vr[2]^2 + vi[2]^2 <= ihdmax^2 * (vr[1]^2 + vi[1]^2))
end
""
function constraint_voltage_ihd_limit(pm::SOC_DHHC, n::Int, i, ihdmax)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, [ihdmax * 1.0; vcat(vr, vi)] in JuMP.SecondOrderCone()) # TODO: change 1.0 to input vmagfund
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
""
function constraint_load_current_fixed_angle(pm::AbstractIVRModel, n::Int, l)
    crd = var(pm, n, :crd, l)
    cid = var(pm, n, :cid, l)
    cmd = var(pm, n, :cmd, l)

    JuMP.@constraint(pm.model, 0.0 <= crd)
    JuMP.@constraint(pm.model, 0.0 <= cid)

    JuMP.@constraint(pm.model, crd <= 1.0)                                      # @H: this needs to be deprecated, c_rating in the variable definition should set better bounds                          
    JuMP.@constraint(pm.model, cid <= 1.0)                                      # @H: this needs to be deprecated, c_rating in the variable definition should set better bounds

    JuMP.@constraint(pm.model, cid == crd)

    JuMP.@constraint(pm.model, cmd^2 <= crd^2 + cid^2)
end
""
function constraint_load_current_variable_angle(pm::AbstractIVRModel, n::Int, l, angmin, angmax)
    crd = var(pm, n, :crd, l)
    cid = var(pm, n, :cid, l)
    cmd = var(pm, n, :cmd, l)

    JuMP.@constraint(pm.model, cmd * min(sin(angmin), sin(angmax)) <= cid)
    JuMP.@constraint(pm.model, cid <= cmd * max(sin(angmin), sin(angmax)))

    JuMP.@constraint(pm.model, cmd * min(cos(angmin), cos(angmax)) <= crd)
    JuMP.@constraint(pm.model, crd <= cmd * max(cos(angmin), cos(angmax)))

    JuMP.@constraint(pm.model, cmd^2 <= crd^2 + cid^2)
end
""
function constraint_load_current_variable_angle(pm::SOC_DHHC, n::Int, l, angmin, angmax)
    crd = var(pm, n, :crd, l)
    cid = var(pm, n, :cid, l)
    cmd = var(pm, n, :cmd, l)

    JuMP.@constraint(pm.model, cmd * min(sin(angmin), sin(angmax)) <= cid)
    JuMP.@constraint(pm.model, cid <= cmd * max(sin(angmin), sin(angmax)))

    JuMP.@constraint(pm.model, cmd * min(cos(angmin), cos(angmax)) <= crd)
    JuMP.@constraint(pm.model, crd <= cmd * max(cos(angmin), cos(angmax)))

    # NB: This changes the direction of the inequality.
    JuMP.@constraint(pm.model, [1/sqrt(2) * cmd; 1/sqrt(2) * cmd; vcat(crd, cid)] in JuMP.RotatedSecondOrderCone())
end
""
function constraint_load_current_variable_angle_relative(pm::AbstractIVRModel, n::Int, l, bus_idx, c1, c2)
    crd = var(pm, n, :crd, l)
    cid = var(pm, n, :cid, l)
    cmd = var(pm, n, :cmd, l)
    vr = var(pm, n, :vr, bus_idx)
    vi = var(pm, n, :vi, bus_idx)

    JuMP.@NLconstraint(pm.model, (vr^2 + vi^2) * crd^2 == cmd^2 * (vr * c1 - vi * c2)^2)
    JuMP.@NLconstraint(pm.model, (vr^2 + vi^2) * cid^2 == cmd^2 * (vi * c1 + vr * c2)^2)
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

    JuMP.@constraint(pm.model,  csrt_fr 
                                + tr * csrt_to 
                                + ti * csit_to 
                                    == 
                                cert 
                                + ert / rsh
                    )
    JuMP.@constraint(pm.model,  csit_fr 
                                + tr * csit_to 
                                - ti * csrt_to
                                    == 
                                ceit 
                                + eit / rsh
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
    if is_zero_sequence(n) && cnf in ['D'] && r ≠ 0.0
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit - vrt / r)
        JuMP.@constraint(pm.model, cit == csit - g_sh * vit - b_sh * vrt - vit / r)
    end
end
""
function constraint_transformer_winding_current_rms_limit(pm::AbstractIVRModel, idx, c_rating)
    crt =  [var(pm, n, :crt, idx) for n in sorted_nw_ids(pm)]
    cit =  [var(pm, n, :cit, idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crt.^2 + cit.^2) <= c_rating^2)
end