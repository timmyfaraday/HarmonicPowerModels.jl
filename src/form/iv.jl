################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
################################################################################

## variables
# bus
""
function variable_bus_voltage(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# branch 
""
function variable_branch_current(pm::_PMs.AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# filter
""
function variable_filter_current(pm::_PMs.AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_filter_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_filter_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# xfmr
""
function variable_transformer_voltage(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_transformer_voltage_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_voltage_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_transformer_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_transformer_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_transformer_current_magnetizing_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_transformer_current_magnetizing_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# generator
""
function variable_gen_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# load 
""
function variable_load_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
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
function objective_voltage_distortion_minimization(pm::_PMs.AbstractIVRModel; bus_id=1) 
    vr = [_PMs.var(pm, n, :vr, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]
    vi = [_PMs.var(pm, n, :vi, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]

    JuMP.@objective(pm.model, Min, sum(vr.^2 + vi.^2))
end
""
function objective_maximum_hosting_capacity(pm::_PMs.AbstractIVRModel)
    cmd = [_PMs.var(pm, n, :cmd, l)  for n in _PMs.nw_ids(pm) 
                                for l in _PMs.ids(pm, :load, nw=n) 
                                if n ≠ fundamental(pm)]
    
    JuMP.@objective(pm.model, Max, sum(cmd))
end

## constraints
# ref bus
""
function constraint_ref_bus(pm::_PMs.AbstractIVRModel, n::Int, i::Int, vref)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, vr == vref)
    JuMP.@constraint(pm.model, vi == 0.0)
end

# bus
""
function constraint_voltage_rms_limit(pm::_PMs.AbstractIVRModel, i, vminrms, vmaxrms)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, vminrms^2 <= sum(vr.^2 + vi.^2)               )
    JuMP.@constraint(pm.model,              sum(vr.^2 + vi.^2)  <= vmaxrms^2 )
end
""
function constraint_voltage_rms_limit(pm::dHHC_SOC, i, vmaxrms, vmfund)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(vmaxrms^2 - vmfund^2); vcat(vr, vi)] in JuMP.SecondOrderCone())
end
""
function constraint_voltage_thd_limit(pm::_PMs.AbstractIVRModel, i, thdmax)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(vr[2:end].^2 + vi[2:end].^2) <= thdmax^2 * (vr[1]^2 + vi[1]^2))
end
""
function constraint_voltage_thd_limit(pm::dHHC_SOC, i, thdmax, vmfund)
    vr = [_PMs.var(pm, n, :vr, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    vi = [_PMs.var(pm, n, :vi, i) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [thdmax * vmfund; vcat(vr, vi)] in JuMP.SecondOrderCone())
end
""
function constraint_voltage_ihd_limit(pm::_PMs.AbstractIVRModel, n::Int, i, ihdmax)
    vr = [_PMs.var(pm, 1, :vr, i), _PMs.var(pm, n, :vr, i)] 
    vi = [_PMs.var(pm, 1, :vi, i), _PMs.var(pm, n, :vi, i)]

    JuMP.@constraint(pm.model, (vr[2]^2 + vi[2]^2) <= ihdmax^2 * (vr[1]^2 + vi[1]^2))
end
""
function constraint_voltage_ihd_limit(pm::dHHC_SOC, n::Int, i, ihdmax, vmfund)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    JuMP.@constraint(pm.model, [ihdmax * vmfund; vcat(vr, vi)] in JuMP.SecondOrderCone())
end
""
function constraint_current_balance(pm::_PMs.AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_xfmr, bus_filters, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    cr = _PMs.var(pm, n, :cr)
    ci = _PMs.var(pm, n, :ci)
    crt = _PMs.var(pm, n, :crt)
    cit = _PMs.var(pm, n, :cit)

    crf = _PMs.var(pm, n, :crf)
    cif = _PMs.var(pm, n, :cif)
    crg = _PMs.var(pm, n, :crg)
    cig = _PMs.var(pm, n, :cig)
    crd = _PMs.var(pm, n, :crd)
    cid = _PMs.var(pm, n, :cid)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crt[t] for t in bus_arcs_xfmr)
                                ==
                                sum(crf[f] for f in bus_filters)
                                + sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr 
                                + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cit[t] for t in bus_arcs_xfmr)
                                ==
                                sum(cif[f] for f in bus_filters)
                                + sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi 
                                - sum(bs for bs in values(bus_bs))*vr
                                )
end

# branch
""
function constraint_current_rms_limit(pm::_PMs.AbstractIVRModel, f_idx, t_idx, c_rating)
    crf =  [_PMs.var(pm, n, :cr, f_idx) for n in sorted_nw_ids(pm)]
    cif =  [_PMs.var(pm, n, :ci, f_idx) for n in sorted_nw_ids(pm)]

    crt =  [_PMs.var(pm, n, :cr, t_idx) for n in sorted_nw_ids(pm)]
    cit =  [_PMs.var(pm, n, :ci, t_idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crf.^2 + cif.^2) <= c_rating^2)
    JuMP.@constraint(pm.model, sum(crt.^2 + cit.^2) <= c_rating^2)
end
""
function constraint_current_rms_limit(pm::dHHC_SOC, f_idx, t_idx, c_rating, cm_fund_fr, cm_fund_to)
    crf =  [_PMs.var(pm, n, :cr, f_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cif =  [_PMs.var(pm, n, :ci, f_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    crt =  [_PMs.var(pm, n, :cr, t_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cit =  [_PMs.var(pm, n, :ci, t_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, sum(crf.^2 + cif.^2) <= c_rating^2 - cm_fund_fr^2)
    JuMP.@constraint(pm.model, sum(crt.^2 + cit.^2) <= c_rating^2 - cm_fund_to^2)
end

# filter
""
function constraint_active_filter(pm::_PMs.AbstractIVRModel, n::Int, f, i)
    vr = [_PMs.var(pm, nw, :vr, i) for nw in sorted_nw_ids(pm)]
    vi = [_PMs.var(pm, nw, :vi, i) for nw in sorted_nw_ids(pm)]
    
    crf = [_PMs.var(pm, nw, :crf, f) for nw in sorted_nw_ids(pm)]
    cif = [_PMs.var(pm, nw, :cif, f) for nw in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model,  vr[fundamental(pm)]*crf[fundamental(pm)] 
                                + vi[fundamental(pm)]*cif[fundamental(pm)] 
                                    == 
                                0.0
                    )
    JuMP.@constraint(pm.model,  sum(vr[n]*crf[n] + vi[n]*cif[n] 
                                        for n in 2:lastindex(vr)) 
                                    == 
                                0.0
                    )
end

# load
""
function constraint_load_constant_power(pm::_PMs.AbstractIVRModel, n::Int, l, i, pd, qd)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)
    crd = _PMs.var(pm, n, :crd, l)
    cid = _PMs.var(pm, n, :cid, l)

    JuMP.@constraint(pm.model, pd == vr*crd  + vi*cid)
    JuMP.@constraint(pm.model, qd == vi*crd  - vr*cid)
end
""
function constraint_load_constant_current(pm::_PMs.AbstractIVRModel, n::Int, l, mult)
    crd = _PMs.var(pm, n, :crd, l)
    cid = _PMs.var(pm, n, :cid, l)
    fund_crd = _PMs.var(pm, 1, :crd, l)
    fund_cid = _PMs.var(pm, 1, :cid, l)

    JuMP.@constraint(pm.model, crd == mult * fund_crd)
    JuMP.@constraint(pm.model, cid == mult * fund_cid)
end
""
function constraint_load_current_angle(pm::_PMs.AbstractIVRModel, n::Int, l, angmin, angmax)
    crd = _PMs.var(pm, n, :crd, l)
    cid = _PMs.var(pm, n, :cid, l)
    cmd = _PMs.var(pm, n, :cmd, l)

    JuMP.@constraint(pm.model, cmd * min(sin(angmin), sin(angmax)) <= cid)
    JuMP.@constraint(pm.model, cid <= cmd * max(sin(angmin), sin(angmax)))

    JuMP.@constraint(pm.model, cmd * min(cos(angmin), cos(angmax)) <= crd)
    JuMP.@constraint(pm.model, crd <= cmd * max(cos(angmin), cos(angmax)))

    JuMP.@constraint(pm.model, cmd^2 <= crd^2 + cid^2)
end
""
function constraint_load_current_angle(pm::dHHC_SOC, n::Int, l, angmin, angmax)
    crd = _PMs.var(pm, n, :crd, l)
    cid = _PMs.var(pm, n, :cid, l)
    cmd = _PMs.var(pm, n, :cmd, l)

    JuMP.@constraint(pm.model, cmd * min(sin(angmin), sin(angmax)) <= cid)
    JuMP.@constraint(pm.model, cid <= cmd * max(sin(angmin), sin(angmax)))

    JuMP.@constraint(pm.model, cmd * min(cos(angmin), cos(angmax)) <= crd)
    JuMP.@constraint(pm.model, crd <= cmd * max(cos(angmin), cos(angmax)))

    # NB: This changes the direction of the inequality.
    JuMP.@constraint(pm.model, [1/sqrt(2) * cmd; 1/sqrt(2) * cmd; vcat(crd, cid)] in JuMP.RotatedSecondOrderCone())
end

# xfmr
""
function constraint_transformer_core_magnetization(pm::_PMs.AbstractIVRModel, n::Int, t, int_a, int_b)
    cmrt = _PMs.var(pm, n, :cmrt, t)
    cmit = _PMs.var(pm, n, :cmit, t)

    et = reduce(vcat,[[_PMs.var(pm, nw, :ert, t), _PMs.var(pm, nw, :eit, t)] 
                                for nw in _PMs.ref(pm, n, :xfmr, t, "Hᴱ")])

    sym_exc_a = Symbol("exc_a_",n,"_",t)
    sym_exc_b = Symbol("exc_b_",n,"_",t)

    JuMP.register(pm.model, sym_exc_a, length(et), int_a; autodiff=true)
    JuMP.register(pm.model, sym_exc_b, length(et), int_b; autodiff=true)

    JuMP.add_nonlinear_constraint(pm.model, :($(cmrt) == $(sym_exc_a)($(et...))))
    JuMP.add_nonlinear_constraint(pm.model, :($(cmit) == $(sym_exc_b)($(et...))))
end
""
function constraint_transformer_core_magnetization(pm::_PMs.AbstractIVRModel, n::Int, t)
    cmrt = _PMs.var(pm, n, :cmrt, t)
    cmit = _PMs.var(pm, n, :cmit, t)

    JuMP.@constraint(pm.model, cmrt == 0.0)
    JuMP.@constraint(pm.model, cmit == 0.0)
end
""
function constraint_transformer_core_voltage_drop(pm::_PMs.AbstractIVRModel, n::Int, t, f_idx, xsc)
    ert = _PMs.var(pm, n, :ert, t)
    eit = _PMs.var(pm, n, :eit, t)

    vrt = _PMs.var(pm, n, :vrt, f_idx)
    vit = _PMs.var(pm, n, :vit, f_idx)

    csrt = _PMs.var(pm, n, :csrt, f_idx)
    csit = _PMs.var(pm, n, :csit, f_idx)
    
    JuMP.@constraint(pm.model, vrt == ert - xsc * csit)
    JuMP.@constraint(pm.model, vit == eit + xsc * csrt)
end
""
function constraint_transformer_core_voltage_balance(pm::_PMs.AbstractIVRModel, n::Int, t, t_idx, tr, ti)
    ert = _PMs.var(pm, n, :ert, t)
    eit = _PMs.var(pm, n, :eit, t)

    vrt = _PMs.var(pm, n, :vrt, t_idx)
    vit = _PMs.var(pm, n, :vit, t_idx)

    JuMP.@constraint(pm.model, vrt == tr * ert - ti * eit)
    JuMP.@constraint(pm.model, vit == tr * eit + ti * ert)
end
""
function constraint_transformer_core_current_balance(pm::_PMs.AbstractIVRModel, n::Int, t, f_idx, t_idx, tr, ti, rsh)
    cmrt = _PMs.var(pm, n, :cmrt, t)
    cmit = _PMs.var(pm, n, :cmit, t)

    csrt_fr = _PMs.var(pm, n, :csrt, f_idx)
    csit_fr = _PMs.var(pm, n, :csit, f_idx)

    csrt_to = _PMs.var(pm, n, :csrt, t_idx)
    csit_to = _PMs.var(pm, n, :csit, t_idx)

    ert = _PMs.var(pm, n, :ert, t)
    eit = _PMs.var(pm, n, :eit, t)

    JuMP.@constraint(pm.model,  csrt_fr 
                                + tr * csrt_to 
                                + ti * csit_to 
                                    == 
                                cmrt 
                                + ert / rsh
                    )
    JuMP.@constraint(pm.model,  csit_fr 
                                + tr * csit_to 
                                - ti * csrt_to
                                    == 
                                cmit 
                                + eit / rsh
                    )
end
""
function constraint_transformer_winding_config(pm::_PMs.AbstractIVRModel, n::Int, i, idx, r, re, xe, gnd)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    vrt = _PMs.var(pm, n, :vrt, idx)
    vit = _PMs.var(pm, n, :vit, idx)

    crt = _PMs.var(pm, n, :crt, idx)
    cit = _PMs.var(pm, n, :cit, idx)

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
function constraint_transformer_winding_current_balance(pm::_PMs.AbstractIVRModel, n::Int, idx, r, b_sh, g_sh, cnf)
    vrt = _PMs.var(pm, n, :vrt, idx)
    vit = _PMs.var(pm, n, :vit, idx)
    
    crt = _PMs.var(pm, n, :crt, idx)
    cit = _PMs.var(pm, n, :cit, idx)

    csrt = _PMs.var(pm, n, :csrt, idx)
    csit = _PMs.var(pm, n, :csit, idx)

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
function constraint_transformer_winding_current_rms_limit(pm::_PMs.AbstractIVRModel, idx, c_rating)
    crt =  [_PMs.var(pm, n, :crt, idx) for n in sorted_nw_ids(pm)]
    cit =  [_PMs.var(pm, n, :cit, idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crt.^2 + cit.^2) <= c_rating^2)
end
""
function constraint_transformer_winding_current_rms_limit(pm::dHHC_SOC, idx, c_rating, cm_fund)
    crt =  [_PMs.var(pm, n, :crt, idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cit =  [_PMs.var(pm, n, :cit, idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, sum(crt.^2 + cit.^2) <= c_rating^2 - cm_fund^2) 
end