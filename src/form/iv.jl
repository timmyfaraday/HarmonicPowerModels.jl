################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
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
function variable_branch_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# filter
""
function variable_filter_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_filter_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_filter_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# xfmr
""
function variable_xfmr_voltage(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_xfmr_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_xfmr_voltage_excitation_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_voltage_excitation_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_xfmr_current(pm::_PMs.AbstractIVRModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_xfmr_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_xfmr_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_xfmr_current_magnetizing_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_xfmr_current_magnetizing_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
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
    if nw ≠ fundamental(pm)
        variable_load_current_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    end
end

## objective
""
function objective_power_flow(pm::_PMs.AbstractIVRModel)
    JuMP.@objective(pm.model, Min, 0.0)
end
""
function objective_voltage_distortion_minimization(pm::_PMs.AbstractIVRModel) 
    bus_id = pm.data["bus_id"]

    vr = [_PMs.var(pm, n, :vr, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]
    vi = [_PMs.var(pm, n, :vi, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]

    JuMP.@objective(pm.model, Min, sum(vr.^2 + vi.^2))
end
""
function objective_maximum_hosting_capacity(pm::_PMs.AbstractIVRModel)
    # maximum efficiency
    if pm.data["principle"] == "maximum efficiency"
        cmd = [_PMs.var(pm, n, :cmd, l) for n in _PMs.nw_ids(pm) 
                                        for l in _PMs.ids(pm, :load, nw=n) 
                                        if n ≠ fundamental(pm)]
    
        JuMP.@objective(pm.model, Max, sum(cmd))
    end

    # absolute equality
    if pm.data["principle"] == "absolute equality"
        cmd = [_PMs.var(pm, n, :cmd, l) for n in _PMs.nw_ids(pm) 
                                        for l in _PMs.ids(pm, :load, nw=n) 
                                        if n ≠ fundamental(pm)]
    
        JuMP.@objective(pm.model, Max, sum(cmd)) 
    end

    # maximin
    if pm.data["principle"] == "maximin"
        cmh = [_PMs.var(pm, n, :cmh) for n in _PMs.nw_ids(pm)
                                     if n ≠ fundamental(pm)]

        JuMP.@objective(pm.model, Max, sum(cmh))
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["principle"] == "Kalai-Smorodinsky bargaining"
        fh = [_PMs.var(pm, n, :fh) for n in _PMs.nw_ids(pm)
                                   if n ≠ fundamental(pm)]

        JuMP.@objective(pm.model, Max, sum(fh))
    end
end

## constraints
# ref bus
""
function constraint_voltage_ref_bus(pm::_PMs.AbstractIVRModel, n::Int, i::Int, vref)
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
    crx = _PMs.var(pm, n, :crx)
    cix = _PMs.var(pm, n, :cix)

    crf = _PMs.var(pm, n, :crf)
    cif = _PMs.var(pm, n, :cif)
    crg = _PMs.var(pm, n, :crg)
    cig = _PMs.var(pm, n, :cig)
    crd = _PMs.var(pm, n, :crd)
    cid = _PMs.var(pm, n, :cid)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crx[t] for t in bus_arcs_xfmr)
                                ==
                                sum(crf[f] for f in bus_filters)
                                + sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr 
                                + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cix[t] for t in bus_arcs_xfmr)
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

    crx =  [_PMs.var(pm, n, :cr, t_idx) for n in sorted_nw_ids(pm)]
    cix =  [_PMs.var(pm, n, :ci, t_idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crf.^2 + cif.^2) <= c_rating^2)
    JuMP.@constraint(pm.model, sum(crx.^2 + cix.^2) <= c_rating^2)
end
""
function constraint_current_rms_limit(pm::dHHC_SOC, f_idx, t_idx, c_rating, cm_fund_fr, cm_fund_to)
    crf =  [_PMs.var(pm, n, :cr, f_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cif =  [_PMs.var(pm, n, :ci, f_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    crx =  [_PMs.var(pm, n, :cr, t_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cix =  [_PMs.var(pm, n, :ci, t_idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(c_rating^2 - cm_fund_fr^2); vcat(crf, cif)] in JuMP.SecondOrderCone())
    JuMP.@constraint(pm.model, [sqrt(c_rating^2 - cm_fund_to^2); vcat(crx, cix)] in JuMP.SecondOrderCone())
end

# fairness principle
""
function constraint_fairness_principle(pm::_PMs.AbstractIVRModel, n, load_ids)
    # maximum efficiency
    if pm.data["principle"] == "maximum efficiency"
        # no additional constraints
    end

    # absolute equality
    if pm.data["principle"] == "absolute equality"
        cmd = [_PMs.var(pm, n, :cmd, l) for l in load_ids]

        for l in load_ids[2:end]
            JuMP.@constraint(pm.model, cmd[first(load_ids)] .== cmd[l])
    end end

    # maximin
    if pm.data["principle"] == "maximin"
        cmh = _PMs.var(pm, n, :cmh)

        for l in load_ids
            cmd = _PMs.var(pm, n, :cmd, l)

            JuMP.@constraint(pm.model, cmh <= cmd)
    end end

    # Kalai-Smorodinsky bargaining
    if pm.data["principle"] == "Kalai-Smorodinsky bargaining"
        fh = _PMs.var(pm, n, :fh)

        for l in load_ids
            cmd = _PMs.var(pm, n, :cmd, l)
            cmdmax = _PMs.ref(pm, n, :load, l, "cmdmax")

            JuMP.@constraint(pm.model, cmd == fh * cmdmax)
    end end
end

# filter
""
function constraint_active_filter_current(pm::_PMs.AbstractIVRModel, f, i)
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
"" # needs work towards v0.2.1
function constraint_load_constant_power(pm::_PMs.AbstractIVRModel, n::Int, l, i, pd, qd)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)
    crd = _PMs.var(pm, n, :crd, l)
    cid = _PMs.var(pm, n, :cid, l)

    JuMP.@constraint(pm.model, pd == vr*crd  + vi*cid)
    JuMP.@constraint(pm.model, qd == vi*crd  - vr*cid)
end
"" # needs work towards v0.2.1
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
"" # needs work towards v0.2.1
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
"" # needs work towards v0.2.1
function constraint_load_constant_current(pm::_PMs.AbstractIVRModel, n::Int, l, mult)
    crd = _PMs.var(pm, n, :crd, l)
    cid = _PMs.var(pm, n, :cid, l)
    fund_crd = _PMs.var(pm, 1, :crd, l)
    fund_cid = _PMs.var(pm, 1, :cid, l)

    JuMP.@constraint(pm.model, crd == mult * fund_crd)
    JuMP.@constraint(pm.model, cid == mult * fund_cid)
end

# xfmr
""
function constraint_xfmr_core_magnetization(pm::_PMs.AbstractIVRModel, n::Int, x, int_a, int_b)
    cmrx = _PMs.var(pm, n, :cmrx, x)
    cmix = _PMs.var(pm, n, :cmix, x)

    et = reduce(vcat,[[_PMs.var(pm, nw, :erx, x), _PMs.var(pm, nw, :eix, x)] 
                                for nw in _PMs.ref(pm, n, :xfmr, x, "Hᴱ")])

    sym_exc_a = Symbol("exc_a_", n, "_", x)
    sym_exc_b = Symbol("exc_b_", n, "_", x)

    JuMP.register(pm.model, sym_exc_a, length(et), int_a; autodiff=true)
    JuMP.register(pm.model, sym_exc_b, length(et), int_b; autodiff=true)

    JuMP.add_nonlinear_constraint(pm.model, :($(cmrx) == $(sym_exc_a)($(et...))))
    JuMP.add_nonlinear_constraint(pm.model, :($(cmix) == $(sym_exc_b)($(et...))))
end
""
function constraint_xfmr_core_magnetization(pm::_PMs.AbstractIVRModel, n::Int, x)
    cmrx = _PMs.var(pm, n, :cmrx, x)
    cmix = _PMs.var(pm, n, :cmix, x)

    JuMP.@constraint(pm.model, cmrx == 0.0)
    JuMP.@constraint(pm.model, cmix == 0.0)
end
""
function constraint_xfmr_core_voltage_drop(pm::_PMs.AbstractIVRModel, n::Int, x, f_idx, xsc)
    erx = _PMs.var(pm, n, :erx, x)
    eix = _PMs.var(pm, n, :eix, x)

    vrx = _PMs.var(pm, n, :vrx, f_idx)
    vix = _PMs.var(pm, n, :vix, f_idx)

    csrx = _PMs.var(pm, n, :csrx, f_idx)
    csix = _PMs.var(pm, n, :csix, f_idx)
    
    JuMP.@constraint(pm.model, vrx == erx - xsc * csix)
    JuMP.@constraint(pm.model, vix == eix + xsc * csrx)
end
"""
first principles: uᵢ = tₓᵢⱼ * uⱼ
eₓₕ = tₓᵢⱼₕ * vₓⱼᵢₕ
eʳₓₕ + j eⁱₓₕ = (tʳₓᵢⱼₕ + j tⁱₓᵢⱼₕ) * (vʳₓⱼᵢₕ + j vⁱₓⱼᵢₕ)
eʳₓₕ + j eⁱₓₕ = tʳₓᵢⱼₕ vʳₓⱼᵢₕ + j tʳₓᵢⱼₕ vⁱₓⱼᵢₕ + j tⁱₓᵢⱼₕ vʳₓⱼᵢₕ + j² tⁱₓᵢⱼₕ vⁱₓⱼᵢₕ
eʳₓₕ + j eⁱₓₕ = tʳₓᵢⱼₕ vʳₓⱼᵢₕ + j tʳₓᵢⱼₕ vⁱₓⱼᵢₕ + j tⁱₓᵢⱼₕ vʳₓⱼᵢₕ - tⁱₓᵢⱼₕ vⁱₓⱼᵢₕ

Re: eʳₓₕ = tʳₓᵢⱼₕ vʳₓⱼᵢₕ - tⁱₓᵢⱼₕ vⁱₓⱼᵢₕ
Im: eⁱₓₕ = tʳₓᵢⱼₕ vⁱₓⱼᵢₕ + tⁱₓᵢⱼₕ vʳₓⱼᵢₕ
"""
function constraint_xfmr_core_voltage_phase_shift(pm::_PMs.AbstractIVRModel, n::Int, x, t_idx, tr, ti)
    erx = _PMs.var(pm, n, :erx, x)
    eix = _PMs.var(pm, n, :eix, x)

    vrx = _PMs.var(pm, n, :vrx, t_idx)
    vix = _PMs.var(pm, n, :vix, t_idx)

    JuMP.@constraint(pm.model, erx == tr * vrx - ti * vix)
    JuMP.@constraint(pm.model, eix == tr * vix + ti * vrx)
end
"""
first principles: conj(tₓᵢⱼ) * iₓᵢⱼ + iₓⱼᵢ = 0
conj(tₓᵢⱼₕ) * (iˢₓᵢⱼₕ - iᵐₓₕ - eₓₕ / rˢʰₓₕ) + iₓⱼᵢₕ = 0
(tʳₓᵢⱼₕ - j tⁱₓᵢⱼₕ) * (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ + j (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ)) + iˢ⁻ʳₓⱼᵢₕ + j iˢ⁻ⁱₓⱼᵢₕ = 0
tʳₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ) + j tʳₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ) - j tⁱₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ) - j² tⁱₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ) + iˢ⁻ʳₓⱼᵢₕ + j iˢ⁻ⁱₓⱼᵢₕ = 0
tʳₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ) + j tʳₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ) - j tⁱₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ) + tⁱₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ) + iˢ⁻ʳₓⱼᵢₕ + j iˢ⁻ⁱₓⱼᵢₕ = 0

Re: tʳₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ) + tⁱₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ) + iˢ⁻ʳₓⱼᵢₕ = 0
Im: tʳₓᵢⱼₕ (iˢ⁻ⁱₓᵢⱼₕ - iᵐ⁻ⁱₓₕ - eⁱₓₕ / rˢʰₓₕ) - tⁱₓᵢⱼₕ (iˢ⁻ʳₓᵢⱼₕ - iᵐ⁻ʳₓₕ - eʳₓₕ / rˢʰₓₕ) + iˢ⁻ⁱₓⱼᵢₕ = 0
"""
function constraint_xfmr_core_current_balance(pm::_PMs.AbstractIVRModel, n::Int, x, f_idx, t_idx, tr, ti, rsh)
    cmrx = _PMs.var(pm, n, :cmrx, x)
    cmix = _PMs.var(pm, n, :cmix, x)

    csrx_fr = _PMs.var(pm, n, :csrx, f_idx)
    csix_fr = _PMs.var(pm, n, :csix, f_idx)

    csrx_to = _PMs.var(pm, n, :csrx, t_idx)
    csix_to = _PMs.var(pm, n, :csix, t_idx)

    erx = _PMs.var(pm, n, :erx, x)
    eix = _PMs.var(pm, n, :eix, x)

    JuMP.@constraint(pm.model,  tr * (csrx_fr - cmrx - erx / rsh)
                                + ti * (csix_fr - cmix - eix / rsh)
                                + csrx_to 
                                    == 
                                0.0
                    )
    JuMP.@constraint(pm.model,  tr * (csix_fr - cmix - eix / rsh)
                                - ti * (csrx_fr - cmrx - erx / rsh)
                                + csix_to
                                    == 
                                0.0
                    )
end
""
function constraint_xfmr_winding_config(pm::_PMs.AbstractIVRModel, n::Int, i, idx, r, re, xe, gnd)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    vrx = _PMs.var(pm, n, :vrx, idx)
    vix = _PMs.var(pm, n, :vix, idx)

    crx = _PMs.var(pm, n, :crx, idx)
    cix = _PMs.var(pm, n, :cix, idx)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(n)
        JuMP.@constraint(pm.model, vrx == vr - r * crx)
        JuMP.@constraint(pm.model, vix == vi - r * cix)
    end

    # h ∈ 𝓗⁰, gnd == true -> cnf ∈ {Ye, Ze}
    if is_zero_sequence(n) && gnd == 1
        JuMP.@constraint(pm.model, vrx == vr - (r + 3re) * crx + 3xe * cix)
        JuMP.@constraint(pm.model, vix == vi - (r + 3re) * cix - 3xe * crx)
    end

    # h ∈ 𝓗⁰, gnd == false -> cnf ∈ {D, Y, Z}
    if is_zero_sequence(n) && gnd != 1
        JuMP.@constraint(pm.model, crx == 0)
        JuMP.@constraint(pm.model, cix == 0)
    end
end
""
function constraint_xfmr_winding_current_balance(pm::_PMs.AbstractIVRModel, n::Int, idx, r, b_sh, g_sh, cnf)
    vrx = _PMs.var(pm, n, :vrx, idx)
    vix = _PMs.var(pm, n, :vix, idx)
    
    crx = _PMs.var(pm, n, :crx, idx)
    cix = _PMs.var(pm, n, :cix, idx)

    csrx = _PMs.var(pm, n, :csrx, idx)
    csix = _PMs.var(pm, n, :csix, idx)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(n)
        JuMP.@constraint(pm.model, crx == csrx - g_sh * vrx + b_sh * vix)
        JuMP.@constraint(pm.model, cix == csix - g_sh * vix - b_sh * vrx)
    end

    # h ∈ 𝓗⁰, cnf ∈ {Y(e), Z(e)}
    if is_zero_sequence(n) && cnf in ['Y','Z']
        JuMP.@constraint(pm.model, crx == csrx - g_sh * vrx + b_sh * vix)
        JuMP.@constraint(pm.model, cix == csix - g_sh * vix - b_sh * vrx)
    end

    # h ∈ 𝓗⁰, cnf ∈ {D}
    if is_zero_sequence(n) && cnf in ['D'] && r ≠ 0.0
        JuMP.@constraint(pm.model, crx == csrx - g_sh * vrx + b_sh * vix - vrx / r)
        JuMP.@constraint(pm.model, cix == csix - g_sh * vix - b_sh * vrx - vix / r)
    end
end
""
function constraint_xfmr_current_rms_limit(pm::_PMs.AbstractIVRModel, idx, c_rating)
    crx =  [_PMs.var(pm, n, :crx, idx) for n in sorted_nw_ids(pm)]
    cix =  [_PMs.var(pm, n, :cix, idx) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crx.^2 + cix.^2) <= c_rating^2)
end
""
function constraint_xfmr_current_rms_limit(pm::dHHC_SOC, idx, c_rating, cm_fund)
    crx =  [_PMs.var(pm, n, :crx, idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cix =  [_PMs.var(pm, n, :cix, idx) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(c_rating^2 - cm_fund^2); vcat(crx, cix)] in JuMP.SecondOrderCone())
end