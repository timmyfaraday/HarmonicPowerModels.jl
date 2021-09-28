
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
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_arcs_xfmr, bus_gens, bus_loads, bus_gs, bus_bs)
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
    crd = var(pm, n, :crd)
    cid = var(pm, n, :cid)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crdc[d] for d in bus_arcs_dc)
                                + sum(crt[t] for t in bus_arcs_xfmr)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                + sum(cit[t] for t in bus_arcs_xfmr)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )

end

""
function constraint_transformer_core_excitation(pm::AbstractIVRModel, n::Int, t)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    JuMP.@constraint(pm.model, cert == 0.0)
    JuMP.@constraint(pm.model, ceit == 0.0)
end

""
function constraint_transformer_core_excitation(pm::AbstractIVRModel, n::Int, t, int_a, int_b)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    voltage_harmonics_ntws = _PMs.ref(pm, nw, :xfmr, t, "voltage_harmonics_ntws")

    et = reduce(vcat,[[var(pm, nw, :ert, t),var(pm, nw, :eit, t)] 
                       for nw in voltage_harmonics_ntws])

    sym_exc_a = Symbol("exc_a_",n,"_",t)
    sym_exc_b = Symbol("exc_b_",n,"_",t)

    JuMP.register(pm.model, sym_exc_a, length(et), int_a; autodiff=true)
    JuMP.register(pm.model, sym_exc_b, length(et), int_b; autodiff=true)

    JuMP.add_NL_constraint(pm.model, :($(cert) == $(sym_exc_a)($(et...))))
    JuMP.add_NL_constraint(pm.model, :($(ceit) == $(sym_exc_b)($(et...))))
end

""
function constraint_transformer_core_voltage_drop(pm::AbstractIVRModel, n::Int, t, f_idx, xsc)
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
function constraint_transformer_core_voltage_balance(pm::AbstractIVRModel, n::Int, t, t_idx, tr, ti)
    ert = var(pm, n, :ert, t)
    eit = var(pm, n, :eit, t)

    vrt = var(pm, n, :vrt, t_idx)
    vit = var(pm, n, :vit, t_idx)

    JuMP.@constraint(pm.model, vrt == tr * ert - ti * eit)
    JuMP.@constraint(pm.model, vit == tr * eit + ti * ert)
end

""
function constraint_transformer_core_current_balance(pm::AbstractIVRModel, n::Int, t, f_idx, t_idx, tr, ti)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    csrt_fr = var(pm, n, :csrt, f_idx)
    csit_fr = var(pm, n, :csit, f_idx)

    csrt_to = var(pm, n, :csrt, t_idx)
    csit_to = var(pm, n, :csit, t_idx)

    JuMP.@constraint(pm.model, csrt_fr + tr * csrt_to - ti * csit_to 
                                == cert 
                    )
    JuMP.@constraint(pm.model, csit_fr + tr * csit_to + ti * csrt_to
                                == ceit
                    )
end

""
function constraint_transformer_winding_config(pm::AbstractIVRModel, n::Int, nh, i, idx, r, re, xe, gnd)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    vrt = var(pm, n, :vrt, idx)
    vit = var(pm, n, :vit, idx)

    crt = var(pm, n, :crt, idx)
    cit = var(pm, n, :cit, idx)

    @show nh, i, idx, r, re, xe, gnd, is_zero_sequence(nh)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(nh)
        JuMP.@constraint(pm.model, vrt == vr - r * crt)
        JuMP.@constraint(pm.model, vit == vi - r * cit)
    end

    # h ∈ 𝓗⁰, gnd == true -> cnf ∈ {Ye, Ze}
    if is_zero_sequence(nh) && gnd == 1
        JuMP.@constraint(pm.model, vrt == vr - (r + 3re) * crt + 3xe * cit)
        JuMP.@constraint(pm.model, vit == vi - (r + 3re) * cit - 3xe * crt)
    end

    # h ∈ 𝓗⁰, gnd == false -> cnf ∈ {D, Y, Z}
    if is_zero_sequence(nh) && gnd != 1
        JuMP.@constraint(pm.model, crt == 0)
        JuMP.@constraint(pm.model, cit == 0)
    end
end

""
function constraint_transformer_winding_current_balance(pm::AbstractIVRModel, n::Int, nh, idx, r, b_sh, g_sh, cnf)
    vrt = var(pm, n, :vrt, idx)
    vit = var(pm, n, :vit, idx)
    
    crt = var(pm, n, :crt, idx)
    cit = var(pm, n, :cit, idx)

    csrt = var(pm, n, :csrt, idx)
    csit = var(pm, n, :csit, idx)

    # h ∈ 𝓗⁺ ⋃ 𝓗⁻
    if !is_zero_sequence(nh)
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit)
        JuMP.@constraint(pm.model, cit == csit - g_sh * vit - b_sh * vrt)
    end

    # h ∈ 𝓗⁰, cnf ∈ {Y(e), Z(e)}
    if is_zero_sequence(nh) && cnf in ['Y','Z']
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit)
        JuMP.@constraint(pm.model, cit == csit - g_sh * vit - b_sh * vrt)
    end

    # h ∈ 𝓗⁰, cnf ∈ {D}
    if is_zero_sequence(nh) && cnf in ['D']
        JuMP.@constraint(pm.model, crt == csrt - g_sh * vrt + b_sh * vit - vrt / r)
        JuMP.@constraint(pm.model, crt == csit - g_sh * vit - b_sh * vrt - vit / r)
    end
end

""
function constraint_voltage_magnitude_rms(pm::AbstractIVRModel, i, vminrms, vmaxrms, nharmonics)
    vr = [var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm)]
    vi = [var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm)]
    @assert vminrms>0
    @assert vmaxrms>vminrms

    JuMP.@constraint(pm.model, vminrms^2 <= sum(vr.^2 + vi.^2)/nharmonics               )
    JuMP.@constraint(pm.model,              sum(vr.^2 + vi.^2)/nharmonics  <= vmaxrms^2 )
end


""
function constraint_voltage_thd(pm::AbstractIVRModel, i, fundamental, thdmax)
    harmonics = Set(_PMs.nw_ids(pm))
    nonfundamentalharmonics = setdiff(harmonics, fundamental)
    vr = [var(pm, nw, :vr, i) for nw in nonfundamentalharmonics]
    vi = [var(pm, nw, :vi, i) for nw in nonfundamentalharmonics]

    vrfun = var(pm, fundamental, :vr, fundamental)
    vifun = var(pm, fundamental, :vi, fundamental)

    JuMP.@constraint(pm.model, sum(vr.^2 + vi.^2) <= thdmax^2*(vrfun^2 + vifun^2))
end


function constraint_load_constant_power(pm::AbstractIVRModel, n::Int, i, bus, pd, qd)
    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    crd = var(pm, n, :crd, i)
    cid = var(pm, n, :cid, i)
    ccmd = var(pm, n, :ccmd, i)

    JuMP.@constraint(pm.model, pd == vr*crd  + vi*cid)
    JuMP.@constraint(pm.model, qd == vi*crd  - vr*cid)

    JuMP.@constraint(pm.model, ccmd == crd^2  + cid^2)
end


function constraint_load_constant_current(pm::AbstractIVRModel, n::Int, i, bus, multiplier)
    crd = var(pm, n, :crd, i)
    cid = var(pm, n, :cid, i)
    crd_fund = var(pm, 1, :crd, i)
    cid_fund = var(pm, 1, :cid, i)

    JuMP.@constraint(pm.model, crd == multiplier * crd_fund)
    JuMP.@constraint(pm.model, cid == multiplier * cid_fund)

    
    # ccmd = var(pm, n, :ccmd, i)
    #current magnitude squared of fundamental
    # ccmdfundamental = var(pm, 1, :ccmd, i)

    # JuMP.@constraint(pm.model, ccmd == crd^2  + cid^2)
    # JuMP.@constraint(pm.model, ccmd == (multiplier)^2 * ccmdfundamental)
end


function constraint_vm_auxiliary_variable(pm::AbstractIVRModel, n::Int, i)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
    vm = var(pm, n, :vm, i)

    JuMP.@constraint(pm.model, vm^2 == vr^2  + vi^2)
end


"reference bus angle constraint"
function constraint_ref_bus(pm::AbstractIVRModel, n::Int, i::Int)
    if n == 1 #fundamental frequency, fix reference angle
        JuMP.@constraint(pm.model, var(pm, n, :vi)[i] == 0.0)
        JuMP.@constraint(pm.model, var(pm, n, :vr)[i] >= 0.9)
    else #fix harmonic voltage at reference bus to 0+j0
        JuMP.@constraint(pm.model, var(pm, n, :vi)[i] == 0.0)
        JuMP.@constraint(pm.model, var(pm, n, :vr)[i] == 0.0)
    end
end

function constraint_current_limit_rms(pm::AbstractIVRModel, f_idx, c_rating, nharmonics)
    (l, f_bus, t_bus) = f_idx
    t_idx = (l, t_bus, f_bus)

    crf =  [var(pm, n, :cr, f_idx) for n in _PMs.nw_ids(pm)]
    cif =  [var(pm, n, :ci, f_idx) for n in _PMs.nw_ids(pm)]

    crt =  [var(pm, n, :cr, t_idx) for n in _PMs.nw_ids(pm)]
    cit =  [var(pm, n, :ci, t_idx) for n in _PMs.nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crf^2 + cif^2)/nharmonics <= c_rating^2)
    JuMP.@constraint(pm.model, sum(crt^2 + cit^2)/nharmonics <= c_rating^2)
end
