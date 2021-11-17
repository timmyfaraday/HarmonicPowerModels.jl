
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

    expression_transformer_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    # expression_transformer_series_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_transformer_excitation_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
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
function constraint_transformer_core_excitation(pm::AbstractIVRModel, n::Int, t, int_a, int_b, grad_a, grad_b)
    cert = var(pm, n, :cert, t)
    ceit = var(pm, n, :ceit, t)

    voltage_harmonics_ntws = _PMs.ref(pm, n, :xfmr, t, "voltage_harmonics_ntws")

    et = reduce(vcat,[[var(pm, nw, :ert, t),var(pm, nw, :eit, t)] 
                       for nw in voltage_harmonics_ntws])

    sym_exc_a = Symbol("exc_a_",n,"_",t)
    sym_exc_b = Symbol("exc_b_",n,"_",t)

    JuMP.register(pm.model, sym_exc_a, length(et), int_a, grad_a)
    JuMP.register(pm.model, sym_exc_b, length(et), int_b, grad_b)

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

    JuMP.@constraint(pm.model, csrt_fr + tr * csrt_to + ti * csit_to 
                                == cert 
                    )
    JuMP.@constraint(pm.model, csit_fr + tr * csit_to - ti * csrt_to
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
function constraint_voltage_magnitude_rms(pm::AbstractIVRModel, i, vminrms, vmaxrms)
    w = [var(pm, nw, :w, i) for nw in _PMs.nw_ids(pm)]
    @assert vminrms>0
    @assert vmaxrms>vminrms

    JuMP.@constraint(pm.model, vminrms^2 <= sum(w)               )
    JuMP.@constraint(pm.model,              sum(w)  <= vmaxrms^2 )
end


""
function constraint_voltage_thd(pm::AbstractIVRModel, i, fundamental, thdmax)
    harmonics = Set(_PMs.nw_ids(pm))
    nonfundamentalharmonics = setdiff(harmonics, [fundamental])
    w = [var(pm, nw, :w, i) for nw in nonfundamentalharmonics]
    wfun = var(pm, fundamental, :w, i)

    JuMP.@constraint(pm.model, sum(w) <= thdmax^2*(wfun))
end

function constraint_voltage_harmonics_relative_magnitude(pm::AbstractIVRModel, n::Int, i, rm, fundamental)
    w  = var(pm, n, :w, i)
    wfun  = var(pm, fundamental, :w, i)

    JuMP.@constraint(pm.model, w <= rm*wfun)
end

function constraint_load_constant_power(pm::AbstractIVRModel, n::Int, i, bus, pd, qd)
    vr = var(pm, n, :vr, bus)
    vi = var(pm, n, :vi, bus)
    crd = var(pm, n, :crd, i)
    cid = var(pm, n, :cid, i)

    JuMP.@constraint(pm.model, pd == vr*crd  + vi*cid)
    JuMP.@constraint(pm.model, qd == vi*crd  - vr*cid)

end


function constraint_load_constant_current(pm::AbstractIVRModel, n::Int, i, bus, multiplier)
    crd = var(pm, n, :crd, i)
    cid = var(pm, n, :cid, i)
    fundamental = 1 
    crd_fund = var(pm, fundamental, :crd, i)
    cid_fund = var(pm, fundamental, :cid, i)

    JuMP.@constraint(pm.model, crd == multiplier * crd_fund)
    JuMP.@constraint(pm.model, cid == multiplier * cid_fund)
end


function constraint_vm_auxiliary_variable(pm::AbstractIVRModel, n::Int, i)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)
    w  = var(pm, n, :w, i)
    #TODO choose whether this is relaxation or exact
    JuMP.@constraint(pm.model, w >= vr^2  + vi^2)
end


"reference bus angle constraint"
function constraint_ref_bus(pm::AbstractIVRModel, n::Int, i::Int)
    if n == 1 #fundamental frequency, fix reference angle
        JuMP.@constraint(pm.model, var(pm, n, :vi)[i] == 0.0)
        #TODO you should be able to set this more freely, but it helps a lot with stability, keeping it in for now.
        JuMP.@constraint(pm.model, var(pm, n, :vr)[i] == 1.0)
    else #fix harmonic voltage at reference bus to 0+j0
        JuMP.@constraint(pm.model, var(pm, n, :vi)[i] == 0.0)
        JuMP.@constraint(pm.model, var(pm, n, :vr)[i] == 0.0)
    end
end

function constraint_current_limit_rms(pm::AbstractIVRModel, f_idx, c_rating)
    (l, f_bus, t_bus) = f_idx
    t_idx = (l, t_bus, f_bus)

    crf =  [var(pm, n, :cr, f_idx) for n in _PMs.nw_ids(pm)]
    cif =  [var(pm, n, :ci, f_idx) for n in _PMs.nw_ids(pm)]

    crt =  [var(pm, n, :cr, t_idx) for n in _PMs.nw_ids(pm)]
    cit =  [var(pm, n, :ci, t_idx) for n in _PMs.nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(crf^2 + cif^2) <= c_rating^2)
    JuMP.@constraint(pm.model, sum(crt^2 + cit^2) <= c_rating^2)
end


function constraint_active_filter(pm::AbstractIVRModel, i, fundamental)
    pgfun = var(pm, fundamental, :pg, i)
    pg =  [var(pm, n, :pg, i) for n in _PMs.nw_ids(pm)]

    JuMP.@NLconstraint(pm.model, pgfun == 0)
    JuMP.@NLconstraint(pm.model, sum(pg[n] for n in _PMs.nw_ids(pm)) == 0)
end


""
function objective_current_distortion_minimization(pm::AbstractIVRModel; gen_id=1, fundamental=1)
    harmonics = Set(_PMs.nw_ids(pm))
    nonfundamentalharmonics = setdiff(harmonics, [fundamental])

    crg = [var(pm, n, :crg, gen_id) for n in nonfundamentalharmonics]
    cig = [var(pm, n, :cig, gen_id) for n in nonfundamentalharmonics]

    pg = var(pm, fundamental, :pg, gen_id)

    JuMP.@NLconstraint(pm.model, pg <= 1.001*0.643386)

    #minimize magnitude of nonfundamental harmonics
    JuMP.@objective(pm.model, Min, sum(crg.^2 + cig.^2))
end


""
function objective_voltage_distortion_minimization(pm::AbstractIVRModel; bus_id=6, fundamental=1, gen_id=1)
    harmonics = Set(_PMs.nw_ids(pm))
    nonfundamentalharmonics = setdiff(harmonics, [fundamental])

    vr = [var(pm, n, :vr, bus_id) for n in nonfundamentalharmonics]
    vi = [var(pm, n, :vi, bus_id) for n in nonfundamentalharmonics]

    pg = var(pm, fundamental, :pg, gen_id)

    JuMP.@NLconstraint(pm.model, pg <= 1.001*0.643386)

    #minimize magnitude of nonfundamental harmonics
    JuMP.@objective(pm.model, Min, sum(vr.^2 + vi.^2))
end


#solution_processors=[sol_data_model!]
function sol_data_model!(pm::AbstractIVRModel, solution::Dict)
    _PMs.apply_pm!(_sol_data_model_ivr!, solution)
end


""
function _sol_data_model_ivr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vr") && haskey(bus, "vi")
                bus["vm"] = hypot(bus["vr"], bus["vi"])
                bus["va"] = atan(bus["vi"], bus["vr"])*180/pi
            end
        end
    end

    if haskey(solution, "branch")
        for (i, branch) in solution["branch"]
            if haskey(branch, "pf") && haskey(branch, "pt")
                branch["ploss"] = branch["pf"] + branch["pt"]
            end
            if haskey(branch, "qf") && haskey(branch, "qt")
                branch["qloss"] = branch["qf"] + branch["qt"]
            end
        end
    end

    if haskey(solution, "xfmr")
        for (i, xfmr) in solution["xfmr"]
            if haskey(xfmr, "pt_fr") && haskey(xfmr, "pt_to")
                xfmr["ptloss"] = xfmr["pt_fr"] + xfmr["pt_to"]
            end
            if haskey(xfmr, "qt_fr") && haskey(xfmr, "qt_to")
                xfmr["qtloss"] = xfmr["qt_fr"] + xfmr["qt_to"]
            end
        end
    end

    if haskey(solution, "gen")
        for (i, gen) in solution["gen"]
            if haskey(gen, "crg") && haskey(gen, "cig")
                gen["cm"] = hypot(gen["crg"], gen["cig"])
            end
        end
    end

    
    if haskey(solution, "load")
        for (i, load) in solution["load"]
            if haskey(load, "crd") && haskey(load, "cid")
                load["cm"] = hypot(load["crd"], load["cid"])
            end
        end
    end
end

function append_indicators!(result, hdata)
    solu = result["solution"]
    fundamental = "1"
    harmonics = Set(n for (n,nw) in solu["nw"])
    nonfundamentalharmonics = setdiff(harmonics, [fundamental])

    for (i,bus) in solu["nw"][fundamental]["bus"]
        vfun = solu["nw"][fundamental]["bus"][i]["vr"] + im*solu["nw"][fundamental]["bus"][i]["vi"]
        v   = [solu["nw"][n]["bus"][i]["vr"] + im*solu["nw"][n]["bus"][i]["vi"]  for n in harmonics]
        vnonfun = [solu["nw"][n]["bus"][i]["vr"] + im*solu["nw"][n]["bus"][i]["vi"]  for n in nonfundamentalharmonics]
        
        rms = sqrt(sum(abs.(v).^2))
        thd = sqrt(sum(abs.(vnonfun).^2))/abs(vfun)

        delete!(solu["nw"][fundamental]["bus"][i], "w")
        delete!(solu["nw"][fundamental]["bus"][i], "vr")
        delete!(solu["nw"][fundamental]["bus"][i], "vi")

        solu["nw"][fundamental]["bus"][i]["vrms"] = rms
        solu["nw"][fundamental]["bus"][i]["vthd"] = thd
        solu["nw"][fundamental]["bus"][i]["vmaxrms"] = hdata["nw"][fundamental]["bus"][i]["vmaxrms"]
        solu["nw"][fundamental]["bus"][i]["vminrms"] = hdata["nw"][fundamental]["bus"][i]["vminrms"]
        solu["nw"][fundamental]["bus"][i]["thdmax"] = hdata["nw"][fundamental]["bus"][i]["thdmax"]

    end

    for (i,gen) in solu["nw"][fundamental]["gen"]
        gencost = hdata["nw"][fundamental]["gen"][i]["cost"]
        pg = gen["pg"]
        if length(gencost) >0
            gen["totcost"] = sum(gencost[c]*(pg)^c for c in 1:length(gencost))
        end

        cfun = solu["nw"][fundamental]["gen"][i]["crg"] + im*solu["nw"][fundamental]["gen"][i]["cig"]
        call   = [solu["nw"][n]["gen"][i]["crg"] + im*solu["nw"][n]["gen"][i]["cig"]  for n in harmonics]
        cnonfun = [solu["nw"][n]["gen"][i]["crg"] + im*solu["nw"][n]["gen"][i]["cig"]  for n in nonfundamentalharmonics]

        rms = sqrt(sum(abs.(call).^2))
        thd = sqrt(sum(abs.(cnonfun).^2))/abs(cfun)
        
        solu["nw"][fundamental]["gen"][i]["crms"] = rms
        solu["nw"][fundamental]["gen"][i]["cthd"] = thd
        solu["nw"][fundamental]["gen"][i]["c_rating"] = hdata["nw"][fundamental]["gen"][i]["c_rating"]
        solu["nw"][fundamental]["gen"][i]["pmax"] = hdata["nw"][fundamental]["gen"][i]["pmax"]
    end

    if haskey(solu["nw"][fundamental], "load")
        for (i,load) in solu["nw"][fundamental]["load"]

            cfun = solu["nw"][fundamental]["load"][i]["crd"] + im*solu["nw"][fundamental]["load"][i]["cid"]
            call   = [solu["nw"][n]["load"][i]["crd"] + im*solu["nw"][n]["load"][i]["cid"]  for n in harmonics]
            cnonfun = [solu["nw"][n]["load"][i]["crd"] + im*solu["nw"][n]["load"][i]["cid"]  for n in nonfundamentalharmonics]

            rms = sqrt(sum(abs.(call).^2))
            thd = sqrt(sum(abs.(cnonfun).^2))/abs(cfun)
            
            solu["nw"][fundamental]["load"][i]["crms"] = rms
            solu["nw"][fundamental]["load"][i]["cthd"] = thd
        end 
    end
end