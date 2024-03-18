################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.1 -  add function to solve hhc problem iteratively using cuts for the   #
#           load current angle (TVA)                                           #
################################################################################

function solve_model_with_current_angle_cuts(hdata, model_type::Type, optimizer; 
                                             ref_extensions=[],
                                             solution_processors=[], kwargs...)
    # set iteration and start timer
    itr = 1
    obj = Float64[]
    start_time = time()
    
    # solve first iteration
    pm  = _PMs.instantiate_model(hdata, model_type, build_hhc,
                                    ref_extensions=ref_extensions)
    res = _PMs.optimize_model!(pm, optimizer=optimizer, 
                                    solution_processors=solution_processors)
    push!(obj, res["objective"])

    # build arrays for aref and asol in pm.data
    for n in _PMs.nw_ids(pm), l in _PMs.ids(pm, :load, nw=n) if n ≠ 1
        pm.data["nw"]["$n"]["load"]["$l"]["aref"] = 
                [hdata["nw"]["$n"]["load"]["$l"]["angle_ref"]]
        pm.data["nw"]["$n"]["load"]["$l"]["arng"] = 
                [hdata["nw"]["$n"]["load"]["$l"]["angle_rng"]]
        pm.data["nw"]["$n"]["load"]["$l"]["asol"] = 
                [res["solution"]["nw"]["$n"]["load"]["$l"]["ca"]]
    end end

    # update the linear inequality ar * crd + ai * cid >= cmd
    while itr < 20
        for n in _PMs.nw_ids(pm), l in _PMs.ids(pm, :load, nw=n) if n ≠ 1
            # get the input and solution of the last iteration
            aref = last(_PMs.ref(pm, n, :load, l, "aref"))
            arng = last(_PMs.ref(pm, n, :load, l, "arng"))
            asol = last(_PMs.ref(pm, n, :load, l, "asol"))
            
            # update the angle range
            arng /= 2
            push!(_PMs.ref(pm, n, :load, l, "arng"), arng)

            # update the reference angle
            adif = aref - asol
            if isapprox(cos(adif), 1.0; atol=1e-10)
                aref_itr = aref
            elseif !isapprox(cos(adif), 1.0; atol=1e-10) && sin(adif) < 0.0
                aref_itr = aref + arng / 2
            elseif !isapprox(cos(adif), 1.0; atol=1e-10) && sin(adif) > 0.0
                aref_itr = aref - arng / 2
            end
            push!(_PMs.ref(pm, n, :load, l, "aref"), aref_itr)
            
            # update the coefficients of the linear inequality
            amin = aref_itr - arng / 2
            amax = aref_itr + arng / 2

            crd = _PMs.var(pm, n, :crd, l)
            cid = _PMs.var(pm, n, :cid, l)
            cstr = JuMP.constraint_by_name(pm.model, "cstr_$(n)_$(l)")

            JuMP.set_normalized_coefficient(cstr, crd, (sin(amin) - sin(amax)) / sin(amin-amax))
		    JuMP.set_normalized_coefficient(cstr, cid, (cos(amax) - cos(amin)) / sin(amin-amax))
        end end    

        # update the iteration count
        itr += 1

        # solve next iteration
        res = _PMs.optimize_model!(pm, optimizer=optimizer, 
                                    solution_processors=solution_processors)
        push!(obj, res["objective"])

        # update the solution angle
        for n in _PMs.nw_ids(pm), l in _PMs.ids(pm, :load, nw=n) if n ≠ 1
            asol = res["solution"]["nw"]["$n"]["load"]["$l"]["ca"]
            push!(_PMs.ref(pm, n, :load, l, "asol"), asol)
        end end
    end

    # set the solution
    for n in _PMs.nw_ids(pm), l in _PMs.ids(pm, :load, nw=n) if n ≠ 1
        res["solution"]["nw"]["$n"]["load"]["$l"]["aref"] = _PMs.ref(pm, n, :load, l, "aref")
        res["solution"]["nw"]["$n"]["load"]["$l"]["arng"] = _PMs.ref(pm, n, :load, l, "arng")
        res["solution"]["nw"]["$n"]["load"]["$l"]["asol"] = _PMs.ref(pm, n, :load, l, "asol")
    end end

    # set solve properties
    res["objective"]  = obj
    res["solve_time"] = time() - start_time    
    res["iterations"] = itr

    return res
end