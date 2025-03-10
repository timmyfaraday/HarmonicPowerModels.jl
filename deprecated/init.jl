################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Hakan Ergun, Tom Van Acker                                          #
################################################################################
# Changelog:                                                                   #
# v0.2.0 -  reviewed TVA                                                       #
# v0.2.1 -  update_hdata_with_fairness_principle_data: optimized to only build #
#           the model once and update for each load                            #
################################################################################

""
function update_hdata_with_fairness_principle_data!(hdata, model_type::Type, hhc_optimizer)
    if hdata["principle"] == "Kalai-Smorodinsky bargaining"
        if !haskey(hdata, "skip KS precalc") || hdata["skip KS precalc"] == false
            # make a deepcopy of the hdata
            hdata_temp = deepcopy(hdata)

            # set the principle to maximum efficiency
            hdata_temp["principle"] = "maximum efficiency"

            # instantiate with all loads
            pm = _PMs.instantiate_model(hdata_temp, model_type, build_hhc; ref_extensions=[ref_add_filter!, ref_add_xfmr!])
        
            # solve the maximum efficiency harmonic hosting capacity problem for 
            # each individual load
            for (l,load) in hdata["nw"]["1"]["load"]
                for (nw,ntw) in hdata_temp["nw"], (nl,load) in ntw["load"] 
                    nw_ = parse(Int, nw)
                    nl_ = parse(Int, nl)
                    c_rating = load["c_rating"]
                    if nw_ ≠ 1
                        if nl ≠ l
                            JuMP.set_lower_bound(_PMs.var(pm, nw_, :crd, nl_), 0.0)
                            JuMP.set_upper_bound(_PMs.var(pm, nw_, :crd, nl_), 0.0)
                            JuMP.set_lower_bound(_PMs.var(pm, nw_, :cid, nl_), 0.0)
                            JuMP.set_upper_bound(_PMs.var(pm, nw_, :cid, nl_), 0.0)
                            JuMP.set_lower_bound(_PMs.var(pm, nw_, :cmd, nl_), 0.0)
                            JuMP.set_upper_bound(_PMs.var(pm, nw_, :cmd, nl_), 0.0)
                        else 
                            JuMP.set_lower_bound(_PMs.var(pm, nw_, :crd, nl_), -c_rating)
                            JuMP.set_upper_bound(_PMs.var(pm, nw_, :crd, nl_),  c_rating)
                            JuMP.set_lower_bound(_PMs.var(pm, nw_, :cid, nl_), -c_rating)
                            JuMP.set_upper_bound(_PMs.var(pm, nw_, :cid, nl_),  c_rating)
                            JuMP.set_lower_bound(_PMs.var(pm, nw_, :cmd, nl_),  0.0)
                            JuMP.set_upper_bound(_PMs.var(pm, nw_, :cmd, nl_),  c_rating)
                        end 
                    end 
                end

                # solve the harmonic hosting capacity problem for the single load
                results_hhc_temp = _PMs.optimize_model!(pm, optimizer = hhc_optimizer)

                # write away the solution for each network
                for (nw,ntw) in results_hhc_temp["solution"]["nw"] 
                    if nw ≠ "1"
                        hdata["nw"]["$nw"]["load"]["$l"]["cmdmax"] = ntw["load"]["$l"]["cmd"]
                    end 
                end
            end
        end
    end
end

""
function update_hdata_with_fundamental_hpf_results!(hdata, model_type::Type, optimizer)
    # remove all but the fundamental network
    hpf_data = deepcopy(hdata)
    for n in keys(hpf_data["nw"])
        if n ≠ "1"
            delete!(hpf_data["nw"], n)
        end
    end

    # solve hpf problem for the fundamental harmonic only
    hpf_results = solve_hpf(hpf_data, model_type, optimizer)

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["vm"] = floor(hpf_results["solution"]["nw"]["1"]["bus"][i]["vm"], digits=10)
        bus["va"] = hpf_results["solution"]["nw"]["1"]["bus"][i]["va"]
    end
    for (b, branch) in hdata["nw"]["1"]["branch"]
        branch["cm_fr"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cr_fr"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["ci_fr"]^2)
        branch["cm_to"] = sqrt( hpf_results["solution"]["nw"]["1"]["branch"][b]["cr_to"]^2 
                                + hpf_results["solution"]["nw"]["1"]["branch"][b]["ci_to"]^2)
    end
    for (g, gen) in hdata["nw"]["1"]["gen"]
        gen["cm"] = sqrt(   hpf_results["solution"]["nw"]["1"]["gen"][g]["crg"]^2 
                            + hpf_results["solution"]["nw"]["1"]["gen"][g]["cig"]^2)
    end
    if haskey(hdata["nw"]["1"], "xfmr")
        for (x, xfmr) in hdata["nw"]["1"]["xfmr"]
            xfmr["ctm_fr"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["crx_fr"]^2 
                                    + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cix_fr"]^2)
            xfmr["ctm_to"] = sqrt(  hpf_results["solution"]["nw"]["1"]["xfmr"][x]["crx_to"]^2 
                                    + hpf_results["solution"]["nw"]["1"]["xfmr"][x]["cix_to"]^2)
        end
    end
end