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
        # make a deepcopy of the hdata
        hdata_temp = deepcopy(hdata)

        # set the principle to maximum efficiency
        hdata_temp["principle"] = "maximum efficiency"

        # instantiate with all loads
        pm = _PMs.instantiate_model(hdata_temp, model_type, build_hhc)
    
        # solve the maximum efficiency harmonic hosting capacity problem for 
        # each individual load
        for (l,load) in hdata["nw"]["1"]["load"]
            for (nw,ntw) in hdata_temp["nw"], (nl,load) in ntw["load"] 
                c_rating = load["c_rating"]
                if nl ≠ l
                    JuMP.set_lower_bound(_PMs.var(pm, nw, :crd, nl), 0.0)
                    JuMP.set_upper_bound(_PMs.var(pm, nw, :crd, nl), 0.0)
                    JuMP.set_lower_bound(_PMs.var(pm, nw, :cid, nl), 0.0)
                    JuMP.set_upper_bound(_PMs.var(pm, nw, :cid, nl), 0.0)
                    JuMP.set_lower_bound(_PMs.var(pm, nw, :cmd, nl), 0.0)
                    JuMP.set_upper_bound(_PMs.var(pm, nw, :cmd, nl), 0.0)
                else 
                    JuMP.set_lower_bound(_PMs.var(pm, nw, :crd, nl), -c_rating)
                    JuMP.set_upper_bound(_PMs.var(pm, nw, :crd, nl),  c_rating)
                    JuMP.set_lower_bound(_PMs.var(pm, nw, :cid, nl), -c_rating)
                    JuMP.set_upper_bound(_PMs.var(pm, nw, :cid, nl),  c_rating)
                    JuMP.set_lower_bound(_PMs.var(pm, nw, :cmd, nl),  0.0)
                    JuMP.set_upper_bound(_PMs.var(pm, nw, :cmd, nl),  c_rating)
            end end

            # solve the harmonic hosting capacity problem for the single load
            results_hhc_temp = optimize_model!(pm, hhc_optimizer)

            # write away the solution for each network
            for (nw,ntw) in results_hhc_temp["solution"]["nw"] if nw ≠ "1"
                hdata["nw"]["$nw"]["load"]["$l"]["cmdmax"] = ntw["load"]["$l"]["cmd"]
            end end
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
    status = hpf_results["termination_status"]
    status in (_MOI.OPTIMAL, _MOI.LOCALLY_SOLVED, _MOI.ALMOST_LOCALLY_SOLVED) ||
        error("update_hdata_with_fundamental_hpf_results!: the fundamental " *
              "harmonic power flow did not solve (termination status $status); " *
              "the hosting capacity would otherwise use a flat-start profile.")

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["vm"] = floor(hpf_results["solution"]["nw"]["1"]["bus"][i]["vm"], digits=10)
        bus["va"] = hpf_results["solution"]["nw"]["1"]["bus"][i]["va"]
    end
    update_hdata_with_angle_offsets!(hdata)
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

"""
    update_hdata_with_angle_offsets!(hdata)

Set `load["angle_offset"]` [deg] for every load in every harmonic network,
following `hdata["angle_reference"]`:

- `"global"` (default): offset 0 — harmonic angles are given in the network
  frame, i.e., relative to the reference bus voltage.
- `"local"`: harmonic angles are given relative to the unit's own fundamental
  bus voltage, as a converter model produces them. A time shift of the
  fundamental by θ_i rotates the h-th harmonic phasor by h·θ_i, so the angle in
  the network frame gains `angle_offset = h·θ_i`, with θ_i the fundamental
  voltage angle at the unit's bus.

Requires the fundamental bus voltage angles `bus["va"]` [deg], as set by
`update_hdata_with_fundamental_hpf_results!`.
"""
function update_hdata_with_angle_offsets!(hdata)
    mode = get(hdata, "angle_reference", "global")
    mode in ("global", "local") || error(
        "update_hdata_with_angle_offsets!: angle_reference = \"$mode\" is not " *
        "recognised; use \"global\" (default) or \"local\".")

    fund_bus = hdata["nw"]["1"]["bus"]
    for (n, ntw) in hdata["nw"]
        n == "1" && continue
        nh = parse(Int, n)
        for load in values(ntw["load"])
            θ = mode == "local" ? fund_bus["$(load["load_bus"])"]["va"] : 0.0
            load["angle_offset"] = nh * θ
        end
    end
end