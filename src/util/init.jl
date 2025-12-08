################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Hakan Ergun                                                         #
################################################################################
# Changelog:                                                                   #
#                                                                              #
################################################################################
function update_hdata_with_fundamental_hpf_results!(hdata, model_type::Type, optimizer; init = "pf")
    # remove all but the fundamental network
    hpf_data = deepcopy(hdata)
    for n in keys(hpf_data["nw"])
        if n ≠ "1"
            delete!(hpf_data["nw"], n)
        end
    end

    if init .== "pf"
        # solve hpf problem for the fundamental harmonic
        fpf = solve_hpf(hpf_data, HarmonicPowerModel, optimizer)["solution"]["nw"]["1"]
    elseif init .== "opf"
        fpf = solve_hopf(hpf_data, HarmonicPowerModel, optimizer)["solution"]["nw"]["1"]
    end

    # update hdata with the results of the hpf problem
    for (i, bus) in hdata["nw"]["1"]["bus"]
        bus["v_fund_magn"] = fpf["bus"][i]["vbm"]
    #    bus["v_rms_max"] = max(fpf["bus"][i]["vbm"], bus["v_rms_max"])
    #    bus["v_rms_min"] = min(fpf["bus"][i]["vbm"], bus["v_rms_min"])
    end
    for (b, branch) in hdata["nw"]["1"]["branch"]
        branch["i_fund_magn"] = [fpf["branch"][b]["cbm_fr"], fpf["branch"][b]["cbm_to"]]
        #branch["i_rms_max"] = max(maximum(branch["i_fund_magn"]),  branch["i_rms_max"])
    end
    for (x, xfmr) in hdata["nw"]["1"]["xfmr"]
        xfmr["i_fund_magn"] = [fpf["xfmr"][x]["cxm_$nb"] for nb in xfmr["bus"]]
        # xfmr["i_rms_max"] = [max(xfmr["i_rms_max"], maximum(xfmr["i_fund_magn"])), max(xfmr["i_rms_max"], maximum(xfmr["i_fund_magn"]))]
    end
    for (g, gen) in hdata["nw"]["1"]["gen"]
        gen["i_fund_magn"] = fpf["gen"][g]["cgm"]
        # gen["i_rms_max"] = max(gen["i_rms_max"], gen["i_fund_magn"])
    end
end

function calculate_maximum_harmonic_source_current_injection!(data, Zh)    
    for (nw, network) in data["nw"]
      if network !== "nw"
          h = parse(Int, nw)
          for (r, hsrc) in network["hsrc"]
              bus = data["hsrc"]["nw"]["1"]["bus"]
              v_ihd_max = network["bus"]["$bus"]["v_ihd_max"] 
              hb_idx = network["bus"]["$bus"]["index"]
              load["cmdmax"] = v_ihd_max / sqrt(Zh["$hb_idx"][h]["re"]^2 + Zh["$hb_idx"][h]["im"]^2)
          end
      end
    end
end



################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
# v0.3.0 - adapted for extended graph representation                           #
################################################################################

""
function solve_hhc_init(hdata, model_type::Type, optimizer; optimizer_soc = optimizer, kwargs...)
    
    return solve_model(hdata, model_type, optimizer, build_hhc_init; 
                            solution_processors=[_HPM.sol_data_model!], 
                            multinetwork=true, 
                            kwargs...)
end

""
function build_hhc_init(pm::HarmonicPowerModel)
    # variables 
    for n in _PMs.nw_ids(pm)
        ## fairness variable
        if n ≠ fundamental(pm)
            variable_fairness_principle(pm, nw = n, bounded=true)
        end
        add_hhc_variables(pm, n)
    end

    # objective 
    #objective_maximum_hosting_capacity(pm)

    # constraints 
    ## overall or fundamental constraints
    ### node 
    for i in ids(pm, :bus)
        constraint_bus_voltage_rms_limit(pm, i)
        #constraint_bus_voltage_thd_limit(pm, i)
    end

    ### branch
    for b in ids(pm, :branch)
        constraint_branch_current_rms_limit(pm, b)
    end
    
    ## xfmr 
    for x in ids(pm, :xfmr)
        constraint_xfmr_winding_current_rms_limit(pm, x)
    end
    
    ### filter
    for f in ids(pm, :filter)
        constraint_filter_current(pm, f)
        constraint_filter_current_rms_limit(pm, f)
    end
    ### generator
    for g in ids(pm, :gen)
        constraint_gen_current_rms_limit(pm, g)
        constraint_gen_power_active_fundamental_limit(pm, g)
        constraint_gen_power_reactive_fundamental_limit(pm, g)
    end

    ## harmonic constraints
    for n in _PMs.nw_ids(pm)
        ### fairness principle
        # if n ≠ fundamental(pm)
        #     constraint_fairness_principle(pm, nw=n)
        # end

        ### reference bus
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_ref_voltage(pm, i, nw=n)
        end
        
        # ### clean bus 
        for i in _PMs.ids(pm, :clean_buses, nw=n)
            constraint_clean_voltage(pm, i, nw=n)
        end

        # ### bus
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_bus_current_balance(pm, i, nw=n)
            constraint_bus_voltage_ihd_limit(pm, i, nw=n)
        end

        # ### branch
        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_branch_current_from(pm, b, nw=n)
            constraint_branch_current_to(pm, b, nw=n)
            constraint_branch_voltage_drop(pm, b, nw=n)
        end

        # ### xfmr
        for x in _PMs.ids(pm, :xfmr, nw=n)
            constraint_xfmr_core_magnetization(pm, x, nw=n)
            constraint_xfmr_core_voltage_drop(pm, x, nw=n)
            constraint_xfmr_core_voltage_phase_shift(pm, x, nw=n)
            constraint_xfmr_core_current_balance(pm, x, nw=n)
            
            constraint_xfmr_winding_current_balance(pm, x, nw=n)
            constraint_xfmr_winding_voltage_drop(pm, x, nw=n)
            constraint_xfmr_winding_zero_seq_current_blocking(pm, x, nw=n) # check this constraint!
        end

        ### generator
        for g in _PMs.ids(pm, :gen, nw=n)
            constraint_gen_current(pm, g, nw=n)
        end

        ## harmonic source
        for r in _PMs.ids(pm, :hsrc, nw=n)
            constraint_hsrc_current(pm, r, nw=n)
        end

        ## harmonic load
        for l in _PMs.ids(pm, :hload, nw=n)
            constraint_hload_current(pm, l, nw=n)
        end

        ### shunt
        for s in _PMs.ids(pm, :shunt, nw=n)
            constraint_shunt_current(pm, s, nw=n)
        end
    end
end
