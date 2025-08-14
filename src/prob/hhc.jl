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
function solve_hhc(hdata, model_type::Type, optimizer; optimizer_soc = optimizer, kwargs...)
    if model_type == dHHCPowerModel
        update_hdata_with_fundamental_hpf_results!(hdata, HarmonicPowerModel, optimizer)
    end
    
    return solve_model(hdata, model_type, optimizer_soc, build_hhc; 
                            solution_processors=[_HPM.sol_data_model!], 
                            multinetwork=true, 
                            kwargs...)
end

""
function add_hhc_variables(pm, n)
    ## voltage variables 
    variable_bus_voltage(pm, nw=n, bounded=true)
    variable_xfmr_voltage(pm, nw=n, bounded=true)

    ## edge current variables
    variable_branch_current(pm, nw=n, bounded=true)
    variable_xfmr_current(pm, nw=n, bounded=true)

    ## unit current variables
    variable_filter_current(pm, nw=n, bounded=true)
    variable_gen_current(pm, nw=n, bounded=true)
    variable_hload_current(pm, nw=n, bounded=true)                         # empty variable
    variable_hsrc_current(pm, nw=n, bounded=true)
    variable_shunt_current(pm, nw=n, bounded=true)

    # edge power variables
    variable_branch_power(pm, nw=n, bounded=false)
    variable_xfmr_power(pm, nw=n, bounded=false)

    # unit power variables
    variable_gen_power(pm, nw=n, bounded=false)
end

""
function build_hhc(pm::HarmonicPowerModel)
    # variables 
    for n in _PMs.nw_ids(pm)
        ## fairness variable
        if n ≠ fundamental(pm)
            variable_fairness_principle(pm, nw=n, bounded=true)
        end
        add_hhc_variables(pm, n) 
    end

    # objective 
    objective_maximum_hosting_capacity(pm)

    # constraints 
    ## overall or fundamental constraints
    ### node 
    for i in ids(pm, :bus)
        constraint_bus_voltage_rms_limit(pm, i)
        constraint_bus_voltage_thd_limit(pm, i)
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
        if n ≠ fundamental(pm)
            constraint_fairness_principle(pm, nw=n)
        end

        # ### reference bus
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

        ### harmonic load
        for l in _PMs.ids(pm, :hload, nw=n)
            constraint_hload_power(pm, l, nw=n)
        end

        ### shunt
        for s in _PMs.ids(pm, :shunt, nw=n)
            constraint_shunt_current(pm, s, nw=n)
        end
    end
end

""
function build_hhc(pm::dHHCPowerModel)
    # add SOCtoNonConvexQuadBridge
    JuMP.add_bridge(pm.model, _MOI.Bridges.Constraint.SOCtoNonConvexQuadBridge)

    # variables 
    for n in _PMs.nw_ids(pm) if n ≠ fundamental(pm)
        ## fairness variable 
        variable_fairness_principle(pm, nw=n, bounded=true)
        add_hhc_variables(pm, n)    
    end end

    # objective 
    objective_maximum_hosting_capacity(pm)

    # constraints 
    ## overall constraints
    ### node
    for i in ids(pm, :bus)
        constraint_bus_voltage_rms_limit(pm, i)
        constraint_bus_voltage_thd_limit(pm, i)
    end

    ### branch
    for b in ids(pm, :branch)
        constraint_branch_current_rms_limit(pm, b)
    end
    ### xfmr 
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
    end

    ## harmonic constraints
    for n in _PMs.nw_ids(pm) if n ≠ fundamental(pm)
        ### fairness principle
        constraint_fairness_principle(pm, nw=n)
        
        ### reference bus
        for i in _PMs.ids(pm, :ref_buses, nw=n)
                constraint_ref_voltage(pm, i, nw=n)
        end
        
        ### clean bus 
        for i in _PMs.ids(pm, :clean_buses, nw=n)
            constraint_clean_voltage(pm, i, nw=n)
        end

        ### bus
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_bus_current_balance(pm, i, nw=n)
            constraint_bus_voltage_ihd_limit(pm, i, nw=n)
        end

        ### branch
        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_branch_current_from(pm, b, nw=n)
            constraint_branch_current_to(pm, b, nw=n)
            constraint_branch_voltage_drop(pm, b, nw=n)
        end

        ### xfmr
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

        ### harmonic source
        for r in _PMs.ids(pm, :hsrc, nw=n)
            constraint_hsrc_current(pm, r, nw = n)
        end

        ### harmonic load
        for l in _PMs.ids(pm, :hload, nw=n)
            constraint_hload_power(pm, l, nw=n)
        end

        ### shunt 
        for s in _PMs.ids(pm, :shunt, nw=n)
            constraint_shunt_current(pm, s, nw=n)
        end
end end end

