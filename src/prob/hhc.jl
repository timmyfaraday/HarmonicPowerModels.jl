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
function solve_hhc(hdata, model_type::Type, optimizer; kwargs...)
    # update hdata for chosen fairness principle
    update_hdata_with_fairness_principle_data!(hdata, model_type, optimizer) 

    # solve non-linear harmonic hosting capacity problem
    return _PMs.solve_model(hdata, model_type, optimizer, build_hhc; 
                                ref_extensions=[ref_add_filter!,
                                                ref_add_xfmr!], 
                                solution_processors=[ _HPM.sol_data_model!], 
                                multinetwork=true, kwargs...)
end
""
function solve_hhc(hdata, model_type::Type, hhc_optimizer, hpf_optimizer; kwargs...)
    # solve fundamental harmonic power flow problem and update hdata
    #update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, hpf_optimizer)

    # update hdata for chosen fairness principle
    update_hdata_with_fairness_principle_data!(hdata, dHHC_SOC, hhc_optimizer) 

    # solve second order cone harmonic hosting capacity problem
    return _PMs.solve_model(hdata, model_type, hhc_optimizer, build_hhc; 
                                ref_extensions=[ref_add_filter!,
                                                ref_add_xfmr!], 
                                solution_processors=[ _HPM.sol_data_model!], 
                                multinetwork=true, kwargs...)
end

""
function build_hhc(pm::HarmonicPowerModel)
    # variables 
    for n in _PMs.nw_ids(pm)
        ## fairness variable
        if n ≠ fundamental(pm)
            variable_fairness_principle(pm, nw=n, bounded=true)
        end 

        ## voltage variables 
        variable_bus_voltage(pm, nw=n, bounded=true)
        variable_xfmr_voltage(pm, nw=n, bounded=true)

        ## edge current variables
        variable_branch_current(pm, nw=n, bounded=true)
        variable_xfmr_current(pm, nw=n, bounded=true)

        ## unit current variables
        variable_filter_current(pm, nw=n, bounded=false)
        variable_gen_current(pm, nw=n, bounded=true)
        variable_hsrc_current(pm, nw=n, bounded=true)
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
    ### generator
    for g in ids(pm, :gen)
        constraint_gen_current_rms_limit(pm, g)

        _PMs.constraint_gen_active_bounds(pm, g, nw=fundamental(pm))   # CHECK IF NEEDED
        _PMs.constraint_gen_reactive_bounds(pm, g, nw=fundamental(pm))    # CHECK IF NEEDED
    end
    ### xfmr 
    for x in ids(pm, :xfmr)
        constraint_xfmr_current_rms_limit(pm, x)
    end

    ## harmonic constraints
    for n in _PMs.nw_ids(pm)
        ### fairness principle
        if n ≠ fundamental(pm)
            constraint_fairness_principle(pm, nw=n)
        end

        ### reference bus
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_ref_voltage(pm, i, nw=n)
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

        ### generator
        for g in _PMs.ids(pm, :gen, nw=n)
            constraint_gen_current(pm, g, nw=n)
        end

        ### harmonic source
        for s in _PMs.ids(pm, :source, nw=n)
            constraint_hsrc_current(pm, s, nw = n)
        end

        ### xfmr
        for x in _PMs.ids(pm, :xfmr, nw=n)
            constraint_xfmr_core_magnetization(pm, x, nw=n)
            constraint_xfmr_core_voltage_drop(pm, x, nw=n)
            constraint_xfmr_core_voltage_phase_shift(pm, x, nw=n)
            constraint_xfmr_core_current_balance(pm, x, nw=n)
            
            constraint_xfmr_winding_config(pm, x, nw=n)
            constraint_xfmr_winding_current_balance(pm, x, nw=n)
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
        if n ≠ fundamental(pm)
            variable_fairness_principle(pm, nw=n, bounded=true)
        end 

        ## voltage variables 
        variable_bus_voltage(pm, nw=n, bounded = true)
        variable_xfmr_voltage(pm, nw=n, bounded = true)

        ## edge current variables
        variable_branch_current(pm, nw=n, bounded = true)
        variable_xfmr_current(pm, nw=n, bounded = true)

        ## node current variables
        variable_filter_current(pm, nw=n, bounded = true)
        variable_hsrc_current(pm, nw=n, bounded = true)
        variable_gen_current(pm, nw=n, bounded = true)
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
        constraint_current_rms_limit(pm, b)
    end
    ### xfmr 
    for x in ids(pm, :xfmr)
        constraint_xfmr_current_rms_limit(pm, x)
    end
    ### generator 
    for g in ids(pm, :gen)
        constraint_gen_current_rms_limit(pm, g)
    end

    ## harmonic constraints
    for n in _PMs.nw_ids(pm) if n ≠ fundamental(pm)
        ### fairness principle
        constraint_fairness_principle(pm, nw=n)
        
        ### reference node
        if !haskey(pm.setting, "fix_refbus_angle") || pm.setting["fix_refbus_angle"] == true 
            for i in _PMs.ids(pm, :ref_buses, nw=n)
                constraint_ref_voltage(pm, i, nw=n)
            end
        end

        ### node
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_bus_current_balance(pm, i, nw=n)
            constraint_bus_voltage_ihd_limit(pm, i, nw=n)
        end

        ### branch
        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)

            _PMs.constraint_voltage_drop(pm, b, nw=n)
        end

        ### xfmr
        for x in _PMs.ids(pm, :xfmr, nw=n)
            constraint_xfmr_core_magnetization(pm, x, nw=n)
            constraint_xfmr_core_voltage_drop(pm, x, nw=n)
            constraint_xfmr_core_voltage_phase_shift(pm, x, nw=n)
            constraint_xfmr_core_current_balance(pm, x, nw=n)
            
            constraint_xfmr_winding_config(pm, x, nw=n)
            constraint_xfmr_winding_current_balance(pm, x, nw=n)
        end

        ### harmonic source
        for s in _PMs.ids(pm, :source, nw=n)
            constraint_hsrc_current(pm, s, nw = n)
        end
    end end
end