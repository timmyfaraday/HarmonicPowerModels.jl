################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
################################################################################

""
function solve_hpf(hdata, model_type::Type, optimizer; kwargs...)
    return _PMs.solve_model(hdata, model_type, optimizer, build_hpf; 
                                ref_extensions=[ref_add_filter!,
                                                ref_add_xfmr!], 
                                solution_processors=[_HPM.sol_data_model!], 
                                multinetwork=true, kwargs...)
end

""
function build_hpf(pm::_PMs.AbstractIVRModel)
    # variables
    for n in _PMs.nw_ids(pm)
        ## voltage variables
        variable_bus_voltage(pm, nw=n, bounded=false)
        variable_transformer_voltage(pm, nw=n, bounded=false)
        
        ## edge current variables
        variable_branch_current(pm, nw=n, bounded=false)
        variable_transformer_current(pm, nw=n, bounded=false)

        ## unit current variables
        variable_filter_current(pm, nw=n, bounded=false)
        variable_gen_current(pm, nw=n, bounded=false)
        variable_load_current(pm, nw=n, bounded=false)        
    end 

    # objective
    objective_power_flow(pm)

    # constraint
    ## overall or fundamental constraints
    ### filter
    for f in _PMs.ids(pm, :filter, nw=fundamental(pm))
        constraint_active_filter(pm, f, nw=fundamental(pm))
    end

    ## harmonic constraints
    for n in _PMs.nw_ids(pm)
        ### reference node
        for i in _PMs.ids(pm, :ref_buses, nw=n) 
            constraint_ref_bus(pm, i, nw=n)
        end

        ### node
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end

        ### branch 
        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)
            
            _PMs.constraint_voltage_drop(pm, b, nw=n)
        end

        ### xfmr 
        for t in _PMs.ids(pm, :xfmr, nw=n)
            constraint_transformer_core_magnetization(pm, t, nw=n)
            constraint_transformer_core_voltage_drop(pm, t, nw=n)
            constraint_transformer_core_voltage_balance(pm, t, nw=n)
            constraint_transformer_core_current_balance(pm, t, nw=n)
            
            constraint_transformer_winding_config(pm, t, nw=n)
            constraint_transformer_winding_current_balance(pm, t, nw=n)
        end

        ### harmonic unit
        for l in _PMs.ids(pm, :load, nw=n)
            constraint_load_power(pm, l, nw=n)
        end
    end
end
