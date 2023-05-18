################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# solve functions
""
function solve_hopf(file, model_type::Type, optimizer; kwargs...)
    return _PMs.solve_model(file, model_type, optimizer, build_hopf; ref_extensions=[ref_add_xfmr!],  solution_processors=[ _HPM.sol_data_model!], multinetwork=true, kwargs...)
end

""
function solve_mc_hopf(file, model_type::Type, optimizer; kwargs...)
    return _PMD.solve_mc_model(file, model_type, optimizer, build_mc_hopf; multinetwork=true, kwargs...)
end

# build functions
""
function build_hopf(pm::_PMs.AbstractIVRModel)
    ## variables
    for n in _PMs.nw_ids(pm)
        _PMs.variable_bus_voltage(pm, nw=n, bounded=false)
        _PMs.variable_bus_voltage_magnitude_sqr(pm, nw=n, bounded=false)
        variable_transformer_voltage(pm, nw=n, bounded=false)
        
        _PMs.variable_branch_current(pm, nw=n, bounded=false)
        variable_transformer_current(pm, nw=n, bounded=false)

        variable_load_current(pm, nw=n, bounded=false)
        variable_gen_current(pm, nw=n, bounded=false) 
    end 

    ## objective
    objective_voltage_distortion_minimization(pm)

    ## constraint
    # overall constraints
    for i in _PMs.ids(pm, :bus, nw=1)
        constraint_bus_voltage_rms_limit(pm, i, nw=1)
        constraint_bus_voltage_thd_limit(pm, i, nw=1)
    end

    for b in _PMs.ids(pm, :branch, nw=1)
        constraint_branch_current_rms_limit(pm, b, nw=1)
    end

    for g in _PMs.ids(pm, :gen, nw=1)                                           # @F: would it not be cleaner to make a new data entry for filters
        constraint_active_filter(pm, g, nw=1)
    end

    # harmonic constraints
    for n in _PMs.nw_ids(pm)
        for i in _PMs.ids(pm, :ref_buses, nw=n) 
            constraint_voltage_reference(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_bus_voltage_ihd_limit(pm, i, nw=n)
            constraint_bus_voltage_magnitude_sqr(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)
            
            _PMs.constraint_voltage_drop(pm, b, nw=n)

            _PMs.constraint_current_limit(pm, b, nw=n)
        end
    
        for g in _PMs.ids(pm, :gen, nw=n)
            _PMs.constraint_gen_active_bounds(pm, g, nw=n)
            _PMs.constraint_gen_reactive_bounds(pm, g, nw=n)
        end

        for l in _PMs.ids(pm, :load, nw=n)
            constraint_load_power(pm, l, nw=n)
        end
        
        for t in _PMs.ids(pm, :xfmr, nw=n)
            constraint_transformer_core_excitation(pm, t, nw=n)
            constraint_transformer_core_voltage_drop(pm, t, nw=n)
            constraint_transformer_core_voltage_balance(pm, t, nw=n)
            constraint_transformer_core_current_balance(pm, t, nw=n)
            
            constraint_transformer_winding_config(pm, t, nw=n)
            constraint_transformer_winding_current_balance(pm, t, nw=n)
        end
    end
end

""
function build_mc_hopf(pm::_PMD.AbstractExplicitNeutralIVRModel)
    ## variables
    for n in _PMD.nw_ids(pm)
        _PMD.variable_mc_bus_voltage(pm; nw=n)
        _PMD.variable_mc_bus_voltage_magnitude_sqr(pm; nw=n)

        _PMD.variable_mc_branch_current(pm; nw=n)
        _PMD.variable_mc_generator_current(pm; nw=n)
        _PMD.variable_mc_load_current(pm; nw=m)
        _PMD.variable_mc_transformer_current(pm; nw=n)

        _PMD.variable_mc_load_power(pm; nw=n)
        _PMD.variable_mc_generator_power(pm; nw=n)
        _PMD.variable_mc_transformer_power(pm; nw=n)
    end

    ## objective
    objective_mc_voltage_distortion_minimization(pm)

    ## constraints
    for i in _PMD.ids(pm, :bus, nw=1)
        constraint_mc_bus_voltage_rms_limit(pm, i, nw=1)
        constraint_mc_bus_voltage_thd_limit(pm, i, nw=1)
    end

    for b in _PMD.ids(pm, :branch, nw=1)
        constraint_mc_branch_current_rms_limit(pm, b, nw=1)
    end

    for g in _PMD.ids(pm, :gen, nw=1)
        constraint_mc_active_filter(pm, g, nw=1)
    end

    # harmonic constraints
    for n in _PMD.nw_ids(pm)
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            _PMD.constraint_mc_voltage_reference(pm, i, nw=n)
        end

        for i in _PMD.ids(pm, :bus, nw=n)
            constraint_mc_bus_voltage_ihd_limit(pm, i, nw=n)
            constraint_mc_bus_voltage_magnitude_sqr(pm, i, nw=n)
        end

        for b in _PMD.ids(pm, n, :branch)
            _PMD.constraint_mc_current_from(pm, b, nw=n)
            _PMD.constraint_mc_current_to(pm, b, nw=n)

            _PMD.constraint_mc_bus_voltage_drop(pm, b, nw=n)
        end

        for g in _PMD.ids(pm, n, :gen)
            _PMD.constraint_mc_generator_power(pm, g, nw=n)
            _PMD.constraint_mc_generator_current(pm, g, nw=n)
        end

        for l in _PMD.ids(pm, n, :load)
            constraint_mc_load_power(pm, l, nw=n)
            constraint_mc_load_current(pm, l, nw=n)
        end

        for t in _PMD.ids(pm, n, :transformer)
            _PMD.constraint_mc_transformer_voltage(pm, t, nw=n)
            _PMD.constraint_mc_transformer_current(pm, t, nw=n)
        end

        for i in _PMD.ids(pm, n, :bus)
            _PMD.constraint_mc_current_balance(pm, i, nw=n)
        end
    end
end