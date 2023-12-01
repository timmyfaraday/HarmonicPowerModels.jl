################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

""
function solve_hhc(file, model_type::Type, optimizer; kwargs...)
    return _PMs.solve_model(file, model_type, optimizer, build_hhc; ref_extensions=[ref_add_xfmr!],  solution_processors=[ _HPM.sol_data_model!], multinetwork=true, kwargs...)
end

""
function build_hhc(pm::_PMs.AbstractIVRModel)
    ## variables 
    for n in _PMs.nw_ids(pm)
        # voltage variables 
        _PMs.variable_bus_voltage(pm, nw=n, bounded=false)
        _PMs.variable_bus_voltage_magnitude_sqr(pm, nw=n, bounded=false)
        variable_transformer_voltage(pm, nw=n, bounded=false)

        # edge current variables
        _PMs.variable_branch_current(pm, nw=n, bounded=false)
        variable_transformer_current(pm, nw=n, bounded=false)

        # node current variables
        variable_load_current(pm, nw=n, bounded=true)
        variable_gen_current(pm, nw=n, bounded=false)
    end

    ## objective 
    objective_maximum_hosting_capacity(pm)

    ## constraints 
    # overall constraints
    for i in _PMs.ids(pm, :bus, nw=1)
        constraint_voltage_rms_limit(pm, i, nw=1)
        constraint_voltage_thd_limit(pm, i, nw=1)
    end

    for b in _PMs.ids(pm, :branch, nw=1)
        constraint_current_rms_limit(pm, b, nw=1)
    end

    # harmonic constraints
    for n in _PMs.nw_ids(pm)
        # node
        for i in _PMs.ids(pm, :ref_buses, nw=1)
            constraint_ref_bus(pm, i, nw=n)
        end
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_voltage_ihd_limit(pm, i, nw=n)
            constraint_voltage_magnitude_sqr(pm, i, nw=n)
        end

        # edge
        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)
            
            _PMs.constraint_voltage_drop(pm, b, nw=n)

            _PMs.constraint_current_limit(pm, b, nw=n)
        end
        for t in _PMs.ids(pm, :xfmr, nw=n)
            constraint_transformer_core_excitation(pm, t, nw=n)
            constraint_transformer_core_voltage_drop(pm, t, nw=n)
            constraint_transformer_core_voltage_balance(pm, t, nw=n)
            constraint_transformer_core_current_balance(pm, t, nw=n)
            
            constraint_transformer_winding_config(pm, t, nw=n)
            constraint_transformer_winding_current_balance(pm, t, nw=n)
        end

        # unit
        for g in _PMs.ids(pm, :gen, nw=n)
            _PMs.constraint_gen_active_bounds(pm, g, nw=n)
            _PMs.constraint_gen_reactive_bounds(pm, g, nw=n)
        end
        for l in _PMs.ids(pm, :load, nw=n)
            constraint_load_current(pm, l, nw = n)
        end
    end
end