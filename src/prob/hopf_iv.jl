""
function run_hopf_iv(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_hopf_iv; multinetwork=true, kwargs...)
end

""
function build_hopf_iv(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        _PMs.variable_bus_voltage(pm, nw=n)
        variable_transformer_voltage(pm, nw=n)
        
        _PMs.variable_gen_current(pm, nw=n)
        _PMs.variable_branch_current(pm, nw=n)
        _PMs.variable_dcline_current(pm, nw=n)
        variable_transformer_current(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            _PMs.constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end

        for b in ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)
            
            _PMs.constraint_voltage_drop(pm, b, nw=n)
            _PMs.constraint_voltage_angle_difference(pm, b, nw=n)

            _PMs.constraint_thermal_limit_from(pm, b, nw=n)
            _PMs.constraint_thermal_limit_to(pm, b, nw=n)
        end
        
        for t in ids(pm, :xfmr, nw=n)
            constraint_transformer_core_voltage_drop(pm, t, nw=n)
            constraint_transformer_core_voltage_balance(pm, t, nw=n)
            constraint_transformer_core_current_balance(pm, t, nw=n)
            
            constraint_transformer_winding_config(pm, t, nw=n)
            constraint_transformer_winding_current_balance(pm, t, nw=n)
        end

        for d in ids(pm, :dcline, nw=n)
            _PMs.constraint_dcline_power_losses(pm, d, nw=n)
        end
    end

    for b in ids(pm, :bus)
        constraint_voltage_magnitude_rms(pm, i, nw=n)
    end

    _PMs.objective_min_fuel_and_flow_cost(pm)
end