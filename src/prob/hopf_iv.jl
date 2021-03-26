""
function run_harmonic_opf_iv(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_harmonic_opf_iv; multinetwork=true, kwargs...)
end

""
function build_harmonic_opf_iv(pm::AbstractPowerModel)
    for n in nws(pm)
        _PMs.variable_bus_voltage(pm, nw=n)
        variable_transformer_voltage(pm, nw=n)
        
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
        
        for t in ids(pm, :transformer, nw=n)
            constraint_current_transformer(pm, t, nw=n)
            constraint_voltage_drop_transformer(pm, t, nw=n)
        end

        for d in ids(pm, :dcline, nw=n)
            _PMs.constraint_dcline_power_losses(pm, d, nw=n)
        end
    end

    for i in ids(pm, :bus)
        constraint_voltage_magnitude_rms(pm, i)
    end

    for t in ids(pm, :transformer)
        constraint_voltage_transformer(pm, t)
    end

    objective_min_fuel_and_flow_cost(pm)
end