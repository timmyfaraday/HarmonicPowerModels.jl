""
function run_hopf_iv(file, model_type::Type, optimizer; kwargs...)
    return _PMs.run_model(file, model_type, optimizer, build_hopf_iv; ref_extensions=[ref_add_xfmr!], multinetwork=true, kwargs...)
end

""
function build_hopf_iv(pm::AbstractPowerModel)
    bounded= true

    for (n, network) in _PMs.nws(pm)
        _PMs.variable_bus_voltage(pm, nw=n, bounded=bounded)
        _PMs.variable_bus_voltage_magnitude_sqr(pm, nw=n, bounded=bounded)
        variable_transformer_voltage(pm, nw=n, bounded=bounded)
        
        _PMs.variable_branch_current(pm, nw=n)
        _PMs.variable_dcline_current(pm, nw=n)
        variable_transformer_current(pm, nw=n, bounded=bounded)

        variable_load_current(pm, nw=n, bounded=bounded)
        variable_gen_current(pm, nw=n, bounded=bounded)
        
    end 

    for (n, network) in _PMs.nws(pm)
        for i in _PMs.ids(pm, :ref_buses, nw=n) 
            constraint_ref_bus(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_vm_auxiliary_variable(pm, i, nw=n)
        end
    
        for g in _PMs.ids(pm, :gen, nw=n)
            _PMs.constraint_gen_active_bounds(pm, g, nw=n)
            _PMs.constraint_gen_reactive_bounds(pm, g, nw=n)
        end

        for i in _PMs.ids(pm, :load, nw=n)
            constraint_load_power(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)
            
            _PMs.constraint_voltage_drop(pm, b, nw=n)

            #TODO add current magnitude constraints instead
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

        # for d in _PMs.ids(pm, :dcline, nw=n)
        #     _PMs.constraint_dcline_power_losses(pm, d, nw=n)
        # end
    end
    
    #constraints across harmonics
    fundamental = 1
    for i in _PMs.ids(pm, :bus, nw=fundamental)
        # constraint_voltage_magnitude_rms(pm, i)
        # constraint_voltage_thd(pm, i, fundamental=fundamental)
    end

    for b in _PMs.ids(pm, :branch, nw=fundamental)
        # constraint_current_limit_rms(pm, b)
    end

    _PMs.objective_min_fuel_and_flow_cost(pm)
end

function ref_add_xfmr!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any}) ## data not actually needed!
    _PMs.apply_pm!(_ref_add_xfmr!, ref, data)
end

function _ref_add_xfmr!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if !haskey(ref, :xfmr)
        # error(_LOGGER, "required xfmr data not found")
        ref[:xfmr] = Dict()
        ref[:xfmr_arcs_from] = Dict()
        ref[:xfmr_arcs_to] = Dict()
        ref[:xfmr_arcs] = Dict()
        ref[:bus_arcs_xfmr] = Dict((i, []) for (i,bus) in ref[:bus])

    else
        ref[:xfmr] = Dict(x for x in ref[:xfmr] if  x.second["f_bus"] in keys(ref[:bus]) &&
                                                    x.second["t_bus"] in keys(ref[:bus])
            )
        
        ref[:xfmr_arcs_from] = [(t,xfmr["f_bus"],xfmr["t_bus"]) for (t,xfmr) in ref[:xfmr]]
        ref[:xfmr_arcs_to]   = [(t,xfmr["t_bus"],xfmr["f_bus"]) for (t,xfmr) in ref[:xfmr]]

        ref[:xfmr_arcs] = [ref[:xfmr_arcs_from]; ref[:xfmr_arcs_to]]

        bus_arcs_xfmr = Dict((i, []) for (i,bus) in ref[:bus])
        for (t,i,j) in ref[:xfmr_arcs]
            push!(bus_arcs_xfmr[i], (t,i,j))
        end
        ref[:bus_arcs_xfmr] = bus_arcs_xfmr
    end
end