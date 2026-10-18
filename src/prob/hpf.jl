################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
# v0.3.0 - reviewed TVA                                                        #
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
        variable_xfmr_voltage(pm, nw=n, bounded=false)
        
        ## edge current variables
        variable_branch_current(pm, nw=n, bounded=false)
        variable_xfmr_current(pm, nw=n, bounded=false)

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
    for f in ids(pm, :filter)
        constraint_active_filter(pm, f, nw=fundamental(pm))
    end
    ### generator — fundamental power bounds (P/Q based on the bus voltage)
    for g in ids(pm, :gen)
        _PMs.constraint_gen_active_bounds(pm, g, nw=fundamental(pm))
        _PMs.constraint_gen_reactive_bounds(pm, g, nw=fundamental(pm))
    end

    ## harmonic constraints
    for n in _PMs.nw_ids(pm)
        ### reference / clean buses — the fundamental slack fixes a voltage; at
        ### the harmonics only clean buses are forced to V=0 (ideal ground), all
        ### other buses are grounded through the generator subtransient impedance
        ### and the line/bus shunts.
        if n == fundamental(pm)
            for i in _PMs.ids(pm, :ref_buses, nw=n)
                constraint_voltage_ref_bus(pm, i, nw=n)
            end
        else
            for i in _PMs.ids(pm, :bus, nw=n)
                constraint_voltage_clean_bus(pm, i, nw=n)
            end
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

        ### generator — harmonic subtransient-impedance shunt to ground
        for g in _PMs.ids(pm, :gen, nw=n)
            constraint_gen_current(pm, g, nw=n)
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

        ### harmonic unit
        for l in _PMs.ids(pm, :load, nw=n)
            constraint_load_power(pm, l, nw=n)
        end
    end
end
