################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# variables ####################################################################
## fairness principle variables ################################################
""
function variable_fairness_principle(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    # maximum efficiency
    if pm.data["hhc_principle"] == "maximum efficiency"
        # no additional variables
    end

    # absolute equality
    if pm.data["hhc_principle"] == "absolute equality"
        # no additional variables 
    end

    # maximin
    if pm.data["hhc_principle"] == "maximin"
        cmh = _PMs.var(pm, nw)[:cmh] = 
                JuMP.@variable( pm.model, 
                                base_name="$(nw)_cmh",
                                start=0.0)

        if bounded
            JuMP.set_lower_bound(cmh, 0.0)
        end

        report && (_PMs.sol(pm, nw, :fairness)[:cmh] = cmh)
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["hhc_principle"] == "Kalai-Smorodinsky bargaining"
        fh  = _PMs.var(pm, nw)[:fh] = 
                JuMP.@variable( pm.model, 
                                base_name="$(nw)_fh",
                                start=0.0)

        if bounded
            JuMP.set_lower_bound(fh, 0.0)
            JuMP.set_upper_bound(fh, 1.0)
        end

        report && (_PMs.sol(pm, nw, :fairness)[:fh] = fh)
end end

# constraints ##################################################################
## fairness principle constraints ##############################################
""
function constraint_fairness_principle(pm::AbstractHarmonicModel; nw::Int=fundamental(pm))
    ids = sort(collect(_PMs.ids(pm, :hsrc, nw=nw)))

    constraint_fairness_principle(pm, nw, ids)
end
""
function constraint_fairness_principle(pm::AbstractHarmonicModel, n, ids)
    # maximum efficiency
    if pm.data["hhc_principle"] == "maximum efficiency"
        # no additional constraints
    end

    # absolute equality
    if pm.data["hhc_principle"] == "absolute equality"
        crm = [_PMs.var(pm, n, :crm, r) for r in ids]

        # for r in ids[2:end]
        #     JuMP.@constraint(pm.model, crm[first(ids)] == crm[r])
        # end 
    end

    # maximin
    if pm.data["hhc_principle"] == "maximin"
        cmh = _PMs.var(pm, n, :cmh)
        for r in ids
            crm = _PMs.var(pm, n, :crm, r)
            JuMP.@constraint(pm.model, cmh <= crm)
        end 
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["hhc_principle"] == "Kalai-Smorodinsky bargaining"
        fh = _PMs.var(pm, n, :fh)

        for r in ids
            crm = _PMs.var(pm, n, :crm, r)
            crmax = _PMs.ref(pm, n, :hsrc, r, "crmax")

            JuMP.@constraint(pm.model, crm == fh * crmax)
        end 
    end
end

# objective ####################################################################
## harmonic power flow objective ###############################################
""
function objective_power_flow(pm::_PMs.AbstractIVRModel)
    JuMP.@objective(pm.model, Min, 0.0)
end

## voltage distortion minimization objective ###################################
""
function objective_voltage_distortion_minimization(pm::_PMs.AbstractIVRModel) 
    bus_id = pm.data["bus_id"]

    vbr = [_PMs.var(pm, n, :vbr, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]
    vbi = [_PMs.var(pm, n, :vbi, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]

    JuMP.@objective(pm.model, Min, sum(vbr.^2 + vbi.^2))
end

## harmonic hosting capacity objective #########################################
""
function objective_maximum_hosting_capacity(pm::_PMs.AbstractIVRModel)
    # maximum efficiency
    if pm.data["hhc_principle"] == "maximum efficiency"
        crm = [_PMs.var(pm, n, :crm, r) for n in _PMs.nw_ids(pm) 
                                        for r in _PMs.ids(pm, :hsrc, nw=n) 
                                        if n ≠ fundamental(pm)]
    
        JuMP.@objective(pm.model, Max, sum(crm)) 
    end

    # absolute equality
    if pm.data["hhc_principle"] == "absolute equality"
        crm = [_PMs.var(pm, n, :crm, r) for n in _PMs.nw_ids(pm) 
                                        for r in _PMs.ids(pm, :hsrc, nw=n) 
                                        if n ≠ fundamental(pm)]
    
        JuMP.@objective(pm.model, Max, sum(crm)) 
    end

    # maximin
    if pm.data["hhc_principle"] == "maximin"
        cmh = [_PMs.var(pm, n, :cmh) for n in _PMs.nw_ids(pm)
                                     if n ≠ fundamental(pm)]

        JuMP.@objective(pm.model, Max, sum(cmh))
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["hhc_principle"] == "Kalai-Smorodinsky bargaining"
        fh = [_PMs.var(pm, n, :fh) for n in _PMs.nw_ids(pm)
                                   if n ≠ fundamental(pm)]

        JuMP.@objective(pm.model, Max, sum(fh))
    end
end