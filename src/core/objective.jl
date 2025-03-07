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

# variable fairness principle ################################################
""
function variable_fairness_principle(pm::_PMs.AbstractPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    # maximum efficiency
    if pm.data["principle"] == "maximum efficiency"
        # no additional variables
    end

    # absolute equality
    if pm.data["principle"] == "absolute equality"
        # no additional variables 
    end

    # maximin
    if pm.data["principle"] == "maximin"
        cmh = _PMs.var(pm, nw)[:cmh] = JuMP.@variable(pm.model, base_name="$(nw)_cmh",
                start = 0.0
        )

        if bounded
            JuMP.set_lower_bound(cmh, 0.0)
            # JuMP.set_upper_bound(cmd[d], c_rating)                            # @Hakan: dit is ook bij :cmd, van waar komt deze c rating, the fundamental component - wat als die er niet is?
        end

        report && (_PMs.sol(pm, nw, :fairness)[:cmh] = cmh)
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["principle"] == "Kalai-Smorodinsky bargaining"
        fh =  _PMs.var(pm, nw)[:fh] = JuMP.@variable(pm.model, base_name="$(nw)_fh",
                start = 0.0
        )

        if bounded
            JuMP.set_lower_bound(fh, 0.0)
            JuMP.set_upper_bound(fh, 1.0)
        end

        report && (_PMs.sol(pm, nw, :fairness)[:fh] = fh)
    end
end

# constraint fairness principle ################################################
""
function constraint_fairness_principle(pm::HarmonicPowerModel; nw::Int=fundamental(pm))
    source_ids = sort(collect(_PMs.ids(pm, :source nw=nw)))

    constraint_fairness_principle(pm, nw, source_ids)
end
""
function constraint_fairness_principle(pm::HarmonicPowerModel, n, source_ids)
    # maximum efficiency
    if pm.data["principle"] == "maximum efficiency"
        # no additional constraints
    end

    # absolute equality
    if pm.data["principle"] == "absolute equality"
        csm = [_PMs.var(pm, n, :csm, s) for s in source_ids]

        for s in source_ids[2:end]
            JuMP.@constraint(pm.model, csm[first(source_ids)] == csm[s]) ## to be CHECKED
        end 
    end

    # maximin
    if pm.data["principle"] == "maximin"
        cmh = _PMs.var(pm, n, :cmh)

        for s in source_ids
            csm = _PMs.var(pm, n, :csm, s)

            JuMP.@constraint(pm.model, cmh <= csm)
        end 
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["principle"] == "Kalai-Smorodinsky bargaining"
        fh = _PMs.var(pm, n, :fh)

        for s in source_ids
            csm = _PMs.var(pm, n, :csm, s)
            csmax = _PMs.ref(pm, n, :source, s, "csmax")

            JuMP.@constraint(pm.model, csm == fh * csmax)
        end 
    end
end

# objective harmonic power flow (hpf) ##############################################
""
function objective_power_flow(pm::_PMs.AbstractIVRModel)
    JuMP.@objective(pm.model, Min, 0.0)
end
# objective harmonic optimal power flow (hopf) #####################################
""
function objective_voltage_distortion_minimization(pm::_PMs.AbstractIVRModel) 
    bus_id = pm.data["bus_id"]

    vr = [_PMs.var(pm, n, :vr, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]
    vi = [_PMs.var(pm, n, :vi, bus_id) for n in _PMs.nw_ids(pm) if n ≠ 1]

    JuMP.@objective(pm.model, Min, sum(vr.^2 + vi.^2))
end
# objective harmonic hosting capacity (hhc) ########################################
""
function objective_maximum_hosting_capacity(pm::_PMs.AbstractIVRModel)
    # maximum efficiency
    if pm.data["principle"] == "maximum efficiency"
        csm = [_PMs.var(pm, n, :csm, s) for n in _PMs.nw_ids(pm) 
                                        for l in _PMs.ids(pm, :source, nw=n) 
                                        if n ≠ fundamental(pm)]
    
        JuMP.@objective(pm.model, Max, sum(csm))
    end

    # absolute equality
    if pm.data["principle"] == "absolute equality"
        csm = [_PMs.var(pm, n, :csm, l) for n in _PMs.nw_ids(pm) 
                                        for l in _PMs.ids(pm, :source, nw=n) 
                                        if n ≠ fundamental(pm)]
    
        JuMP.@objective(pm.model, Max, sum(csm)) 
    end

    # maximin
    if pm.data["principle"] == "maximin"
        cmh = [_PMs.var(pm, n, :cmh) for n in _PMs.nw_ids(pm)
                                     if n ≠ fundamental(pm)]

        JuMP.@objective(pm.model, Max, sum(cmh))
    end

    # Kalai-Smorodinsky bargaining
    if pm.data["principle"] == "Kalai-Smorodinsky bargaining"
        fh = [_PMs.var(pm, n, :fh) for n in _PMs.nw_ids(pm)
                                   if n ≠ fundamental(pm)]

        JuMP.@objective(pm.model, Max, sum(fh))
    end
end