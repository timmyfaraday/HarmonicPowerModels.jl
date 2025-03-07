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

# util #########################################################################
function calc_max_hsrc_current(load::Dict{String,Any}, bus::Dict{String,Any})
    bus_id = load["load_bus"]
    return sqrt(load["pd"]^2 + load["qd"]^2) / bus["$bus_id"]["vmax"]
end

function hrsc_multiplier(load::Dict{String,Any}, ntw::Dict{String,Any}, nw::String)
    bus = ntw["bus"]["$(load["source_id"][2])"]
    nh = parse(Int,nw)
                
    return haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 1.0 ;
end

function hsrc_ref_angle(bus::Dict{String,Any}, nw::String)
    nh = parse(Int,nw)
    ref_angle = 0.0
    if haskey(bus, "ref_angle")
        if is_pos_sequence(nh)
            ref_angle = bus["ref_angle"]
        elseif is_neg_sequence(nh)
            ref_angle  = -bus["ref_angle"]
        elseif is_zero_sequence(nh)
            ref_angle  = 0.0
        end
    end 

    return ref_angle
end

# parameters ###################################################################
""
function add_source_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (l, load) in ntw["load"]
        if nw == "1"
        ntw["source"][l] = Dict( "source_bus" => load["load_bus"],
                                "pd"        => load["pd"],
                                "qd"        => load["qd"],
                                "multiplier"=> hrsc_multiplier(load, ntw, nw),
                                "ref_angle" => hsrc_ref_angle(ntw["bus"], nw))

        else
        ntw["source"][l] = Dict( "source_bus" => load["load_bus"],
                                "pd"        => load["pd"],
                                "qd"        => load["qd"],
                                "multiplier"=> hrsc_multiplier(load, ntw, nw),
                                "ref_angle" => hsrc_ref_angle(ntw["bus"], nw)),
                                "csmax"     => calc_max_hsrc_current(load, ntw["bus"])
        end 
    end 
end

# variables ####################################################################
""
function variable_hsrc_current(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_hsrc_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_hsrc_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    if nw ≠ fundamental(pm)
        variable_hsrc_current_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    end
end
""
function variable_hsrc_current_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csr = _PMs.var(pm, nw)[:csr] = JuMP.@variable(pm.model,
            [s in _PMs.ids(pm, nw, :source)], base_name="$(nw)_csr",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :source, s), "csr_start", 0.0)
    )

    if bounded
        for (s, source) in _PMs.ref(pm, nw, :source)
            csmax = source["csmax"]
            JuMP.set_lower_bound(csr[s], -csmax)
            JuMP.set_upper_bound(csr[s],  csmax)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :source, :csr, _PMs.ids(pm, nw, :source), csr)
end
""
function variable_hsrc_current_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csi = _PMs.var(pm, nw)[:csi] = JuMP.@variable(pm.model,
            [s in _PMs.ids(pm, nw, :source)], base_name="$(nw)_csi",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :source, s), "csi_start", 0.0)
    )

    if bounded
        for (s, source) in _PMs.ref(pm, nw, :source)
            csmax = source["csmax"]
            JuMP.set_lower_bound(csi[s], -csmax)
            JuMP.set_upper_bound(csi[s],  csmax)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :source, :csi, _PMs.ids(pm, nw, :source), csi)
end
""
function variable_hsrc_current_magnitude(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csm = _PMs.var(pm, nw)[:csm] = JuMP.@variable(pm.model,
            [s in _PMs.ids(pm, nw, :source)], base_name="$(nw)_csm",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :source, s), "csm_start", 0.0)
    )

    if bounded
        for (s, source) in _PMs.ref(pm, nw, :source)
            csmax = sourece["csmax"]
            JuMP.set_lower_bound(csm[s], 0.0)
            JuMP.set_upper_bound(csm[s], csmax)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :source, :csm, _PMs.ids(pm, nw, :source), csm)
end
""
## constraint harmonic source current ###############################################
function constraint_hsrc_current(pm::HarmonicPowerModel, s::Int; nw::Int=fundamental(pm))
    hsrc = _PMs.ref(pm, nw, :source, s)

    i      = hsrc["source_bus"]
    pd, qd = hsrc["pd"], hrsc["qd"]

    aref   = hsrc["ref_angle"]

    if nw == fundamental(pm)
        constraint_hsrc_constant_power(pm, nw, s, i, pd, qd)
    else
        constraint_hsrc_current_angle(pm, nw, s, aref)
    end  
end
""
function constraint_hsrc_constant_power(pm::HarmonicPowerModel, n::Int, s, i, pd, qd)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)
    csr = _PMs.var(pm, n, :csr, s)
    csi = _PMs.var(pm, n, :csi, s)

    JuMP.@constraint(pm.model, pd == vr*csr  + vi*csi)
    JuMP.@constraint(pm.model, qd == vi*csr  - vr*csi)
end
""
function constraint_hsrc_current_angle(pm::HarmonicPowerModel, n::Int, s, aref)
    csr = _PMs.var(pm, n, :csr, s)
    csi = _PMs.var(pm, n, :csi, s)
    csm = _PMs.var(pm, n, :csm, s)

    JuMP.@constraint(pm.model, csm * sind(aref) == csi)
    JuMP.@constraint(pm.model, csm * cosd(aref) == csr)
end
""
## constraint harmonic source power ###############################################
function constraint_hsrc_power(pm::HarmonicPowerModel, s::Int; nw::Int=fundamental(pm))
    hsrc = _PMs.ref(pm, nw, :source, s)

    i       = hsrc["source_bus"]
    pd, qd  = hsrc["pd"], hsrc["qd"]
    mult    = hsrc["multiplier"]

    if nw == 1
        constraint_hsrc_constant_power(pm, nw, s, i, pd, qd)
    else
        constraint_hrsc_constant_current(pm, nw, s, mult)
    end
end

""
function constraint_hsrc_constant_current(pm::_PMs.AbstractIVRModel, n::Int, s, mult)
    csr = _PMs.var(pm, n, :csr, s)
    csi = _PMs.var(pm, n, :csi, s)
    fund_csr = _PMs.var(pm, 1, :csr, s)
    fund_csi = _PMs.var(pm, 1, :csi, s)

    JuMP.@constraint(pm.model, csr == mult * fund_csr)
    JuMP.@constraint(pm.model, csi == mult * fund_csi)
end