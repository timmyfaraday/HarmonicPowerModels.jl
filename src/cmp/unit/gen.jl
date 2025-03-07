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
""
function calc_gen_conductance(gen::Dict{String, Any}, nw::String)
    nh = parse(Int, nw)
    rx_ratio = 1/20.0
    return (1 / (sqrt(gen["pmax"]^2 + gen["qmax"]^2) * rx_ratio)) / sqrt(nh)
end
""
function calc_gen_admittance(gen::Dict{String, Any}, nw::String)
    nh = parse(Int, nw)
    return (1 / sqrt(gen["pmax"]^2 + gen["qmax"]^2)) / nh
end
""
function calculate_max_gen_current(gen::Dict{String, Any})
    return sqrt(gen["pmax"]^2 + gen["qmax"]^2)
end
""
function calculate_fundamental_gen_current(gen::Dict{String, Any}, ntw::Dict{String, Any})
    gen_bus = gen["gen_bus"]
    bus_voltage = ntw["bus"]["$gen_bus"]["vm"]
    return sqrt(gen["pmax"]^2 + gen["qmax"]^2) / bus_voltage
end
# parameters ###################################################################
""
function add_gen_hdata(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (g, gen) in ntw["gen"]
        b
        if nw == "1"
            ntw["gen"][g] = Dict(  "gen_bus" => gen["gen_bus"],
                                   "pg" => gen["pg"],
                                   "qg" => gen["qg"],
                                   "cost" => gen["cost"],
                                   "pmax" => gen["pmax"],
                                   "pmin" => gen["pmin"],
                                   "qmax" => gen["qmax"],
                                   "qmin" => gen["qmin"],
                                   "ncost" => gen["ncost"],
                                   "inf" => 0,                               # check what to do with the infinite bus
                                   "gsc" => calc_gen_conductance(gen, nw),
                                   "bsc" => calc_gen_admittance(gen, nw),
                                   "c_rating" => calculate_max_gen_current(gen),
                                   "cm" => calculate_fundamental_gen_current(gen, ntw),
                                   "gen_status" => gen["gen_status"])
        
        else
            ntw["gen"][g] = Dict( "gen_bus" => gen["gen_bus"],
                                  "inf" => 0,                               # check what to do with the infinite bus
                                  "gsc" => calc_gen_conductance(gen, nw),
                                  "bsc" => calc_gen_admittance(gen, nw),
                                  "c_rating" => calculate_gen_current(gen))
        end 
    end
end

# variables ####################################################################
""
function variable_gen_current_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cgr = _PMs.var(pm, nw)[:cgr] = JuMP.@variable(pm.model,
            [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_cgr",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :gen, g), "cgr_start", 0.0)
    )

    if bounded
        for (g, gen) in _PMs.ref(pm, nw, :gen)
            c_rating = gen["c_rating"]
            JuMP.set_lower_bound(cgr[g], -c_rating)
            JuMP.set_upper_bound(cgr[g],  c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :cgr, _PMs.ids(pm, nw, :gen), cgr)
end
""
function variable_gen_current_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    cgi = _PMs.var(pm, nw)[:cgi] = JuMP.@variable(pm.model,
            [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_cgi",
            start = _PMs.comp_start_value(_PMs.ref(pm, nw, :gen, g), "cgi_start", 0.0)
    )

    if bounded
        for (g, gen) in _PMs.ref(pm, nw, :gen)
            c_rating = gen["c_rating"]
            JuMP.set_lower_bound(cgi[g], -c_rating)
            JuMP.set_upper_bound(cgi[g],  c_rating)
        end
    end

    report && _PMs.sol_component_value(pm, nw, :gen, :cgi, _PMs.ids(pm, nw, :gen), cgi)
end

# contraint generator current ####################################################################
""
function constraint_gen_current(pm::HarmonicPowerModel, g::Int; nw::Int=fundamental(pm))
    gen = _PMs.ref(pm, fundamental(pm), :gen, g)
    bus = gen["gen_bus"]

    inf = _PMs.ref(pm, nw, :gen, g, "inf")
    
    gsc = _PMs.ref(pm, nw, :gen, g, "gsc")
    bsc = _PMs.ref(pm, nw, :gen, g, "bsc")

    if iszero(inf) && nw ≠ fundamental
        constraint_gen_current(pm, nw, g, bus, gsc, bsc)
    end
end
""
function constraint_gen_current(pm::HarmonicPowerModel, n::Int, g, i, gsc, bsc)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    cgr = _PMs.var(pm, n, :cgr, g)
    cgi = _PMs.var(pm, n, :cgi, g)

    JuMP.@constraint(pm.model, cgr == gsc*vr - bsc*vi)
    JuMP.@constraint(pm.model, cgi == gsc*vi + bsc*vr)
end
""
# contraint generator current rms limit ###########################################################
""
function constraint_gen_current_rms_limit(pm::HarmonicPowerModel, g::Int)
    gen = _PMs.ref(pm, fundamental(pm), :gen, g)

    c_rating = gen["c_rating"]

    constraint_gen_current_rms_limit(pm, g, c_rating)
end
function constraint_gen_current_rms_limit(pm::HarmonicPowerModel, g, c_rating)
    cgr =  [_PMs.var(pm, n, :cgr, g) for n in sorted_nw_ids(pm)]
    cgi =  [_PMs.var(pm, n, :cgi, g) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model, sum(cgr.^2 + cgi.^2) <= c_rating^2)
end
""
# contraint generator current rms limit in SOC #####################################################
""
function constraint_gen_current_rms_limit(pm::dHHCPowerModel, g::Int)
    gen = _PMs.ref(pm, fundamental(pm), :gen, g)
    
    cm_fund = gen["cm"]
    c_rating = gen["c_rating"]

    constraint_gen_current_rms_limit(pm, g, c_rating, cm_fund)
end
""
function constraint_gen_current_rms_limit(pm::dHHCPowerModel, g, c_rating, cm_fund)
    cgr =  [_PMs.var(pm, n, :cgr, g) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]
    cgi =  [_PMs.var(pm, n, :cgi, g) for n in sorted_nw_ids(pm) if n ≠ fundamental(pm)]

    JuMP.@constraint(pm.model, [sqrt(c_rating^2 - cm_fund^2)./1000; vcat(cgr, cgi)./1000] in JuMP.SecondOrderCone())  # Fix the RHS scaling with extra parameter later
end






