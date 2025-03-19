################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# parameters ###################################################################
""
function add_fixed_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (nb, bus) in ntw["bus"] 
        if      nw == "1" && fdata["bus"][nb]["bus_type"] == 5
            bus["v_fund_ref"] = 1.0
        elseif  nw ≠ "1"  && fdata["bus"][nb]["bus_type"] == 5
            bus["v_harm_ref"] = 0.1
end end end

# constraints ##################################################################
## reference bus voltage constraint ############################################
""
function constraint_fixed_voltage(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    if nw == fundamental(pm) 
        v_fund_ref  = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_ref")

        constraint_fixed_voltage_fundamental(pm, nw, i, v_fund_ref)
    else
        v_harm_ref  = _PMs.ref(pm, nw, :bus, i, "v_harm_ref")

        constraint_fixed_voltage_harmonic(pm, nw, i, v_harm_ref)
end end
""
function constraint_fixed_voltage_fundamental(pm::HarmonicPowerModel, n::Int, i, v_fund_ref)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    JuMP.@constraint(pm.model, vbr == v_fund_ref)
    JuMP.@constraint(pm.model, vbi == 0.0)
end
""
function constraint_fixed_voltage_harmonic(pm::HarmonicPowerModel, n::Int, i, v_harm_ref)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    JuMP.@constraint(pm.model, vbr == v_harm_ref)
    JuMP.@constraint(pm.model, vbi == 0.0)
end