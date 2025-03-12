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
function add_clean_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (nb, bus) in ntw["bus"] 
        if nw == "1" && bus["type"] == 4
            bus["v_fund_ref"] = 1.0
end end end

# constraints ##################################################################
## reference bus voltage constraint ############################################
""
function constraint_clean_voltage(pm::HarmonicPowerModel, i::Int; nw::Int=fundamental(pm))
    v_fund_ref  = _PMs.ref(pm, fundamental(pm), :bus, i, "v_fund_ref")

    nw == fundamental(pm) && constraint_clean_voltage_fundamental(pm, nw, i, v_fund_ref)
    nw ≠  fundamental(pm) && constraint_clean_voltage_harmonic(pm, nw, i)
end
""
function constraint_clean_voltage_fundamental(pm::HarmonicPowerModel, n::Int, i, v_fund_ref)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    JuMP.@constraint(pm.model, vbr == v_fund_ref)
    JuMP.@constraint(pm.model, vbi == 0.0)
end
""
function constraint_clean_voltage_harmonic(pm::HarmonicPowerModel, n::Int, i)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    JuMP.@constraint(pm.model, vbr == 0.0)
    JuMP.@constraint(pm.model, vbi == 0.0)
end