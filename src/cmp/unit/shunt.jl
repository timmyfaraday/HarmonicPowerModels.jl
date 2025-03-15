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
function add_shunt_hdata!(hdata::Dict{String,Any}, fdata::Dict{String,Any})
    for (nw, ntw) in hdata["nw"], (s, shunt) in ntw["shunt"]
            h       = parse(Int, nw)
            sdata   = fdata["shunt"][s]
            if nw == "1"
                shunt["id"]         = sdata["index"]
                shunt["bus"]        = sdata["shunt_bus"]
                #-----------------------------------#
                shunt["b"]          = sdata["bs"]
                shunt["g"]          = sdata["gs"]
            else
                shunt["b"]          = sdata["bs"] * h^(sign(sdata["bs"]))
                shunt["g"]          = sdata["gs"] / sqrt(h)
end end end
# variable #####################################################################
""
function variable_shunt_current(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true, kwargs...)
    variable_shunt_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_shunt_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
""
function variable_shunt_current_real(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csr = _PMs.var(pm, nw)[:csr] = 
            JuMP.@variable( pm.model,
                            [s in _PMs.ids(pm, nw, :shunt)], 
                            base_name="$(nw)_csr",
                            start=0.0)

    report && _PMs.sol_component_value(pm, nw, :shunt, :csr, _PMs.ids(pm, nw, :shunt), csr)
end
""
function variable_shunt_current_imaginary(pm::HarmonicPowerModel; nw::Int=fundamental(pm), bounded::Bool=true, report::Bool=true)
    csi = _PMs.var(pm, nw)[:csi] = 
            JuMP.@variable( pm.model,
                            [s in _PMs.ids(pm, nw, :shunt)], 
                            base_name="$(nw)_csi",
                            start=0.0)

    report && _PMs.sol_component_value(pm, nw, :shunt, :csi, _PMs.ids(pm, nw, :shunt), csi)
end

# constraints ##################################################################
""
function constraint_shunt_current(pm::HarmonicPowerModel, s::Int; nw::Int=fundamental(pm))
    i   = _PMs.ref(pm, fundamental(pm), :shunt, s, "bus")

    g   = _PMs.ref(pm, nw, :shunt, s, "g")
    b   = _PMs.ref(pm, nw, :shunt, s, "b")

    if nw ≠ fundamental(pm)
        constraint_shunt_current(pm, nw, s, i, g, b)
end end
""
function constraint_shunt_current(pm::HarmonicPowerModel, n::Int, s, i, g, b)
    vbr = _PMs.var(pm, n, :vbr, i)
    vbi = _PMs.var(pm, n, :vbi, i)

    csr = _PMs.var(pm, n, :csr, s)
    csi = _PMs.var(pm, n, :csi, s)

    JuMP.@constraint(pm.model, csr == g * vbr - b * vbi)
    JuMP.@constraint(pm.model, csi == g * vbi + b * vbr)
end
