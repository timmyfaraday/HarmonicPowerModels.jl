################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

""
function ref_add_xfmr!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any}) ## data not actually needed!
    _PMs.apply_pm!(_ref_add_xfmr!, ref, data)
end
""
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