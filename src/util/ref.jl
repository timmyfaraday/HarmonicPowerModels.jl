################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
################################################################################

# filter
""
function ref_add_filter!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    _PMs.apply_pm!(_ref_add_filter!, ref, data)
end
""
function _ref_add_filter!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if !haskey(ref, :filter)
        ref[:filter] = Dict()
        ref[:bus_filters] = Dict((i, []) for (i,bus) in ref[:bus])
    else
        ref[:filter] = Dict(f for f in ref[:filter] if f.second["bus"] in keys(ref[:bus]))

        bus_filters = Dict((i, Int[]) for (i,bus) in ref[:bus])
        for (i,filter) in ref[:filter]
            push!(bus_filters[filter["bus"]], i)
        end
        ref[:bus_filters] = bus_filters
    end
end

# xfmr
""
function ref_add_xfmr!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    _PMs.apply_pm!(_ref_add_xfmr!, ref, data)
end
""
function _ref_add_xfmr!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if !haskey(ref, :xfmr)
        ref[:xfmr] = Dict()
        ref[:xfmr_arcs_from] = Dict()
        ref[:xfmr_arcs_to] = Dict()
        ref[:xfmr_arcs] = Dict()
        ref[:bus_arcs_xfmr] = Dict((i, []) for (i,bus) in ref[:bus])
    else
        ref[:xfmr] = Dict(x for x in ref[:xfmr] if  x.second["f_bus"] in keys(ref[:bus]) &&
                                                    x.second["t_bus"] in keys(ref[:bus]))
        
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