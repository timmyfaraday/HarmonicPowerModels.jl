################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
# v0.3.0 - adapted for extended graph representation                           #
################################################################################

# component list ###############################################################
const cmp_list = ["bus", "branch", "xfmr", "filter", "gen", "hload", "hsrc"]

# init #########################################################################
""
init_hdata_cmp(fdata::Dict{String,Any}, cmp::String, prob::Symbol) = 
    if haskey(fdata, cmp)
        Dict{String,Any}(nc => Dict{String,Any}() for nc in keys(fdata[cmp]))
    elseif cmp == "hload" && prob in [:hpf, :hopf] 
        Dict{String,Any}(nl => Dict{String,Any}() for nl in keys(fdata["load"]))
    elseif cmp == "hsrc" && prob in [:hhc] 
        Dict{String,Any}(nl => Dict{String,Any}() for nl in keys(fdata["load"]))
    else
        Dict{String,Any}()
    end
""
function init_hdata(fdata::Dict{String,Any}, H::Vector{Int}, prob::Symbol, bus_id)
    hdata = Dict{String,Any}(   "multinetwork"  => true,
                                "name"          => fdata["name"], 
                                "nw"            => Dict{String,Any}(),
                                "per_unit"      => fdata["per_unit"],
                                "s_base_mva"    => fdata["baseMVA"],
                                "bus_id"        => bus_id)

    for h in H
        hdata["nw"]["$h"] = Dict{String,Any}(cmp => init_hdata_cmp(fdata, cmp, prob)
                                                for cmp in cmp_list)
    end
    
    return hdata
end

# build from matpower file #####################################################
""
function build_hdata_from_matpower_file(fdata::Dict{String,Any}; 
                                        H::Vector{Int}=Int[1], 
                                        prob::Symbol=:hpf,
                                        bus_id::Int=1,
                                        xfmr_magn::Dict{String,Any}=Dict{String,Any}())
    hdata = init_hdata(fdata, H, prob, bus_id)

    add_bus_hdata!(hdata, fdata)
    add_ref_hdata!(hdata, fdata)
    add_clean_hdata!(hdata, fdata)

    add_branch_hdata!(hdata, fdata)
    add_xfmr_hdata!(hdata, fdata, xfmr_magn)

    add_filter_hdata!(hdata, fdata)
    add_gen_hdata!(hdata, fdata)
    prob in [:hpf, :hopf]   && add_hload_hdata!(hdata, fdata) 
    prob in [:hhc]          && add_hsrc_hdata!(hdata, fdata)

    return hdata
end