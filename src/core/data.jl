################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
# v0.3.0 - adapted for extended graph representation                           #
################################################################################

# component list ###############################################################
const cmp_list = ["bus", "branch", "xfmr", "hload", "hsrc", "gen"]

# init #########################################################################
""
init_hdata_cmp(fdata::Dict{String,Any}, cmp::String) = 
    if haskey(fdata, cmp)
        Dict{String,Any}(nc => Dict{String,Any}() for nc in keys(fdata[cmp]))
    else
        Dict{String,Any}()
    end
""
function init_hdata(fdata::Dict{String,Any}, H::Vector{Int})
    hdata = Dict{String,Any}(   "multinetwork"  => true,
                                "name"          => fdata["name"], 
                                "nw"            => Dict{String,Any}(),
                                "per_unit"      => fdata["per_unit"],
                                "s_base_mva"    => fdata["baseMVA"])

    for h in H
        hdata["nw"]["$h"] = Dict{String,Any}(cmp => init_hdata_cmp(fdata, cmp)
                                                for cmp in cmp_list)
    end
    
    return hdata
end

# build from matpower file #####################################################
""
function build_hdata_from_matpower_file(fdata::Dict{String,Any}; 
                                        H::Vector{Int}=Int[1], 
                                        xfmr_magn::Dict{String,Any}=Dict{String,Any}())
    hdata = init_hdata(fdata, H)

    add_bus_hdata!(hdata, fdata)
    add_ref_hdata!(hdata, fdata)

    add_branch_hdata!(hdata, fdata)
    add_xfmr_hdata!(hdata, fdata, xfmr_magn)

    add_gen_hdata!(hdata, fdata)
    add_hload_hdata!(hdata, fdata)
    add_hscr_hdata!(hdata, fdata)

    return hdata
end