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
const cmp_list = ["bus", "branch", "xfmr", "hsrc", "gen"]

# init #########################################################################
""
init_hdata_cmp(fdata::Dict{String,Any}, cmp::String) = 
    Dict{String,Any}(nc => Dict{String,Any} for nc in keys(fdata[cmp]))
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
end end

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
    add_hscr_hdata!(hdata, fdata)
end

# ""
# function _HPM.replicate(data::Dict{String, Any}; 
#                             bus_id::Int=1,
#                             H::Array{Int}=Int[], 
#                             xfmr_magn::Dict{String,Any}=Dict{String,Any}())
    

#     ### Add entries to the fundamental data ####################################
#     # set the branch current rating
#     # if haskey(data,"branch") 
#     #     for branch in values(data["branch"])
#     #         f_bus = data["bus"][string(branch["f_bus"])]
#     #         t_bus = data["bus"][string(branch["t_bus"])]
            
#     #         vmmin = min(f_bus["vmin"], t_bus["vmin"])
            
#     #         #branch["c_rating"] = branch["rate_a"] / sqrt(3) #/ vmmin             # @Hakan: klopt dit
#     # end end

#     # # add xfmr current rating
#     # if haskey(data,"xfmr") 
#     #     for xfmr in values(data["xfmr"])
#     #         f_bus = data["bus"][string(xfmr["f_bus"])]
#     #         t_bus = data["bus"][string(xfmr["t_bus"])]
            
#     #         vmmin = min(f_bus["vmin"], t_bus["vmin"])
            
#     #         # xfmr["c_rating"] = xfmr["rateA"] / data["baseMVA"] / vmmin        # @Hakan: klopt dit
#     # end end

#     # # add the thd limits based on standard, if available
#     # for bus in values(data["bus"])
#     #     if haskey(bus, "standard") && bus["standard"] in keys(thd_limits)
#     #         bus["thdmax"] = thd_limits[bus["standard"]]
#     # end end

#     ### create multi-network data structure, only keep ntws in H ###############
#     hdata = _PMs.replicate(data, last(H))
#     for nh in 1:last(H) if nh ∉ H delete!(hdata["nw"],"$nh") end end

#     # add bus_id
#     hdata["bus_id"] = bus_id

#     # rename the case
#     hdata["name"] = data["name"]
#     haskey(data, "principle") ? hdata["principle"] = data["principle"] : ~ ;

#     ### add entries to the harmonic data #######################################
#     # add the xfmr magnetizing current
#     if !isempty(xfmr_magn)
#         sample_magnetizing_current(hdata, xfmr_magn)
#     end

#     # re-evaluate the data for each harmonic 
#     for (nw,ntw) in hdata["nw"]
#         # translate the ntw-id to an Int for the harmonic number 
#         nh = parse(Int,nw)

#         # re-evaluate load power
#         for load in values(ntw["load"])
#             bus = ntw["bus"]["$(load["source_id"][2])"]
                
#             mult = haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 1.0 ;

#             haskey(load, "pd") ? load["pd"] *= mult : ~ ;
#             haskey(load, "qd") ? load["qd"] *= mult : ~ ;
#             load["multiplier"] = mult

#             load["c_rating"] = 1.0 # TODO

#             if haskey(bus, "ref_angle")
#                 if is_pos_sequence(nh)
#                     load["ref_angle"] = bus["ref_angle"]
#                 elseif is_neg_sequence(nh)
#                     load["ref_angle"] = -bus["ref_angle"]
#                 elseif is_zero_sequence(nh)
#                     load["ref_angle"] = 0.0
#             end end

#             if haskey(bus, "angle_range")
#                 load["harmonic_angle_range"] = bus["angle_range"]
#             end
#         end

#         # re-evaluate gen 
#         for gen in values(ntw["gen"])
#             if haskey(gen, "isfilter") && gen["isfilter"] == 1
#                 # do nothing, is handled by constraints constraint_active_filter
#                 gen["pmin"] = -abs(gen["pmax"])
#                 gen["qmin"] = -abs(gen["qmax"])
#             else #is true generator
#                 if nw != "1" #cost of harmonics set to 0 
#                     if haskey(gen, "cost")
#                         gen["cost"] *= 0 
#                     end
#                     #harmonics can be injected/absorbed to match load 
#                     gen["pmin"] = -abs(gen["pmax"])
#                     gen["qmin"] = -abs(gen["qmax"])
#             end end
#             gen_bus = gen["gen_bus"]
#             #gen["c_rating"] = sqrt(max(abs(gen["pmin"]), abs(gen["pmax"]))^2 + max(abs(gen["qmin"]), abs(gen["qmax"]))^2) / (sqrt(3) * ntw["bus"]["$gen_bus"]["vmin"])
#         end

#         # re-evaluate the bus data 
#         for bus in values(ntw["bus"])
#             # use fundamental as limit for rms
#             bus["vminrms"] = bus["vmin"]
#             bus["vmaxrms"] = bus["vmax"]
            
#             # true harmonics don't have minimum voltage
#             if nw != "1" 
#                 bus["vmin"] = 0.0 
#             end

#             # add the ihd limits based on standard, if available
#             if haskey(bus, "standard") 
#                 std = bus["standard"]
#                 if std in keys(ihd_limits)
#                     if nh <= length(ihd_limits[std])
#                         bus["ihdmax"] = ihd_limits[std][nh]
#                     else
#                         @warn "harmonic $nh not included in $std"
#         end end end end

#         # re-evaluate the branch data 
#         for branch in values(ntw["branch"])
#             haskey(branch, "br_r") ? branch["br_r"] *= sqrt(nh) : ~ ;
#             haskey(branch, "br_x") ? branch["br_x"] *= nh : ~ ;
#             haskey(branch, "b_fr") ? branch["b_fr"] *= nh : ~ ;
#             haskey(branch, "b_to") ? branch["b_to"] *= nh : ~ ;
#         end

#         # re-evaluate the gen data
#         for gen in values(ntw["gen"])
#             haskey(gen, "bsc") ? gen["bsc"] /= nh : ~ ;
#             haskey(gen, "gsc") ? gen["gsc"] /= sqrt(nh) : ~ ;
#         end

#         # re-evaluate the transformer data
#         if haskey(ntw, "xfmr") 
#             for xfmr in values(ntw["xfmr"])
#                 haskey(xfmr, "xsc") ? xfmr["xsc"] *= nh : ~ ;
#                 haskey(xfmr, "r1")  ? xfmr["r1"] *= sqrt(nh) : ~ ;
#                 haskey(xfmr, "r2")  ? xfmr["r2"] *= sqrt(nh) : ~ ;

#                 haskey(xfmr, "xe1") ? xfmr["xe1"] *= nh : ~ ;
#                 haskey(xfmr, "xe2") ? xfmr["xe2"] *= nh : ~ ;
#                 haskey(xfmr, "re1") ? xfmr["re1"] *= sqrt(nh) : ~ ;
#                 haskey(xfmr, "re2") ? xfmr["re2"] *= sqrt(nh) : ~ ;

#                 xfmr["cnf1"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][1]) : 'Y' ;
#                 xfmr["cnf2"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][2]) : 'Y' ;

#                 shift = haskey(xfmr, "vg") ? parse(Int, xfmr["vg"][3:end]) : 0 ;
#                 if is_pos_sequence(nh)
#                     xfmr["tr"] = cosd(-30.0 * shift)
#                     xfmr["ti"] = sind(-30.0 * shift)
#                 elseif is_neg_sequence(nh)
#                     xfmr["tr"] = cosd(30.0 * shift)
#                     xfmr["ti"] = sind(30.0 * shift)
#                 elseif is_zero_sequence(nh)
#                     xfmr["tr"] = 1.0
#                     xfmr["ti"] = 0.0
#         end end end
#     end 

#     return hdata
# end
