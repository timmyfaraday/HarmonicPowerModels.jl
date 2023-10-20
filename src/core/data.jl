################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

""
is_pos_sequence(nh::Int) = nh % 3 == 1
is_neg_sequence(nh::Int) = nh % 3 == 2
is_zero_sequence(nh::Int) = nh % 3 == 0

"""
    HarmonicPowerModels.extend_H!(H::Array{Int}, data::Dict{String, Any})

Returns all relevant harmonics `H` based on a PowerModels data dictionary `data`
and transformer magnetizing dictionary `xfmr_magn`.
"""
function extend_H!(H::Array{Int}, data::Dict{String, Any}, 
                   xfmr_magn::Dict{String, Any})
    Hbus = [parse(Int,ni[4:end]) for ni in keys(data["bus"]["1"]) if ni[1:2] == "nh"]
    push!(H, Hbus...)
    
    if haskey(xfmr_magn, "Hᴱ") push!(H, xfmr_magn["Hᴱ"]...) end
    if haskey(xfmr_magn, "Hᴵ") push!(H, xfmr_magn["Hᴵ"]...) end  

    unique!(sort!(H))
end

"""
    HarmonicPowerModels.replicate
"""
function _HPM.replicate(data::Dict{String, Any}; H::Array{Int}=Int[], 
                        xfmr_magn::Dict{String,Any}=Dict{String,Any}())
    # extend the user-provided harmonics `H` based on data input
    extend_H!(H, data, xfmr_magn)

    # create multi-network data structure, only keep the ntws whos id in H
    hdata = _PMs.replicate(data, last(H))
    for nh in 1:last(H) if nh ∉ H delete!(hdata["nw"],"$nh") end end

    # rename the case
    hdata["name"] = data["name"]

    # add the xfmr magnetizing current
    if !isempty(xfmr_magn)
        sample_magnetizing_current(hdata, xfmr_magn)
    end

    # re-evaluate the data for each harmonic 
    for (nw,ntw) in hdata["nw"]
        # translate the ntw-id to an Int for the harmonic number 
        nh = parse(Int,nw)

        # re-evaluate load power
        for load in values(ntw["load"])
            bus = ntw["bus"]["$(load["source_id"][2])"]
            
            mult = haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 0.0 ;
            
            haskey(load, "pd") ? load["pd"] *= mult : ~ ;
            haskey(load, "qd") ? load["qd"] *= mult : ~ ;
            load["multiplier"] = mult
        end

        # re-evaluate gen 
        for gen in values(ntw["gen"])
            if haskey(gen, "isfilter") && gen["isfilter"] == 1
                # do nothing, is handled by constraints constraint_active_filter
                gen["pmin"] = -abs(gen["pmax"])
                gen["qmin"] = -abs(gen["qmax"])
            else #is true generator
                if nw != "1" #cost of harmonics set to 0 
                    gen["cost"] *= 0 
                    #harmonics can be injected/absorbed to match load 
                    gen["pmin"] = -abs(gen["pmax"])
                    gen["qmin"] = -abs(gen["qmax"])
                end
            end
        end

        # re-evaluate the bus data 
        for bus in values(ntw["bus"])
            # use fundamental as limit for rms
            bus["vminrms"] = bus["vmin"]
            bus["vmaxrms"] = bus["vmax"]
            
            # true harmonics don't have minimum voltage
            if nw != "1" 
                bus["vmin"] = 0.0 
            end
        end

        # re-evaluate the branch data 
        for branch in values(ntw["branch"])
            haskey(branch, "br_r") ? branch["br_r"] *= sqrt(nh) : ~ ;
            haskey(branch, "br_x") ? branch["br_x"] *= nh : ~ ;
            haskey(branch, "b_fr") ? branch["b_fr"] *= nh : ~ ;
            haskey(branch, "b_to") ? branch["b_to"] *= nh : ~ ;
        end

        # re-evaluate the transformer data
        if haskey(ntw, "xfmr") 
            for xfmr in values(ntw["xfmr"])
                haskey(xfmr, "xsc") ? xfmr["xsc"] *= nh : ~ ;
                haskey(xfmr, "r1")  ? xfmr["r1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "r2")  ? xfmr["r2"] *= sqrt(nh) : ~ ;

                haskey(xfmr, "xe1") ? xfmr["xe1"] *= nh : ~ ;
                haskey(xfmr, "xe2") ? xfmr["xe2"] *= nh : ~ ;
                haskey(xfmr, "re1") ? xfmr["re1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "re2") ? xfmr["re2"] *= sqrt(nh) : ~ ;

                xfmr["cnf1"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][1]) : "Y" ;
                xfmr["cnf2"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][2]) : "Y" ;

                shift = haskey(xfmr, "vg") ? parse(Int, xfmr["vg"][3]) : 0 ;
                if is_pos_sequence(nh)
                    xfmr["tr"] = cosd(-30.0 * shift)
                    xfmr["ti"] = sind(-30.0 * shift)
                elseif is_neg_sequence(nh)
                    xfmr["tr"] = cosd(30.0 * shift)
                    xfmr["ti"] = sind(30.0 * shift)
                elseif is_zero_sequence(nh)
                    xfmr["tr"] = 1.0
                    xfmr["ti"] = 0.0
                end
            end
        end
    end

    return hdata
end