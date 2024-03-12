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
################################################################################

""
const ihd_limits = Dict(
    "Clean Bus" =>                 [1.0, zeros(49)...],
    "IEC61000-2-4:2002, Cl. 2" =>  [1.00000, 0.02000, 0.05000, 0.01000, 0.06000, 
                                    0.00500, 0.05000, 0.00500, 0.01500, 0.00500, 
                                    0.03500, 0.00458, 0.03000, 0.00429, 0.00400, 
                                    0.00406, 0.02000, 0.00389, 0.01761, 0.00375, 
                                    0.00200, 0.00364, 0.01408, 0.00354, 0.01275, 
                                    0.00346, 0.00200, 0.00339, 0.01061, 0.00333, 
                                    0.00975, 0.00328, 0.00200, 0.00324, 0.00833, 
                                    0.00319, 0.00773, 0.00316, 0.00200, 0.00313, 
                                    0.00671, 0.00310, 0.00627, 0.00307, 0.00200, 
                                    0.00304, 0.00551, 0.00302, 0.00518, 0.00300],
    "IEC61000-3-6:2008" =>         [1.00000, 0.01400, 0.02000, 0.00800, 0.02000,
                                    0.00400, 0.02000, 0.00400, 0.01000, 0.00350,
                                    0.01500, 0.00318, 0.01500, 0.00296, 0.00300,
                                    0.00279, 0.01200, 0.00266, 0.01074, 0.00255,
                                    0.00200, 0.00246, 0.00887, 0.00239, 0.00816,
                                    0.00233, 0.00200, 0.00228, 0.00703, 0.00223,
                                    0.00658, 0.00219, 0.00200, 0.00216, 0.00583,
                                    0.00213, 0.00551, 0.00210, 0.00200, 0.00208,
                                    0.00498, 0.00205, 0.00474, 0.00203, 0.00200,
                                    0.00201, 0.00434, 0.00200, 0.00416, 0.00198])

const thd_limits = Dict(
    "Clean Bus" =>                  0.00000,
    "IEC61000-2-4:2002, Cl. 2" =>   0.08000,
    "IEC61000-3-6:2008" =>          0.08000) 

"""
    HarmonicPowerModels.replicate
"""
function _HPM.replicate(data::Dict{String, Any}; 
                            bus_id::Int=1,
                            H::Array{Int}=Int[], 
                            xfmr_magn::Dict{String,Any}=Dict{String,Any}())
    ### Add entries to the fundamental data ####################################
    # set the branch current rating
    if haskey(data,"branch") 
        for branch in values(data["branch"])
            f_bus = data["bus"][string(branch["f_bus"])]
            t_bus = data["bus"][string(branch["t_bus"])]
            
            vmmin = min(f_bus["vmin"], t_bus["vmin"])
            
            branch["c_rating"] = branch["rate_a"] / sqrt(3) / vmmin  
    end end

    # add xfmr current rating
    if haskey(data,"xfmr") 
        for xfmr in values(data["xfmr"])
            f_bus = data["bus"][string(xfmr["f_bus"])]
            t_bus = data["bus"][string(xfmr["t_bus"])]
            
            vmmin = min(f_bus["vmin"], t_bus["vmin"])
            
            xfmr["c_rating"] = xfmr["rateA"] / data["baseMVA"] / vmmin          # @Hakan: klopt dit
    end end

    # add the thd limits based on standard, if available
    for bus in values(data["bus"])
        if haskey(bus, "standard") && bus["standard"] in keys(thd_limits)
            bus["thdmax"] = thd_limits[bus["standard"]]
    end end

    ### create multi-network data structure, only keep ntws in H ###############
    hdata = _PMs.replicate(data, last(H))
    for nh in 1:last(H) if nh ∉ H delete!(hdata["nw"],"$nh") end end

    # add bus_id
    hdata["bus_id"] = bus_id

    # rename the case
    hdata["name"] = data["name"]
    haskey(data, "principle") ? hdata["principle"] = data["principle"] : ~ ;

    ### add entries to the harmonic data #######################################
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
                
            if nh == 1
                mult = haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 1.0 ;
            else
                mult = haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 1.0 ;
            end

            haskey(load, "pd") ? load["pd"] *= mult : ~ ;
            haskey(load, "qd") ? load["qd"] *= mult : ~ ;
            load["multiplier"] = mult

            load["c_rating"] = 1.0 # TODO

            if haskey(bus, "ref_angle")
                if is_pos_sequence(nh)
                    load["ref_angle"] = bus["ref_angle"]
                elseif is_neg_sequence(nh)
                    load["ref_angle"] = -bus["ref_angle"]
                elseif is_zero_sequence(nh)
                    load["ref_angle"] = 0.0
            end end

            if haskey(bus, "rng_angle")
                load["rng_angle"] = bus["rng_angle"]
            end
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
            end end
            gen_bus = gen["gen_bus"]
            gen["c_rating"] = sqrt(max(abs(gen["pmin"]), abs(gen["pmax"]))^2 + max(abs(gen["qmin"]), abs(gen["qmax"]))^2) / (sqrt(3) * ntw["bus"]["$gen_bus"]["vmin"])
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

            # add the ihd limits based on standard, if available
            if haskey(bus, "standard") 
                std = bus["standard"]
                if std in keys(ihd_limits)
                    if nh <= length(ihd_limits[std])
                        bus["ihdmax"] = ihd_limits[std][nh]
                    else
                        println("harmonic $nh not included in $std")            # change to warn 
            end end end 
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

                xfmr["cnf1"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][1]) : 'Y' ;
                xfmr["cnf2"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][2]) : 'Y' ;

                shift = haskey(xfmr, "vg") ? parse(Int, xfmr["vg"][3:end]) : 0 ;
                if is_pos_sequence(nh)
                    xfmr["tr"] = cosd(-30.0 * shift)
                    xfmr["ti"] = sind(-30.0 * shift)
                elseif is_neg_sequence(nh)
                    xfmr["tr"] = cosd(30.0 * shift)
                    xfmr["ti"] = sind(30.0 * shift)
                elseif is_zero_sequence(nh)
                    xfmr["tr"] = 1.0
                    xfmr["ti"] = 0.0
        end end end
    end 

    return hdata
end
