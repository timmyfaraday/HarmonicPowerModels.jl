################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Frederik Geth                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

""
function sol_data_model!(pm::_PMs.AbstractIVRModel, solution::Dict)
    _PMs.apply_pm!(_sol_data_model_ivr!, solution)
end

""
function _sol_data_model_ivr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vr") && haskey(bus, "vi")
                bus["vm"] = hypot(bus["vr"], bus["vi"])
                bus["va"] = atan(bus["vi"], bus["vr"])*180/pi
            end
        end
    end

    if haskey(solution, "branch")
        for (i, branch) in solution["branch"]
            if haskey(branch, "pf") && haskey(branch, "pt")
                branch["ploss"] = branch["pf"] + branch["pt"]
            end
            if haskey(branch, "qf") && haskey(branch, "qt")
                branch["qloss"] = branch["qf"] + branch["qt"]
            end
        end
    end

    if haskey(solution, "xfmr")
        for (x, xfmr) in solution["xfmr"]
            if haskey(xfmr, "px_fr") && haskey(xfmr, "px_to")
                xfmr["pxloss"] = xfmr["px_fr"] + xfmr["px_to"]
            end
            if haskey(xfmr, "qx_fr") && haskey(xfmr, "qx_to")
                xfmr["qxloss"] = xfmr["qx_fr"] + xfmr["qx_to"]
            end
        end
    end

    if haskey(solution, "gen")
        for (i, gen) in solution["gen"]
            if haskey(gen, "crg") && haskey(gen, "cig")
                gen["cm"] = hypot(gen["crg"], gen["cig"])
            end
        end
    end

    
    if haskey(solution, "load")
        for (i, load) in solution["load"]
            if haskey(load, "crd") && haskey(load, "cid")
                load["cm"] = hypot(load["crd"], load["cid"])
            end
        end
    end
end

""
function append_indicators!(result, hdata)
    solu = result["solution"]
    fundamental = "1"
    harmonics = Set(n for (n,nw) in solu["nw"])
    nonfundamentalharmonics = setdiff(harmonics, [fundamental])

    for (i,bus) in solu["nw"][fundamental]["bus"]
        vfun = solu["nw"][fundamental]["bus"][i]["vr"] + im*solu["nw"][fundamental]["bus"][i]["vi"]
        v   = [solu["nw"][n]["bus"][i]["vr"] + im*solu["nw"][n]["bus"][i]["vi"]  for n in harmonics]
        vnonfun = [solu["nw"][n]["bus"][i]["vr"] + im*solu["nw"][n]["bus"][i]["vi"]  for n in nonfundamentalharmonics]
        
        rms = sqrt(sum(abs.(v).^2))
        thd = sqrt(sum(abs.(vnonfun).^2))/abs(vfun)

        delete!(solu["nw"][fundamental]["bus"][i], "w")
        delete!(solu["nw"][fundamental]["bus"][i], "vr")
        delete!(solu["nw"][fundamental]["bus"][i], "vi")

        solu["nw"][fundamental]["bus"][i]["vrms"] = rms
        solu["nw"][fundamental]["bus"][i]["vthd"] = thd
        solu["nw"][fundamental]["bus"][i]["vmaxrms"] = hdata["nw"][fundamental]["bus"][i]["vmaxrms"]
        solu["nw"][fundamental]["bus"][i]["vminrms"] = hdata["nw"][fundamental]["bus"][i]["vminrms"]
        solu["nw"][fundamental]["bus"][i]["thdmax"] = hdata["nw"][fundamental]["bus"][i]["thdmax"]

    end

    for (i,gen) in solu["nw"][fundamental]["gen"]
        gencost = hdata["nw"][fundamental]["gen"][i]["cost"]
        pg = gen["pg"]
        if length(gencost) >0
            gen["totcost"] = sum(gencost[c]*(pg)^c for c in 1:length(gencost))
        end

        cfun = solu["nw"][fundamental]["gen"][i]["crg"] + im*solu["nw"][fundamental]["gen"][i]["cig"]
        call   = [solu["nw"][n]["gen"][i]["crg"] + im*solu["nw"][n]["gen"][i]["cig"]  for n in harmonics]
        cnonfun = [solu["nw"][n]["gen"][i]["crg"] + im*solu["nw"][n]["gen"][i]["cig"]  for n in nonfundamentalharmonics]

        rms = sqrt(sum(abs.(call).^2))
        thd = sqrt(sum(abs.(cnonfun).^2))/abs(cfun)
        
        solu["nw"][fundamental]["gen"][i]["crms"] = rms
        solu["nw"][fundamental]["gen"][i]["cthd"] = thd
        solu["nw"][fundamental]["gen"][i]["c_rating"] = hdata["nw"][fundamental]["gen"][i]["c_rating"]
        solu["nw"][fundamental]["gen"][i]["pmax"] = hdata["nw"][fundamental]["gen"][i]["pmax"]
    end

    if haskey(solu["nw"][fundamental], "load")
        for (i,load) in solu["nw"][fundamental]["load"]

            cfun = solu["nw"][fundamental]["load"][i]["crd"] + im*solu["nw"][fundamental]["load"][i]["cid"]
            call   = [solu["nw"][n]["load"][i]["crd"] + im*solu["nw"][n]["load"][i]["cid"]  for n in harmonics]
            cnonfun = [solu["nw"][n]["load"][i]["crd"] + im*solu["nw"][n]["load"][i]["cid"]  for n in nonfundamentalharmonics]

            rms = sqrt(sum(abs.(call).^2))
            thd = sqrt(sum(abs.(cnonfun).^2))/abs(cfun)
            
            solu["nw"][fundamental]["load"][i]["crms"] = rms
            solu["nw"][fundamental]["load"][i]["cthd"] = thd
        end 
    end
end