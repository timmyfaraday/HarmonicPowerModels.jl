################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

""
function sol_data_model!(pm::AbstractIVRModel, solution::Dict)
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
        for (i, xfmr) in solution["xfmr"]
            if haskey(xfmr, "pt_fr") && haskey(xfmr, "pt_to")
                xfmr["ptloss"] = xfmr["pt_fr"] + xfmr["pt_to"]
            end
            if haskey(xfmr, "qt_fr") && haskey(xfmr, "qt_to")
                xfmr["qtloss"] = xfmr["qt_fr"] + xfmr["qt_to"]
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