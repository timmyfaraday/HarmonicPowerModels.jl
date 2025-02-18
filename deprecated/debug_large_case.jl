################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.1.0 - reviewed TVA                                                        #
# temp   - testing for cutting approach (TVA)                                  #
################################################################################

# using pkgs
using HarmonicPowerModels, PowerModels
using Gurobi
using Ipopt 
using JSON
using Plots
using StatsPlots

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_soc = Gurobi.Optimizer
solver_nlp = Ipopt.Optimizer

# read-in data 
# case = "pglib_opf_case2848_rte.m"
# case = "case1803_snem.m"
case = "nem2300harmonic"
extension = ".json"
csv_filename = joinpath(HPM.BASE_DIR,"results/voltage_phasors.csv")

if extension == ".m"
    path = joinpath(HPM.BASE_DIR,"test","data","matpower", join([case, extension]))
    data = PMs.parse_file(path)
elseif extension == ".json"
    path = joinpath(HPM.BASE_DIR,"test","data","json", join([case, extension]))
    data = Dict{String, Any}()
    open(path) do f
    dicttxt = read(f,String)  # file information to string
        global data = JSON.parse(dicttxt)  # parse and transform data
    end
else
    Memento.warn(_PM._LOGGER, "This input file is not yet supported!")
end


# Add principle
# {"maximum efficiency", "absolute equality", "maximin", "Kalai-Smorodinsky bargaining"}
data["principle"] = "maximum efficiency"

for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-2-4:2002, Cl. 2"
    bus["vmax"] = 1.1
    bus["ref_angle"] = 0.0
    if b == "2136"
        bus["bus_type"] = 2
    end
end

delete!(data["bus"],"79")
delete!(data["bus"],"80")
delete!(data["xfmr"],"332")
delete!(data["xfmr"],"333")
delete!(data["gen"],"41")
delete!(data["gen"],"42")

for (g,gen) in data["gen"]
    if !haskey(gen, "gsc") 
        gen["gsc"] = 0.0
    end
    if !haskey(gen, "bsc") 
        gen["bsc"] = 0.0
    end
end



# define the set of considered harmonics
H = [1, 3, 5, 7]#, 9, 11, 13, 15, 17, 19, 21]#, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49]



# solve HPF problem -- NLP
hdata_nlp = HPM.replicate(data, H=H)
for (n, nw) in hdata_nlp["nw"]
    for (x, xfmr) in nw["xfmr"]
        xfmr["c_rating"] = xfmr["c_rating"] * 100
    end
end
HPM.update_hdata_with_fundamental_hpf_results!(hdata_nlp, dHHC_NLP, solver_nlp)


for (b, bus) in hdata_nlp["nw"]["1"]["bus"]
    if bus["vm"] > 1.1
        println(b, " ",bus["vm"])
    end
end


for (l, load) in hdata_nlp["nw"]["1"]["gen"]
    if load["gen_bus"] == 80 || load["gen_bus"] == 79
        println(l)
    end
end