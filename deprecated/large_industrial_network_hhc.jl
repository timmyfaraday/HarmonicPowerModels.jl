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
using Revise
using HarmonicPowerModels, PowerModels
using Gurobi 
using Ipopt
using JSON

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_nlp = Ipopt.Optimizer
solver_soc = Gurobi.Optimizer

# read-in data 
case    = "rain"
ext     = ".json"
path    = joinpath(HPM.BASE_DIR,"test","data","json", join([case, ext]))

data = Dict{String, Any}()
open(path) do f
dicttxt = read(f,String)                                                        # file information to string
    global data = JSON.parse(dicttxt)                                           # parse and transform data
end

# add name
data["name"] = case
data["per_unit"] = true
data["baseMVA"] = 100.0

# add principle
data["principle"] = "absolute equality"

# update of the data
for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-2-4:2002, Cl. 2"
    # bus["vmax"] = 2.54
    bus["ref_angle"] = 0.0
    bus["status"] = 1
end

for (b, branch) in data["branch"]
    branch["br_status"] = 1
    branch["tap"] = 1
end

for (l, load) in data["load"]
    load["status"] = 1
end

data["shunt"] = Dict()
data["storage"] = Dict()
data["switch"] = Dict()
data["dcline"] = Dict()

data["gen"] = Dict()
data["gen"][1] = Dict("gen_status" => 1,
                        "source_id" => Any["gen", 1],
                        "index" => 1,
                        "gen_bus" => 13,
                        "cost" => 0.0,
                        "pmin" => 0.0,
                        "pmax" => 1000.0,
                        "qmin" => -1000.0,
                        "qmax" => 1000.0)
for (g,gen) in data["gen"]
    if !haskey(gen, "gsc") 
        gen["gsc"] = 0.0
    end
    if !haskey(gen, "bsc") 
        gen["bsc"] = 0.0
    end
end

# define the set of considered harmonics
H = [1, 3, 5, 7]

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)
results_hhc = HPM.solve_hhc(hdata_soc, dHHC_SOC, solver_soc, solver_nlp)

