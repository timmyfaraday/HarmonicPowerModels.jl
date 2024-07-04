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
    bus["vmax"] = 2.54
    bus["ref_angle"] = 0.0
    if b == "2136"
        bus["bus_type"] = 2
    end
end

for (g,gen) in data["gen"]
    if !haskey(gen, "gsc") 
        gen["gsc"] = 0.0
    end
    if !haskey(gen, "bsc") 
        gen["bsc"] = 0.0
    end
end



# define the set of considered harmonics
H = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]#, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49]

# solve HHC problem -- NLP
# hdata_nlp = HPM.replicate(data, H=H)
# results_hhc_nlp = HPM.solve_hhc(hdata_nlp, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)
for (n, nw) in hdata_soc["nw"]
    for (x, xfmr) in nw["xfmr"]
        xfmr["c_rating"] = xfmr["c_rating"] * 100
    end
end
total_time = @elapsed results_hhc = HPM.solve_hhc(hdata_soc, dHHC_SOC, solver_soc, solver_nlp)

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_", data["principle"],".json"]))
json_string = JSON.json(results_hhc)
open(filename,"w") do f
write(f, json_string)
end


##### Processing output #######


filename = joinpath(HPM.BASE_DIR, "results", join([case, "_maximum efficiency.json"]))
r_me = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_me = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_absolute equality.json"]))
r_ae = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ae = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_maximum efficiency.json"]))
r_me = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_me = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_maximin.json"]))
r_mm = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_mm = JSON.parse(dicttxt)
end


bus_voltage_magnitudes = zeros(length(H), length(data["bus"]))
bus_voltage_angles = zeros(length(H), length(data["bus"]))
harmonic_current_injections = zeros(length(H), length(data["load"]))

for h_idx in 2:length(H)
    h = H[h_idx]
    for (b, bus) in r_ae["solution"]["nw"]["$h"]["bus"]
        bus_voltage_magnitudes[h_idx, parse(Int, b)] = bus["vm"]
        bus_voltage_angles[h_idx, parse(Int, b)] = bus["va"]
    end
    for (l, load) in r_ae["solution"]["nw"]["$h"]["load"]
        harmonic_current_injections[h_idx, parse(Int, l)] = load["cm"]
    end
end

StatsPlots.boxplot(bus_voltage_magnitudes')
StatsPlots.boxplot(bus_voltage_angles')
StatsPlots.boxplot(harmonic_current_injections', outliers = false)

# for n in [3,5,7,9,13], l in 1:4
#     ca_nlp = rad2deg(results_hhc_nlp["solution"]["nw"]["$n"]["load"]["$l"]["ca"])
#     ca_soc = rad2deg(results_hhc_soc["solution"]["nw"]["$n"]["load"]["$l"]["ca"])
#     println("h=$n, l=$l: $(round(ca_nlp,digits=8)) vs $(round(ca_soc,digits=8))")
# end

# HPM.csv_export(results_hhc_nlp, csv_filename)


# for (g, gen) in data["gen"] if gen["bsc"] == nothing print(g, "\n") end end

# data_pf = deepcopy(hdata_soc)
# for (n, nw) in data_pf["nw"]
#     if n != "1"
#         delete!(data_pf["nw"], n)
#     end
# end
# hpf_results = HPM.solve_hpf(data_pf, PMs.IVRPowerModel, solver_nlp)


# for (b, bus) in hpf_results["solution"]["nw"]["1"]["bus"]
#     if  sqrt(bus["vr"]^2 + bus["vi"]^2) > 1.1
#     println(b, " ", bus["vr"], " ", bus["vi"], " ", sqrt(bus["vr"]^2 + bus["vi"]^2))
#     end
# end

# for (x, xfmr) in data["xfmr"]
#     if xfmr["f_bus"] == 640 || xfmr["t_bus"] == 640
#         println(x, " ", xfmr["f_bus"], " ", xfmr["t_bus"])
#     end
# end


# for (b, branch) in data["branch"]
#     if branch["f_bus"] == 640 || branch["t_bus"] == 640
#         println(b, " ", branch["f_bus"], " ", branch["t_bus"])
#     end
# end

# for (g, gen) in data["gen"]
#     if gen["gen_bus"] == 79 || gen["gen_bus"] == 80
#         println(g)
#     end
# end