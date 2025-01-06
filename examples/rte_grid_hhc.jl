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
using SparseArrays

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels 

# Include function to calculate Harmonic impedance
include(joinpath(HPM.BASE_DIR,"src/util/hi.jl"))

# set the solver
solver_soc = Gurobi.Optimizer
solver_nlp = Ipopt.Optimizer

# read-in data 
case = "case1888_rte"
extension = ".m"
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
    Memento.warn(_PM._LOGGER, "This input file format is not yet supported!")
end


# Flag if you want to calculate Harmonic impedance: ~ 150 seconds per freq. and select buses
calculate_zh = false
write_ihdmax = false
number_of_buses = length(data["bus"])
zh_buses = collect(1:number_of_buses)

# Add principle
# {"maximum efficiency", "absolute equality", "maximin", "Kalai-Smorodinsky bargaining"}
data["principle"] = "maximum efficiency"

for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-3-6:2008"
    bus["ref_angle"] = 0.0
    bus["vmax"] = 1.7
end

# Switch on all generators
for (g, gen) in data["gen"]
    gen["gen_status"] = 1
end

# solve an OPF to find a feasible starting point
result_opf = PMs.solve_opf(data, PMs.ACPPowerModel, solver_nlp)

# update voltage and power setpoints
for (g, gen) in data["gen"]
    gen["p"] = result_opf["solution"]["gen"]["$g"]["pg"]
    gen["q"] = result_opf["solution"]["gen"]["$g"]["qg"]
end

for (b, bus) in data["bus"]
    bus["vm"] = result_opf["solution"]["bus"]["$b"]["vm"]
    bus["va"] = result_opf["solution"]["bus"]["$b"]["va"]
end

# Convert data format to include transformers:

data["xfmr"] = Dict{String, Any}()
sc_ratio = 0.15

for (br, branch) in data["branch"]
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    # Branches between different voltage levels are transformers
    if data["bus"]["$f_bus"]["base_kv"] !== data["bus"]["$t_bus"]["base_kv"]
        data["xfmr"][br] = xfmr = Dict{String, Any}()
        xfmr["f_bus"] = branch["f_bus"]
        xfmr["t_bus"] = branch["t_bus"]
        xfmr["xsc"] =  branch["br_x"] * sc_ratio
        xfmr["gsh"] =  0.0
        xfmr["r1"] = branch["br_r"] / 2
        xfmr["r2"] = branch["br_r"] / 2
        # Determine vector group: if Vp >= 30 kV -> star, if Vs >= 30.0 -> star, otherwise delta
        Vp = max(data["bus"]["$f_bus"]["base_kv"], data["bus"]["$t_bus"]["base_kv"])
        Vs = max(data["bus"]["$f_bus"]["base_kv"], data["bus"]["$t_bus"]["base_kv"])
        if Vp >= 30.0
            vgp = "Y"
            gnd1 = 1
            re1 = 0
        end

        if Vs >= 30.0
            vgs = "y"
            shift = "0"
            gnd2 = 1
            re2 = 0
        else
            vgs = "d"
            shift = "11"
            gnd2 = 0
            re2 = 0
        end
        xfmr["vg"] = join([vgp, vgs, shift])
        xfmr["gnd1"] = gnd1
        xfmr["gnd2"] = gnd2
        xfmr["re1"] = re1
        xfmr["re2"] = re2
        xfmr["xe1"] = 0.0
        xfmr["xe2"] = 0.0
        xfmr["rateA"] = branch["rate_a"] * data["baseMVA"]

        delete!(data["branch"], br)
    end
end



##### HARMONIC IMPEDANCE CALCULATION #######
if calculate_zh == true
    # Step 1: Calculate Harmonic impedance fr generators
    Ssc     = 2.1e9
    Sbase   = 100e6
    XRr     = 20.0

    z = Sbase / Ssc
    r = z / sqrt(1 + XRr^2)
    x = sqrt(z^2 - r^2)
    for (ng, gen) in data["gen"]
        gen["rsc"] = r
        gen["xsc"] = x
    end
    #### Add bus indexes for the case of non sequential bus numbers
    idx = 1
    for b in collect(sort(parse.(Int, keys(data["bus"]))))
        data["bus"]["$b"]["hb_idx"] = idx
        global idx = idx + 1
    end

    ##### Sslect frequency range and vector of busses for calculating impedance
    Hf = 50:50.0:2500.0

    # Calculate impedance
    @time Zh = calculate_pos_seq_harmonic_impedance(data, collect(Hf), zh_buses)

    filename = joinpath(HPM.BASE_DIR, "results", join([case, "_Zh.json"]))
    json_string = JSON.json(Zh)
    open(filename,"w") do f
    write(f, json_string)
    end
end

###################

# define the set of considered harmonics
H = [h for h in 1:30]


# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)

# first solve hOPF to have good starting values:
# hdata_opf = deepcopy(hdata_soc)

# for (n, nw) in hdata_opf["nw"]
#     if n !== "1" && n !== "2"
#         delete!(hdata_opf["nw"], n)
#     end
# end

# result_opf = HPM.solve_hopf(hdata_opf, PMs.IVRPowerModel, solver_nlp)




########
for (b, branch) in hdata_soc["nw"]["1"]["branch"]
    branch["c_rating"] = branch["c_rating"] * 10
end

for (x, xfmr) in hdata_soc["nw"]["1"]["xfmr"]
    xfmr["c_rating"] = xfmr["c_rating"] * 10
end



if write_ihdmax == true
    filename = joinpath(HPM.BASE_DIR, "results", join([case, "_Zh.json"]))
    Zh = Dict{String, Any}()
    open(filename) do f
    dicttxt = read(f, String)
    global Zh = JSON.parse(dicttxt)
    end

    for (nw, network) in hdata_soc["nw"]
        if network !== "nw"
            h = parse(Int, nw)
            for (l, load) in network["load"]
                load_bus = load["load_bus"]
                ihdmax = network["bus"]["$load_bus"]["ihdmax"] 
                hb_idx = network["bus"]["$load_bus"]["hb_idx"]
                load["cmdmax"] = ihdmax / sqrt(Zh["$hb_idx"][h]["re"]^2 + Zh["$hb_idx"][h]["im"]^2)
            end
        end
    end
    hdata_soc["skip KS precalc"] = true
end

s = Dict("fix_refbus_angle" => true)
total_time = @elapsed results_hhc = HPM.solve_hhc(hdata_soc, dHHC_SOC, solver_soc, solver_nlp; setting = s)

h_1 = H[1]
h_end = H[end]
harmonic_range = join(["_", "$h_1", "_", "$h_end", "_"])
filename = joinpath(HPM.BASE_DIR, "results", join([case, "_", harmonic_range, data["principle"],".json"]))
json_string = JSON.json(results_hhc)
open(filename,"w") do f
write(f, json_string)
end

# for (n, nw) in hdata_soc["nw"]
#     for (x, xfmr) in nw["xfmr"]
#         xfmr["c_rating"] = xfmr["c_rating"] * 100
#     end
# end

# hpf_data = deepcopy(hdata_soc)
# for n in keys(hpf_data["nw"])
#     if n ≠ "1"
#         delete!(hpf_data["nw"], n)
#     end
# end

# # solve hpf problem for the fundamental harmonic only
# hpf_results = HPM.solve_hpf(hpf_data,  PMs.IVRPowerModel, solver_nlp)
# hpf_results["solution"]["nw"]["1"]["branch"]["730"]
# hpf_results["solution"]["nw"]["1"]["xfmr"]["2498"]

#     for (b, bus) in hpf_results["solution"]["nw"]["1"]["bus"]
#         if bus["vm"] > 1.1
#             println( b, " " , bus["vm"])
#         end
#     end

