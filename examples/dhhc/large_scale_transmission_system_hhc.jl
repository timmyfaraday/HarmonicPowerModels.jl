################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Example considering the harmonic hosting capacity of an industrial power     #
# system taken from: Exact Lower Bound on Equitable Harmonic Hosting Capacity  #
# by T. Van Acker and H. Ergun, pg. 8, § III.B.                                #
# ---------------------------------------------------------------------------- #
# This numerical illustration aims to demonstrate the developed harmonic       #
# hosting capacity optimization model on a large transmission grid. The Rte    #
# 1888 bus test system, inspired by the French transmission grid is used for   #
# that purpose. The test system has twelve voltage levels ranging from 380 kV  #
# to 3 kV, 551 transformers, 1980 branches, and 1000 units. All necessary data #
# are derived from the MATPOWER source file [24], and the power quality        #
# standard IEC61000-3-6:2008 [10]. As the source file does not specify the     #
# transformer vector group, a YNYn0 vector group is assumed for transformers   #
# with a primary voltage greater than or equal to 30 kV, while a YNd11 vector  #
# group is assumed for all others. The considered harmonic set is              #
# H ∈ {2, 3, ..., 50}. The aim of this numerical illustration is to            #
# 1) show the impact of the chosen fairness principle on a power system of     #
#    realistic size, and                                                       #
# 2) show the computational efficiency of the second-order cone model.         #
# Only the second-order cone formulation of the harmonic hosting capacity      #
# model is considered, and is solved using GUROBI V11.0.3.                     #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# using pkgs
using HarmonicPowerModels, PowerModels
using Gurobi
using Ipopt 
using JSON
using Plots
using StatsPlots
using SparseArrays
using JuMP

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels 

# Include function to calculate Harmonic impedance
include(joinpath(HPM.BASE_DIR,"src/util/hi.jl"))

# set the solver
solver_soc = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "BarConvTol" => 1e-4)
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
data["principle"] = "maximin"

# Enforce bus voltage limits and standard
for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-3-6:2008"
    bus["ref_angle"] = 0.0
    bus["vmax"] = 1.1#min(1.1, bus["vmax"])
    bus["vmin"] = 0.9##max(0.90, bus["vmin"]) 
end

# get rid of negative reactances and resistances
for (br, branch) in data["branch"]
    branch["tap"] = 1.0
    branch["shift"] = 0.0
    branch["c_rating"] = (branch["rate_a"] / 0.9) 
end

# Switch on all generators and set their lower bound to zero
for (g, gen) in data["gen"]
    gen[ "gen_status"] = 1
    gen["pmin"] = 0.0
    gen_bus = gen["gen_bus"]
    data["bus"]["$gen_bus"]["bus_type"] = 2
end
# solve an OPF to find a feasible starting point
result_opf = PMs.solve_opf_iv(data, PMs.IVRPowerModel, solver_nlp)

# update voltage and power setpoints
for (g, gen) in data["gen"]
    gen["p"] = result_opf["solution"]["gen"]["$g"]["pg"]
    gen["q"] = result_opf["solution"]["gen"]["$g"]["qg"]
    gen["crg"] = result_opf["solution"]["gen"]["$g"]["crg"]
    gen["cig"] = result_opf["solution"]["gen"]["$g"]["cig"]
    gen["cm"] = sqrt(gen["crg"]^2 + gen["cig"]^2)
    gen_bus = gen["gen_bus"]
    v = data["bus"]["$gen_bus"]["base_kv"]
    xd_ohm = (v * 1000)^2 / (sqrt(gen["pmax"]^2 + gen["qmax"]^2) * data["baseMVA"] * 1e6)
    zbase = (v * 1000)^2 / (data["baseMVA"] * 1e6)
    xd_pu = xd_ohm / zbase
    gen["bg"] = - 1 / xd_pu   # x = Sbase / Snom = 1 / Snom in pu. Empiric value Hapold & Oeding Table A.4
    gen["c_rating"]= gen["cm"] # amske sure current rating is fine
end

for (br, branch) in data["branch"]
    branch["cr_fr"] = result_opf["solution"]["branch"]["$br"]["cr_fr"]
    branch["ci_fr"] = result_opf["solution"]["branch"]["$br"]["ci_fr"]
    branch["cr_to"] = result_opf["solution"]["branch"]["$br"]["cr_to"]
    branch["ci_to"] = result_opf["solution"]["branch"]["$br"]["ci_to"] 
    branch["cm_fr"] = sqrt(branch["cr_fr"]^2 + branch["ci_fr"]^2)
    branch["cm_to"] = sqrt(branch["cr_to"]^2 + branch["ci_to"]^2)
end

for (b, bus) in data["bus"]
    bus["vm"] = sqrt(result_opf["solution"]["bus"]["$b"]["vr"]^2 + result_opf["solution"]["bus"]["$b"]["vi"]^2)
    bus["va"] = atan(result_opf["solution"]["bus"]["$b"]["vi"] / result_opf["solution"]["bus"]["$b"]["vr"])

    bus["vi"] = result_opf["solution"]["bus"]["$b"]["vi"]
    bus["vr"] = result_opf["solution"]["bus"]["$b"]["vr"]
    bus["vmax"] = bus["vmax"] + 1e-5
end


# Convert data format to include transformers:

data["xfmr"] = Dict{String, Any}()
sc_ratio = 1.0

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
        xfmr["ctm_fr"] = branch["cm_fr"]
        xfmr["ctm_to"] = branch["cm_to"]
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
        xfmr["rateA"] = branch["rate_a"] #* data["baseMVA"] # replicate function devides by the base MVA....
        xfmr["c_rating"] = branch["c_rating"]

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

idx = 1
for b in collect(sort(parse.(Int, keys(data["bus"]))))
    data["bus"]["$b"]["hb_idx"] = idx
    global idx = idx + 1
end

# define the set of considered harmonics
H = [h for h in 1:50]

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)


#######
# for (b, branch) in hdata_soc["nw"]["1"]["branch"]
#     branch["c_rating"] = branch["c_rating"] * 10 # this needs to be checked.......
# end


# for (x, xfmr) in hdata_soc["nw"]["1"]["xfmr"]
#     xfmr["c_rating"] = xfmr["c_rating"] * 10
# end

if write_ihdmax == true
    #### Add bus indexes for the case of non sequential bus numbers

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

# h_1 = H[1]
# h_end = H[end]
# harmonic_range = join(["_", "$h_1", "_", "$h_end", "_"])
# filename = joinpath(HPM.BASE_DIR, "results", join([case, "_", harmonic_range, data["principle"],".json"]))
# json_string = JSON.json(results_hhc)
# open(filename,"w") do f
# write(f, json_string)
# end


# zh127 = [sqrt(Zh["127"][h]["re"]^2 + Zh["127"][h]["im"]^2)  for h in 1:50]
# index = 1:50
# mh = hcat(index, zh127)
# header = ["h" "zh"]
# writedlm(joinpath(HPM.BASE_DIR,"results", join(["zh127.csv"])),  [header ; mh], ',')


# plot([bus["vm"] for (i, bus) in hdata_soc["nw"]["1"]["bus"]], label = "voltage magnitude")

# hpf_data = deepcopy(hdata_soc)
# for n in keys(hpf_data["nw"])
#     if n ≠ "1"
#         delete!(hpf_data["nw"], n)
#     end
# end

# # solve hpf problem for the fundamental harmonic only
# hpf_results = HPM.solve_hpf(hpf_data, dHHC_NLP, solver_nlp)

# for (b, bus) in hpf_results["solution"]["nw"]["1"]["bus"]
#     if bus["vm"] > 1.15
#         println( b, " " , bus["vm"])
#     end
# end

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

# first solve hOPF to have good starting values:
# hdata_opf = deepcopy(hdata_soc)

# for (n, nw) in hdata_opf["nw"]
#     if n !== "1" && n !== "2"
#         delete!(hdata_opf["nw"], n)
#     end
# end

# result_opf = HPM.solve_hopf(hdata_opf, PMs.IVRPowerModel, solver_nlp)