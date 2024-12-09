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


# Add principle
# {"maximum efficiency", "absolute equality", "maximin", "Kalai-Smorodinsky bargaining"}
data["principle"] = "Kalai-Smorodinsky bargaining"
data["jeremy"] = true

for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-3-6:2008"
    bus["ref_angle"] = 0.0
end

# Switch on all generators
for (g, gen) in data["gen"]
    gen["gen_status"] = 1
end

# Convert data format to include transformers:

data["xfmr"] = Dict{String, Any}()

for (br, branch) in data["branch"]
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    # Branches between different voltage levels are transformers
    if data["bus"]["$f_bus"]["base_kv"] !== data["bus"]["$t_bus"]["base_kv"]
        data["xfmr"][br] = xfmr = Dict{String, Any}()
        xfmr["f_bus"] = branch["f_bus"]
        xfmr["t_bus"] = branch["t_bus"]
        xfmr["xsc"] =  branch["br_x"]
        xfmr["gsh"] =  0.0
        xfmr["r1"] = branch["br_r"] / 2
        xfmr["r2"] = branch["br_r"] / 2
        # Determine vector group: if Vp >= 30 kV -> star, if Vs >= 30.0 -> star, otherwise delta
        Vp = max(data["bus"]["$f_bus"]["base_kv"], data["bus"]["$t_bus"]["base_kv"])
        Vs = max(data["bus"]["$f_bus"]["base_kv"], data["bus"]["$t_bus"]["base_kv"])
        if Vp >= 30.0
            vgp = "Y"
            gnd1 = 1
            re1 = 1e-5
        else
            vgp = "D"
            gnd1 = 0
            re1 = 0
        end

        if Vs >= 30.0
            vgs = "y"
            shift = "0"
            gnd2 = 1
            re2 = 1e-5
        else
            vgs = "d"
            shift = "3"
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
        xfmr["rateA"] = branch["rate_a"]

        delete!(data["branch"], br)
    end
end



##### HARMONIC IMPEDANCE CALCULATION ########
# Step 1: Calculate Harmonic impedance fr generators
# Ssc     = 2.1e9
# Sbase   = 100e6
# XRr     = 20.0

# z = Sbase / Ssc
# r = z / sqrt(1 + XRr^2)
# x = sqrt(z^2 - r^2)
# for (ng, gen) in data["gen"]
#     gen["rsc"] = r
#     gen["xsc"] = x
# end

#### Add bus indexes for the case of non sequential bus numbers
# idx = 1
# for b in collect(sort(parse.(Int, keys(data["bus"]))))
#     data["bus"]["$b"]["hb_idx"] = idx
#     global idx = idx + 1
# end

##### Sslect frequency range and vector of busses for calculating impedance
# number_of_buses = length(data["bus"])
# Hf = 50:50.0:100.0

# # Calculate impedance
# @time Zh = calculate_pos_seq_harmonic_impedance(data, collect(Hf), collect(1:number_of_buses))


###################

# define the set of considered harmonics
H = [h for h in 1:50]

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)
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



