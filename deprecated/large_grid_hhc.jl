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
# case = "nem2300harmonic"
# extension = ".json"

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
    bus["standard"] = "AS/NZS61000-3-6"
    bus["vmax"] = 1.3
    bus["ref_angle"] = 0.0
    # if b == "2136"
    #     bus["bus_type"] = 2
    # end
    # if b == "79" || b == "80"
    #     bus["vmax"] = 1.8
    # end
end

# delete!(data["bus"],"79")
# delete!(data["bus"],"80")
# delete!(data["xfmr"],"332")
# delete!(data["xfmr"],"333")
# delete!(data["gen"],"41")
# delete!(data["gen"],"42")

# for (g,gen) in data["gen"]
#     if !haskey(gen, "gsc") 
#         gen["gsc"] = 0.0
#     end
#     if !haskey(gen, "bsc") 
#         gen["bsc"] = 0.0
#     end
# end

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

idx = 1
for b in collect(sort(parse.(Int, keys(data["bus"]))))
    data["bus"]["$b"]["hb_idx"] = idx
    global idx = idx + 1
end


# number_of_buses = length(data["bus"])
# Hf = 50:50.0:100.0

# # Calculate impedance
# @time Zh = calculate_pos_seq_harmonic_impedance(data, collect(Hf), collect(1:2))

# # Zh_inv = deepcopy(Zh)

# Y = build_admittance_matrix(data)
# h = 50.0
# Yh = eval.(Y)

# Yhc = zeros(ComplexF64, 1888, 1888)
# for i in 1:1888
#     for j in 1:1888
#         Yhc[i, j] = float(Yh[i,j])
#     end
# end

# Yhs = SparseArrays.sparse(Yhc)



# I = [0.0 + 0.0im for i in 1:1888]
# I[1] = 1.0 + 0.0im

# Zh = (Yhs \ I)

# define the set of considered harmonics
H = [h for h in 1:50]#[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49]

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)
for (n, nw) in hdata_soc["nw"]
    for (x, xfmr) in nw["xfmr"]
        xfmr["c_rating"] = xfmr["c_rating"] * 100
    end
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



