

# # small-case sanity tests ######################################################
# ## single branch (r = 0.01, x=0.1, b / 2 = 0.01)
# f = 750.0

# r = 0.01
# x = 0.1
# b = 0.01

# ### method
# data = Dict{String,Any}("bus"       => Dict{String,Any}("$n" => [] for n in 1:2),
#                         "branch"    => Dict{String,Any}("1" => Dict{String,Any}("f_bus" => 1,
#                                                                                 "t_bus" => 2,
#                                                                                 "br_r"  => r,
#                                                                                 "br_x"  => x,
#                                                                                 "b_fr"  => b,
#                                                                                 "b_to"  => b)),
#                         "gen"       => Dict{String,Any}(),
#                         "xfmr"      => Dict{String,Any}())
# Zc = calculate_pos_seq_harmonic_impedance(data, [f], [1,2])
# ### hand-calc
# R = r
# L = x / 2 / pi / 50
# C = b / 2 / pi / 50

# Zl = R * sqrt(f / 50) + im * 2 * pi * f * L + 1 / (im * 2 * pi * f * C)
# Zr = 1 / (im * 2 * pi * f * C)

# Zm = 1 / (1 / Zl + 1 / Zr)

# ### tests 
# @assert Zc[1][1] ≈ Zm
# @assert Zc[1][1] ≈ Zc[2][1]

# ## single generator (Ssc=2.1GVA, XR=20) #######################################
# f = 425.0

# Ssc     = 2.1e9
# Sbase   = 100e6
# XRr     = 20.0

# z = Sbase / Ssc
# r = z / sqrt(1 + XRr^2)
# x = sqrt(z^2 - r^2)

# @assert x == z / sqrt(1 + (1/XRr)^2)

# ### method
# data = Dict{String,Any}("bus"       => Dict{String,Any}("$n" => [] for n in 1),
#                         "branch"    => Dict{String,Any}(),
#                         "gen"       => Dict{String,Any}("1" => Dict{String,Any}("gen_bus" => 1,
#                                                                                 "rsc" => r,
#                                                                                 "xsc" => x)),
#                         "xfmr"      => Dict{String,Any}())
# Zc = calculate_pos_seq_harmonic_impedance(data, [f], [1])

# ### hand-calc
# R = r
# L = x / 2 / pi / 50

# Zm = R * sqrt(f / 50) + im * 2 * pi * f * L

# ### tests 
# @assert Zc[1][1] ≈ Zm

# ## single transformer ##########################################################
# f = 695.0

# x = 0.208
# g = 0.00001
# r = 2 * 0.004988662

# ### method
# data = Dict{String,Any}("bus"       => Dict{String,Any}("$n" => [] for n in 1:2),
#                         "branch"    => Dict{String,Any}(),
#                         "gen"       => Dict{String,Any}(),
#                         "xfmr"      => Dict{String,Any}("1" => Dict{String,Any}("f_bus" => 1,
#                                                                                 "t_bus" => 2,
#                                                                                 "r1" => r / 2,
#                                                                                 "r2" => r / 2,
#                                                                                 "xsc" => x,
#                                                                                 "gsh" => g)))
# Zc = calculate_pos_seq_harmonic_impedance(data, [f], [1,2])

# ### hand-calc, from node 1
# R = r + 1 / g
# L = x / 2 / pi / 50

# Zm = R * sqrt(f / 50.0) + im * 2 * pi * f * L 

# ### tests 
# @assert Zc[1][1] ≈ Zm

## 30bus System ################################################################
# using HarmonicPowerModels
# using PowerModels
# using Plots

# const PMs = PowerModels
# const HPM = HarmonicPowerModels

# path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")
# data = PMs.parse_file(path)

# bus_id = 24

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

# @time Zc = calculate_pos_seq_harmonic_impedance(data, collect(50.0:10.0:2500.0), collect(1:30))

# # plot 
# plot(50.0:10.0:2500.0, abs.(Zc[bus_id]), yaxis=:log)

# ToDo's

# (1) start inductive only on HV - Done
# (2) One by one add B to HV cables - Done
# (3) Have a figure where the impedance is shown fpor each additional cable for one HV bus - Done
# (4) Calculate HHC for all cases and show that with increaseing HV cables drops: absolute equality - Done
# (5) Standard method: How much can we inject at each node: Only at buses with loads, 
# (6) Sensitivity of the bugdet: Explain in the results that the budget would need to be set to ~50% for this case,
#     but in general this might not be applicable everywhere -> thus importance of the HHC definition as in this paper!



#######  Analysis for PES GM paper  #####

using HarmonicPowerModels
using PowerModels
using Plots

const PMs = PowerModels
const HPM = HarmonicPowerModels

## Include function to calculate Harmonic impedance
include(joinpath(HPM.BASE_DIR,"src/util/hi.jl"))
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")
data = PMs.parse_file(path)
data_cable = deepcopy(data)

number_of_buses = 30
bus_id = 8
cable_branches = [0 1 2 3 4 5 6 7 8 9 10 34 35]
H = 50:50.0:2500

# First set capacitances to zero
for c in cable_branches
    if c !== 0
        data_cable["branch"]["$c"]["br_r"] = data_cable["branch"]["$c"]["br_r"]
        data_cable["branch"]["$c"]["br_x"] = data_cable["branch"]["$c"]["br_x"]
        data_cable["branch"]["$c"]["b_fr"] = 0.0
        data_cable["branch"]["$c"]["b_to"] = 0.0
    end
end

Ssc     = 2.1e9
Sbase   = 100e6
XRr     = 20.0

z = Sbase / Ssc
r = z / sqrt(1 + XRr^2)
x = sqrt(z^2 - r^2)
for (ng, gen) in data_cable["gen"]
    gen["rsc"] = r
    gen["xsc"] = x
end

Zc = Dict{String, Any}(["$i" => nothing for i in cable_branches])

Zc["0"] = calculate_pos_seq_harmonic_impedance(data_cable, collect(H), collect(1:number_of_buses))

for cb in cable_branches 
    if cb !== 0
        data_cable["branch"]["$cb"]["br_r"] = data["branch"]["$cb"]["br_r"]
        data_cable["branch"]["$cb"]["br_x"] = data["branch"]["$cb"]["br_x"]
        data_cable["branch"]["$cb"]["b_fr"] = data["branch"]["$cb"]["b_fr"]
        data_cable["branch"]["$cb"]["b_to"] = data["branch"]["$cb"]["b_to"]

        println("Determining harmonic impedances for adding cable ", "$cb")

        @time Zc["$cb"] = calculate_pos_seq_harmonic_impedance(data_cable, collect(H), collect(1:number_of_buses))
    end
end

# plot impedances
hi = plot(H, abs.(Zc["0"][bus_id]), yaxis=:log, label = "No HV cable", fontfamily = "Computer Modern", xlabel = "f in Hz", ylabel = "\$|Z(f)|\$ in pu", legend=:outertopright)
for cb in cable_branches
    if cb !== 0
        plot!(hi, H, abs.(Zc["$cb"][bus_id]), yaxis=:log, label = join(["With branch ", "$cb", " as HV cable"]), fontfamily = "Computer Modern", xlabel = "f in Hz", ylabel = "\$|Z(f)|\$ in pu",legend=:outertopright)

    end
end
display(hi)

savefig(hi, joinpath(HPM.BASE_DIR,"results/impedance_scan.pdf"))

## Calculation of the HHC

# using pkgs
using Ipopt 
using Revise
using PrettyTables
using Gurobi

include(joinpath(HPM.BASE_DIR,"src/util/hi.jl"))
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")
data = PMs.parse_file(path)
data_cable = deepcopy(data)

# set the solver
solver_nlp = Ipopt.Optimizer

# check to improve with OPF
data_cable["branch"]["22"]["rate_a"] = 0.25 * sqrt(3) / 0.9
data_cable["branch"]["13"]["rate_a"] = 0.25 * sqrt(3) / 0.9

# First get rid of cables, e.g., shunts.
for c in cable_branches
    if c !== 0
        data_cable["branch"]["$c"]["br_r"] = data_cable["branch"]["$c"]["br_r"]
        data_cable["branch"]["$c"]["br_x"] = data_cable["branch"]["$c"]["br_x"]
        data_cable["branch"]["$c"]["b_fr"] = 0.0
        data_cable["branch"]["$c"]["b_to"] = 0.0
    end
end

# define the set of considered harmonics
H = [i for i in 1:50]

# RESULTS ######################################################################
## absolute equality (ae) ######################################################
data_cable["principle"] = "absolute equality"

# solve HHC problem -- SOC
hdata = HPM.replicate(data_cable, H=H)

result_hhc = Dict{String, Any}(["$i" => nothing for i in 0:maximum(cable_branches)])

result_hhc["0"] = HPM.solve_hhc(hdata, dHHC_SOC, solver_nlp, solver_nlp)

for cb in cable_branches 
    if cb !== 0
        data_cable["branch"]["$cb"]["br_r"] = data["branch"]["$cb"]["br_r"]
        data_cable["branch"]["$cb"]["br_x"] = data["branch"]["$cb"]["br_x"]
        data_cable["branch"]["$cb"]["b_fr"] = data["branch"]["$cb"]["b_fr"]
        data_cable["branch"]["$cb"]["b_to"] = data["branch"]["$cb"]["b_to"]

        hdata = HPM.replicate(data_cable, H=H)

        println("Determining HHC for adding cable ", "$cb")

        result_hhc["$cb"] = HPM.solve_hhc(hdata, dHHC_SOC, solver_nlp, solver_nlp)
    end
end

hhc_cable = []
for cb in cable_branches
    push!(hhc_cable, result_hhc["$cb"]["objective"])
end
hhc = plot(Int.(0:length(cable_branches)-1), collect(hhc_cable), xlabel = "Number of HV cables", ylabel = "HHC in pu", fontfamily = "Computer Modern", label = "Hosting capacity optimisation")


###### Calculate harmonic currents based on the standard:
# -> each harmonic voltage needs to be less than 3% of fundamental
# -> THD below 5%

H = 50:50.0:2500
ih = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])
ih_min = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])
ih_lin = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])
ih_min_lin = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])

load_buses = [data["load"]["$l"]["load_bus"] for l in sort(parse.(Int, keys(data["load"])))]
global_contribution = [abs.(data["load"]["$l"]["pd"] + im * data["load"]["$l"]["qd"]) for l in sort(parse.(Int, keys(data["load"])))] ./ sum([abs.(load["pd"] + im * load["qd"]) for (l, load) in data["load"]])


ihd_limits = [1.00000, 0.02000, 0.05000, 0.01000, 0.06000, 
0.00500, 0.05000, 0.00500, 0.01500, 0.00500, 
0.03500, 0.00458, 0.03000, 0.00429, 0.00400, 
0.00406, 0.02000, 0.00389, 0.01761, 0.00375, 
0.00200, 0.00364, 0.01408, 0.00354, 0.01275, 
0.00346, 0.00200, 0.00339, 0.01061, 0.00333, 
0.00975, 0.00328, 0.00200, 0.00324, 0.00833, 
0.00319, 0.00773, 0.00316, 0.00200, 0.00313, 
0.00671, 0.00310, 0.00627, 0.00307, 0.00200, 
0.00304, 0.00551, 0.00302, 0.00518, 0.00300]


# General summation law:
# I = ∑i  J[i,h]^α)^(1/α):
# h < 5 -> α = 1 
# 5 <= h < 10 -> α = 1.4
# 5 <= h < 10 -> α = 2

# Calculation with actual harmonic impedance
for c in cable_branches
    Zh = Zc["$c"]
    for h in 2:length(H)
        if h < 5
            α = 1
        elseif 5 <= h && h <= 10
            α = 1.4
        else
            α = 2
        end
        α = 1
        ih["$c"][h] = sum([abs(ihd_limits[h] / Zh[load_buses[i]][h])^α  * global_contribution[i] for i in 1:length(load_buses)])^(1/α)
        ih_min["$c"][h] = minimum([abs(ihd_limits[h] / Zh[load_buses[i]][h]) * global_contribution[i]  for i in 1:length(load_buses)]')
    end
end

# Calculation with linear impdedance assumption harmonic impedance
for c in cable_branches
    Zh = Zc["$c"]
    for h in 2:length(H)
        if h < 5
            α = 1
        elseif 5 <= h && h <= 10
            α = 1.4
        else
            α = 2
        end
        α = 1
        ih_lin["$c"][h] = sum([abs(ihd_limits[h] / (h * Zh[load_buses[i]][1]))^α  * global_contribution[i] for i in 1:length(load_buses)])^(1/α)
        ih_min_lin["$c"][h] = minimum([abs(ihd_limits[h] / (h * Zh[load_buses[i]][1])) * global_contribution[i] for i in 1:length(load_buses)]')
    end
end

thd_ih = [sqrt(sum(ih["$c"][2:end].^2)) for c in cable_branches]
sum_ih_min = [sum(ih_min["$c"][2:end]) for c in cable_branches] .* length(load_buses)
thd_ih_lin = [sqrt(sum(ih_lin["$c"][2:end].^2)) for c in cable_branches]
sum_ih_lin = [sum(ih_lin["$c"][2:end]) for c in cable_branches]
sum_ih_min_lin = [sum(ih_min_lin["$c"][2:end]) for c in cable_branches] .* length(load_buses)


plot!(hhc, Int.(0:length(cable_branches)-1), thd_ih', xlabel = "Number of HV cables", ylabel = "HHC in pu", fontfamily = "Computer Modern", label = "IEC61000-3-6:2008", ylim = (0,1))


savefig(hhc, joinpath(HPM.BASE_DIR,"results/hosting_capacity.pdf"))



