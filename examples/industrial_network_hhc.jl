################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# using pkgs
using HarmonicPowerModels, Ipopt, JuMP, PowerModels

# pkgs cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver = Ipopt.Optimizer

# set the formulation
form = PowerModels.IVRPowerModel

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
data = PMs.parse_file(path)

# build harmonic data
hdata = HPM.replicate(data, H=[1, 3, 5, 7, 9, 13])

# add the individual harmonic distortion limits                                 
ihdmax = Dict("1" => 1.10, "3" => 0.05, "5" => 0.06, "7" => 0.05, "9" => 0.015, "13" => 0.03)
for (nw,ntw) in hdata["nw"], (nb,bus) in ntw["bus"]
        bus["ihdmax"] = ihdmax[nw]
end 

# harmonic hosting capacity
results_hhc = HPM.solve_hhc(hdata, form, solver)

# print the results
println("Fundamental harmonic:")
print_summary(results_hhc["solution"]["nw"]["1"])
println("Third harmonic:")
print_summary(results_hhc["solution"]["nw"]["3"])
println("Fifth harmonic:")
print_summary(results_hhc["solution"]["nw"]["5"])
println("Seventh harmonic:")
print_summary(results_hhc["solution"]["nw"]["7"])
println("Nineth harmonic:")
print_summary(results_hhc["solution"]["nw"]["9"])
println("Thirteen harmonic:")
print_summary(results_hhc["solution"]["nw"]["13"])
