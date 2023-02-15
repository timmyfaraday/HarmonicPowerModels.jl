################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
"""
Example considering harmonic power flow for a two-bus example network taken 
from:
>   Harmonic Optimal Power Flow with Transformer Excitation by F. Geth and T. 
    Van Acker, pg. 7, § IV.A.
"""

# load pkgs
using HarmonicPowerModels
using Ipopt, JuMP, PowerModels

# set the solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

# read-in data
path = joinpath(HarmonicPowerModels.BASE_DIR,"test/data/matpower/two_bus_example_hpf.m")
data = PowerModels.parse_file(path)

# set the formulation
form = PowerModels.IVRPowerModel

# fundamental power flow
results_fund = PowerModels.solve_pf_iv(data, form, solver)
println("Results for the fundamental power flow:")
print_summary(results_fund["solution"])

# build the harmonic data
hdata = HarmonicPowerModels.replicate(data)

# harmonic power flow
results_harm = HarmonicPowerModels.solve_hpf(hdata, form, solver)
println("Results for the harmonic power flow")
println("Fundamental harmonic:")
print_summary(results_harm["solution"]["nw"]["1"])
println("Third harmonic:")
print_summary(results_harm["solution"]["nw"]["3"])