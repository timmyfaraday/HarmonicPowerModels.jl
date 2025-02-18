################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

"""
Example considering harmonic power flow for a two-bus example network taken 
from:
>   Harmonic Optimal Power Flow with Transformer Excitation by F. Geth and T. 
    Van Acker, pg. 7, § IV.A.
"""

# load pkgs
using HarmonicPowerModels, PowerModels
using Ipopt

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver = Ipopt.Optimizer

# read-in data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/two_bus_example_hpf.m")
data = PMs.parse_file(path)

# solve PF problem
results_fund = PMs.solve_pf_iv(data, PMs.IVRPowerModel, solver)

# define the set of considered harmonics
H = [1, 3, 5, 7, 9, 13]

# solve HPF problem
hdata = HPM.replicate(data, H=H)
results_harm = HarmonicPowerModels.solve_hpf(hdata, PMs.IVRPowerModel, solver)

# results
println("Results for the fundamental power flow:")
print_summary(results_fund["solution"])

println("Results for the harmonic power flow")
println("Fundamental harmonic:")
print_summary(results_harm["solution"]["nw"]["1"])
println("Third harmonic:")
print_summary(results_harm["solution"]["nw"]["3"])
