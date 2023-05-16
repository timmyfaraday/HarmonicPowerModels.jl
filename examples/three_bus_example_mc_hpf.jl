################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
"""
Example considering harmonic power flow for a three-bus example unbalanced 
network taken from:
>   PowerModelsDistribution.jl
"""

# load pkgs
using HarmonicPowerModels
using Ipopt, JuMP, PowerModelsDistribution

# set the solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

# read-in data
path = joinpath(HarmonicPowerModels.BASE_DIR,"test/data/opendss/three_bus_example_mc_hpf.dss")
data = PowerModelsDistribution.parse_file(path)
data_math = transform_data_model(data, kron_reduce=false, phase_project=false)

# set the formulation
form = PowerModelsDistribution.IVRENPowerModel

# fundamental power flow
# results_fund = PowerModelsDistribution.solve_mc_pf(data_math, form, solver)

# build the harmonic data
mn_data = make_multinetwork(data)
mn_data_math = transform_data_model(mn_data, kron_reduce=false, phase_project=false)

# harmonic power flow
results_fund = HarmonicPowerModels.solve_mc_hpf(mn_data_math, form, solver)