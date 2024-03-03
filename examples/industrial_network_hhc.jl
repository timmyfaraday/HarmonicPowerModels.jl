################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

# using pkgs
using HarmonicPowerModels, PowerModels
using Clarabel
using Ipopt 

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_soc = Clarabel.Optimizer
solver_nlp = Ipopt.Optimizer

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
data = PMs.parse_file(path)

# define the set of considered harmonics
H = [1, 3, 5, 7, 9, 13]

# solve HHC problem -- NLP
hdata_nlp = HPM.replicate(data, H=H)
results_hhc_nlp = HPM.solve_hhc(hdata_nlp, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)
results_hhc_soc = HPM.solve_hhc(hdata_soc, dHHC_SOC, solver_soc, solver_nlp)