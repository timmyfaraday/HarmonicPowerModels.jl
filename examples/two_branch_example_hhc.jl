################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# using pkgs
using HarmonicPowerModels, JuMP, Plots, PowerModels, Revise
# using Gurobi
using Ipopt 

# pkgs cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_ipopt    = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
# solver_gurobi   = JuMP.optimizer_with_attributes(Gurobi.Optimizer)

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/two_branch_example_hhc.m")
data = PMs.parse_file(path)

# define the set of considered harmonics
H=[1, 7] #[1, 3, 5, 7, 9, 13]

# solve HHC problem -- NLP
hdata_nlp = HPM.replicate(data, H=H)
results_hhc_nlp = HPM.solve_hhc(hdata_nlp, NLP_DHHC, solver_ipopt)

# solve HHC problem -- QC 
hdata_qc = HPM.replicate(data, H=setdiff(H,1))
results_hhc_qc = HPM.solve_hhc(hdata_qc, QC_DHHC, solver_ipopt)

# # solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=setdiff(H,1))
results_hhc_soc = HPM.solve_hhc(hdata_soc, SOC_DHHC, solver_gurobi)