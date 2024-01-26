################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# using pkgs
using HarmonicPowerModels, JuMP, Plots, PowerModels, Revise
using SCS
using Ipopt 

# pkgs cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_soc = JuMP.optimizer_with_attributes(SCS.Optimizer)
solver_nlp = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
data = PMs.parse_file(path)

# define the set of considered harmonics
H = [1, 3, 5, 7, 9, 13]

# solve HHC problem -- NLP
hdata_nlp = HPM.replicate(data, H=H)
results_hhc_nlp = HPM.solve_hhc(hdata_nlp, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=setdiff(H,1))
results_hhc_soc = HPM.solve_hhc(hdata_soc, dHHC_SOC, solver_soc)

for (n, nw) in results_hhc_soc["solution"]["nw"]
    print("Harmonic Order ,",n, ":","\n")
    for (b, bus) in results_hhc_nlp["solution"]["nw"][n]["bus"]
        vmsoc = results_hhc_soc["solution"]["nw"][n]["bus"][b]["vm"]
        print("Bus ", b,": vm nlp = ", bus["vm"], ", vm soc = ",vmsoc, "\n")
    end
end