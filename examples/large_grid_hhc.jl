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

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_soc = Gurobi.Optimizer
solver_nlp = Ipopt.Optimizer

# read-in data 
# case = "pglib_opf_case2848_rte.m"
# case = "case1803_snem.m"
case = "nem2300harmonic"
extension = ".json"
csv_filename = joinpath(HPM.BASE_DIR,"results/voltage_phasors.csv")

if extension == ".m"
    path = joinpath(HPM.BASE_DIR,"test","data","matpower", join([case, extension]))
elseif extension == ".json"
    path = joinpath(HPM.BASE_DIR,"test","data","json", join([case, extension]))
else
    Memento.warn(_PM._LOGGER, "This input file is not yet supported!")
end


data = PMs.parse_file(path)

# Add principle
data["principle"] = "maximin"

# define the set of considered harmonics
H = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49]
H = H[1:6]

data_tf = HPM.parse_transformers(data)
HPM.add_bus_harmonic_limits(data_tf) #standard = "AS/NZS61000-3-6"

# solve HHC problem -- NLP
# hdata_nlp = HPM.replicate(data, H=H)
# results_hhc_nlp = HPM.solve_hhc(hdata_nlp, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc = HPM.replicate(data, H=H)
total_time = @elapsed results_hhc_soc_1 = HPM.solve_hhc(hdata_soc, dHHC_SOC, solver_soc, solver_nlp)

# for n in [3,5,7,9,13], l in 1:4
#     ca_nlp = rad2deg(results_hhc_nlp["solution"]["nw"]["$n"]["load"]["$l"]["ca"])
#     ca_soc = rad2deg(results_hhc_soc["solution"]["nw"]["$n"]["load"]["$l"]["ca"])
#     println("h=$n, l=$l: $(round(ca_nlp,digits=8)) vs $(round(ca_soc,digits=8))")
# end

# HPM.csv_export(results_hhc_nlp, csv_filename)