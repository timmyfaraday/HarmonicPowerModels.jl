################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# This numerical illustration considers a section of a real-world industrial   #
# phase-balanced three-phase power system across five voltage levels, with     #
# four transformers, three cables, and four harmonic units. On top of the      #
# fundamental harmonic for the non-linear formulation, the considered harmonic #
# set is 𝓗 ∈ {3,5,7}, i.e., one harmonic from each sequence subset.           #
# The aim of this numerical illustration is to                                 #
#   - validate the proposed formulations, e.g., zero sequence blocking;        #
#   - show the equivalence between the non-linear and second-order cone        #
#     models; and                                                              #
#   - show the impact of the chosen fairness principle.                        #
# For the purpose of equivalency, both models are solved using Ipopt, where    #
# the second-order constraints are translated to their equivalent non-convex   #
# quadratic form using the appropriate MathOptInterface constraint bridge.     #
# All results are presented in per-unit using a 100 MVA basis.                 #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

# using pkgs
using HarmonicPowerModels, PowerModels
using Ipopt 
using Revise

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_nlp = Ipopt.Optimizer

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
data = PMs.parse_file(path)

# define the set of considered harmonics
H = [1, 3, 5, 7]

## absolute equality (ae) ######################################################
data["principle"] = "absolute equality"

# solve HHC problem -- NLP
hdata_nlp_ae = HPM.replicate(data, H=H)
results_hhc_nlp_ae = HPM.solve_hhc(hdata_nlp_ae, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc_ae = HPM.replicate(data, H=H)
results_hhc_soc_ae = HPM.solve_hhc(hdata_soc_ae, dHHC_SOC, solver_nlp, solver_nlp)

## maximum efficiency (me) #####################################################
data["principle"] = "maximum efficiency"

# solve HHC problem -- NLP
hdata_nlp_me = HPM.replicate(data, H=H)
results_hhc_nlp_me = HPM.solve_hhc(hdata_nlp_me, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc_me = HPM.replicate(data, H=H)
results_hhc_soc_me = HPM.solve_hhc(hdata_soc_me, dHHC_SOC, solver_nlp, solver_nlp)

## maximin (mm) ################################################################
data["principle"] = "maximin"

# solve HHC problem -- NLP
hdata_nlp_mm = HPM.replicate(data, H=H)
results_hhc_nlp_mm = HPM.solve_hhc(hdata_nlp_mm, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc_mm = HPM.replicate(data, H=H)
results_hhc_soc_mm = HPM.solve_hhc(hdata_soc_mm, dHHC_SOC, solver_nlp, solver_nlp)

## Kalai-Smorodinsky bargaining (ks) ###########################################
data["principle"] = "Kalai-Smorodinsky bargaining"

# solve HHC problem -- NLP
hdata_nlp_ks = HPM.replicate(data, H=H)
results_hhc_nlp_ks = HPM.solve_hhc(hdata_nlp_ks, dHHC_NLP, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc_ks = HPM.replicate(data, H=H)
results_hhc_soc_ks = HPM.solve_hhc(hdata_soc_ks, dHHC_SOC, solver_nlp, solver_nlp)

## voltage
THD     = [sqrt(sum([results_hhc_nlp_mm["solution"]["nw"]["$nh"]["bus"]["$nb"]["vm"] for nh ∈ H if nh ≠ 1].^2) / 
                results_hhc_nlp_mm["solution"]["nw"]["1"]["bus"]["$nb"]["vm"]^2) for nb in 0:8]
IHD₃    = [sqrt(results_hhc_nlp_mm["solution"]["nw"]["3"]["bus"]["$nb"]["vm"]^2 / 
                results_hhc_nlp_mm["solution"]["nw"]["1"]["bus"]["$nb"]["vm"]^2) for nb in 0:8]
IHD₅    = [sqrt(results_hhc_nlp_mm["solution"]["nw"]["5"]["bus"]["$nb"]["vm"]^2 / 
                results_hhc_nlp_mm["solution"]["nw"]["1"]["bus"]["$nb"]["vm"]^2) for nb in 0:8]
IHD₇    = [sqrt(results_hhc_nlp_mm["solution"]["nw"]["7"]["bus"]["$nb"]["vm"]^2 / 
                results_hhc_nlp_mm["solution"]["nw"]["1"]["bus"]["$nb"]["vm"]^2) for nb in 0:8]
RMSᵁ    = [sqrt(sum([results_hhc_nlp_mm["solution"]["nw"]["$nh"]["bus"]["$nb"]["vm"] for nh ∈ H].^2)) for nb in 0:8]

## current
RMSᴮ⁻ᶠ  = [sqrt(sum([results_hhc_nlp_ae["solution"]["nw"]["$nh"]["branch"]["$nb"]["cm_fr"] for nh ∈ H].^2)) for nb in 1:4]
RMSᴮ⁻ᵗ  = [sqrt(sum([results_hhc_nlp_ae["solution"]["nw"]["$nh"]["branch"]["$nb"]["cm_to"] for nh ∈ H].^2)) for nb in 1:4]
RMSˣ⁻ᶠ  = [sqrt(sum([results_hhc_nlp_ae["solution"]["nw"]["$nh"]["xfmr"]["$nb"]["cmx_fr"] for nh ∈ H].^2)) for nb in 1:4]
RMSˣ⁻ᵗ  = [sqrt(sum([results_hhc_nlp_ae["solution"]["nw"]["$nh"]["xfmr"]["$nb"]["cmx_to"] for nh ∈ H].^2)) for nb in 1:4]

for nh in [3,5,7] 
    for nb in 0:8 
        vr = round(results_hhc_nlp_mm["solution"]["nw"]["$nh"]["bus"]["$nb"]["vr"], digits=5)
        vi = round(results_hhc_nlp_mm["solution"]["nw"]["$nh"]["bus"]["$nb"]["vi"], digits=5)
        println("coordinate (U$nb)    at ($vr,$vi);") 
    end

    for nx in 1:4
        er = round(results_hhc_nlp_mm["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["erx"], digits=5)
        ei = round(results_hhc_nlp_mm["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["eix"], digits=5)
        println("coordinate (E$nx)    at ($er,$ei);") 
    end

    println("-----")
end

# 
er_1 = results_hhc_nlp_mm["solution"]["nw"]["5"]["xfmr"]["1"]["erx"]
ei_1 = results_hhc_nlp_mm["solution"]["nw"]["5"]["xfmr"]["1"]["eix"]

ea_1 = atand(ei_1/er_1)

x    = 0.01 * cosd(ea_1)
y    = 0.01 * sind(ea_1)
va_2 = results_hhc_nlp_mm["solution"]["nw"]["5"]["bus"]["2"]["va"]

er_3 = results_hhc_nlp_mm["solution"]["nw"]["5"]["xfmr"]["3"]["erx"]
ei_3 = results_hhc_nlp_mm["solution"]["nw"]["5"]["xfmr"]["3"]["eix"]

ea_3 = atand(ei_3/er_3)

x    = 0.02 * cosd(ea_3)
y    = 0.02 * sind(ea_3)
va_6 = results_hhc_nlp_mm["solution"]["nw"]["5"]["bus"]["6"]["va"]


er_1 = results_hhc_nlp_mm["solution"]["nw"]["7"]["xfmr"]["1"]["erx"]
ei_1 = results_hhc_nlp_mm["solution"]["nw"]["7"]["xfmr"]["1"]["eix"]

ea_1 = atand(ei_1/er_1) + 180

x    = 0.02 * cosd(ea_1)
y    = 0.02 * sind(ea_1)
va_2 = results_hhc_nlp_mm["solution"]["nw"]["7"]["bus"]["2"]["va"]

er_3 = results_hhc_nlp_mm["solution"]["nw"]["7"]["xfmr"]["3"]["erx"]
ei_3 = results_hhc_nlp_mm["solution"]["nw"]["7"]["xfmr"]["3"]["eix"]

ea_3 = atand(ei_3/er_3) + 180

x    = 0.02 * cosd(ea_3)
y    = 0.02 * sind(ea_3)
va_6 = results_hhc_nlp_mm["solution"]["nw"]["7"]["bus"]["6"]["va"]