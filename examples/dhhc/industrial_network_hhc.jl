################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Example considering the harmonic hosting capacity of an industrial power     #
# system taken from: Exact Lower Bound on Equitable Harmonic Hosting Capacity  #
# by T. Van Acker and H. Ergun, pg. 6, § III.A.                                #
# ---------------------------------------------------------------------------- #
# This numerical illustration considers a section of a real-world industrial   #
# phase-balanced three-phase power system across five voltage levels, with     #
# four transformers, three cables, and four harmonic units.                    #
# On top of the fundamental harmonic, the considered harmonic set is           # 
# H ∈ {3, 5, 7}, i.e., one harmonic from each sequence subset. The aim of this # 
# numerical illustration is to                                                 #
#   1) validate the proposed formulations; and                                 #
#   2) show the equivalence between the non-linear and second-order cone       #
#      models.                                                                 #
# For the purpose of proving equivalency, both models are solved using Ipopt,  # 
# where the second-order constraints are translated to their equivalent        #
# non-convex quadratic form using the appropriate MATHOPTINTERFACE constraint  #
# bridge.                                                                      #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.3.0 - example template                                                    # 
################################################################################

# INPUT ########################################################################
# using pkgs
using HarmonicPowerModels, PowerModels
using Ipopt 
using PrettyTables
using Clarabel
using Gurobi

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_nlp = Ipopt.Optimizer
solver_soc = Clarabel.Optimizer

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
data = PMs.parse_file(path)

# define the set of considered harmonics
H = [1, 3, 5, 7]

# TODO: Check the hsrc constraint -> casues convergence to infeasible point in OPF.....

# COMPUTATION ##################################################################
# absolute equality (ae) ######################################################

# solve HHC problem -- NLP
hdata_nlp_ae = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "absolute equality")
results_hhc_nlp_ae = HPM.solve_hhc(hdata_nlp_ae, HarmonicPowerModel, solver_nlp)

# solve HHC problem -- SOC
hdata_soc_ae = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "absolute equality")
results_hhc_soc_ae = HPM.solve_hhc(hdata_soc_ae, dHHCPowerModel, solver_nlp; optimizer_soc = solver_soc)

## maximum efficiency (me) #####################################################

# solve HHC problem -- NLP
hdata_nlp_me = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "maximum efficiency")
results_hhc_nlp_me = HPM.solve_hhc(hdata_nlp_me, HarmonicPowerModel, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc_me = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "maximum efficiency")
results_hhc_soc_me = HPM.solve_hhc(hdata_soc_me, dHHCPowerModel, solver_nlp; optimizer_soc = solver_soc)

## maximin (mm) ################################################################

# solve HHC problem -- NLP
hdata_nlp_mm = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "maximin")
results_hhc_nlp_mm = HPM.solve_hhc(hdata_nlp_mm, HarmonicPowerModel, solver_nlp)

# solve HHC problem -- SOC 
hdata_soc_mm = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "maximin")
results_hhc_soc_mm = HPM.solve_hhc(hdata_soc_mm, dHHCPowerModel, solver_nlp; optimizer_soc = solver_soc)

## Kalai-Smorodinsky bargaining (ks) ###########################################
# # Calculate impedance
# Zh = calculate_pos_seq_harmonic_impedance(data, collect(50:50.0:H[end] * 50.0), collect(1:length(data["bus"])))

# # solve HHC problem -- NLP
# hdata_nlp_ks = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "Kalai-Smorodinsky bargaining")
# HPM.calculate_maximum_harmonic_source_current_injection!(hdata_nlp_ks, Zh)
# results_hhc_nlp_ks = HPM.solve_hhc(hdata_nlp_ks, HarmonicPowerModel, solver_nlp)

# # solve HHC problem -- SOC 
# hdata_soc_ks = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, Zh)
# results_hhc_soc_ks = HPM.solve_hhc(hdata_soc_ks, dHHCPowerModel, solver_nlp)

# RESULTS ######################################################################
# TABLE: Harmonic current injection from the non-linear model with the 
# Kalai-Smorodinsky bargaining fairness objective 

# Im(nh,nl) = round(results_hhc_nlp_ks["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"], digits=7);
# Ia(nh,nl) = round(atand(results_hhc_nlp_ks["solution"]["nw"]["$nh"]["load"]["$nl"]["cid"],
#                         results_hhc_nlp_ks["solution"]["nw"]["$nh"]["load"]["$nl"]["crd"]), digits=2);

# header  = (
#             ["", "h=3", "h=3", "h=5", "h=5", "h=7", "h=7"],
#             ["unit", "|Ī| [pu]", "∠Ī [°]", "|Ī| [pu]", "∠Ī [°]", "|Ī| [pu]", "∠Ī [°]"]
#           );
# data    = vcat( hcat("u₁", [Im(3,1), Ia(3,1), Im(5,1), Ia(5,1), Im(7,1), Ia(7,1)]'),
#                 hcat("u₂", [Im(3,2), Ia(3,2), Im(5,2), Ia(5,2), Im(7,2), Ia(7,2)]'),
#                 hcat("u₃", [Im(3,3), Ia(3,3), Im(5,3), Ia(5,3), Im(7,3), Ia(7,3)]'),
#                 hcat("u₄", [Im(3,4), Ia(3,4), Im(5,4), Ia(5,4), Im(7,4), Ia(7,4)]'));

# pretty_table(data, header=header)

# FIGURE: The bus and excitation voltage phasors throughout the industrial
# power system for (a) zero, (b) negative, and (c) positive sequence,
# respectively, for the non-linear model with the Kalai-Smorodinsky bargaining
# fairness objective

# TODO

# for nh in [3,5,7] 
#     for nb in 0:8 
#         vr = round(results_hhc_nlp_ks["solution"]["nw"]["$nh"]["bus"]["$nb"]["vr"], digits=5)
#         vi = round(results_hhc_nlp_ks["solution"]["nw"]["$nh"]["bus"]["$nb"]["vi"], digits=5)
#         println("\\coordinate (U$nb)    at ($vr,$vi);") 
#     end
#     for nx in 1:4
#         er = round(results_hhc_nlp_ks["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["erx"], digits=5)
#         ei = round(results_hhc_nlp_ks["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["eix"], digits=5)
#         println("\\coordinate (E$nx)    at ($er,$ei);") 
#     end
#     println("-----")
# end

# TABLE: Objective value and solve times of non-linear and second-order cone 
# models, for each considered fairness objective 

header  = (
            ["obj.", "mod.", "obj. value", "solve time [s]"]
          );
table_data    = vcat( ["abs. eq." "nl" results_hhc_nlp_ae["objective"] results_hhc_nlp_ae["solve_time"]],
                ["abs. eq." "soc" results_hhc_soc_ae["objective"] results_hhc_soc_ae["solve_time"]],
                ["max. eff." "nl" results_hhc_nlp_me["objective"] results_hhc_nlp_me["solve_time"]],
                ["max. eff." "soc" results_hhc_soc_me["objective"] results_hhc_soc_me["solve_time"]],
                ["maximin" "nl" results_hhc_nlp_mm["objective"] results_hhc_nlp_mm["solve_time"]],
                ["maximin" "soc" results_hhc_soc_mm["objective"] results_hhc_soc_mm["solve_time"]],
                #["KS barg." "nl" results_hhc_nlp_ks["objective"] results_hhc_nlp_ks["solve_time"]],
                #["KS barg." "soc" results_hhc_soc_ks["objective"] results_hhc_soc_ks["solve_time"]]
                )

pretty_table(table_data, header=header)



#   hpf_data = deepcopy(hdata_nlp_ae)
#   for n in keys(hpf_data["nw"])
#       if n ≠ "1"
#           delete!(hpf_data, n)
#       end
#   end

#   # solve hpf problem for the fundamental harmonic only
#   hpf_results = solve_hopf(hpf_data, HarmonicPowerModel, solver_nlp)

# for (b, bus) in results_hhc_nlp_ae["solution"]["nw"]["1"]["bus"]
#    println(b, " " ,sqrt(bus["vbr"]^2 + bus["vbi"]^2))
#    println(b, " " ,sqrt(hpf_results["solution"]["nw"]["1"]["bus"][b]["vbr"]^2 + hpf_results["solution"]["nw"]["1"]["bus"][b]["vbi"]^2))
# end