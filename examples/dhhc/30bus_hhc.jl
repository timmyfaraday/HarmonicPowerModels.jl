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
using MadNLP, Clarabel, Ipopt
using PrettyTables
using Gurobi
using Plots

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_nlp = Ipopt.Optimizer
solver_soc = Clarabel.Optimizer

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")
data = PMs.parse_file(path)


for (b, bus) in data["bus"]
    bus["std"] = "IEC61000-2-4:2002, Cl. 2"
    bus["ref_angle"] = 0.0
    bus["vmax"]  = 1.1
    bus["vmin"]  = 0.9
end

for (g,gen) in data["gen"]
  gen["inf"] = 1
  gen["gsc"] = 0.0
  gen["bsc"] = 0.0
  gen["xr_ratio"] = 30.0
end

# data["hload"] = data["load"]


# function keep_fundamental_load!(hdata)
#     for  n in keys(hdata["nw"])
#       if n ≠ "1"
#         for (l, load) in hdata["nw"][n]["hload"]
#           delete!(hdata["nw"][n]["hload"], l)
#         end
#       end
#     end
#     return hdata
# end
# results_opf = PowerModels.solve_opf(data, ACPPowerModel, solver_nlp)["solution"]

# for (b, bus) in data["bus"]
#     bus["vm"] = results_opf["bus"]["$b"]["vm"]
#     bus["va"] = results_opf["bus"]["$b"]["va"]
# end

# for (g, gen) in data["gen"]
#   gen["pg"] = results_opf["gen"]["$g"]["pg"]
#   gen["qg"] = results_opf["gen"]["$g"]["qg"]
# end

# for (s, shunt) in data["shunt"]
#   delete!(data["shunt"], s)
# end

# set the ref bus to a clean bus
data["bus"]["1"]["bus_type"] = 4

# define the set of considered harmonics
H = [1, 3, 5, 7]

# COMPUTATION ##################################################################
# absolute equality (ae) ######################################################
hdata_ae = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "absolute equality")
results_ae_nlp = HPM.solve_hhc(hdata_ae, HarmonicPowerModel, solver_nlp)
results_ae_soc = HPM.solve_hhc(hdata_ae, dHHCPowerModel, solver_nlp; optimizer_soc = solver_soc)

## maximum efficiency (me) #####################################################
hdata_me = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "maximum efficiency")
results_me_nlp = HPM.solve_hhc(hdata_me, HarmonicPowerModel, solver_nlp)
results_me_soc = HPM.solve_hhc(hdata_me, dHHCPowerModel, solver_nlp; optimizer_soc = solver_soc)

## maximin (mm) ####################################################`############
hdata_mm = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "maximin")
results_mm_nlp = HPM.solve_hhc(hdata_mm, HarmonicPowerModel, solver_nlp)
results_mm_soc = HPM.solve_hhc(hdata_mm, dHHCPowerModel, solver_nlp; optimizer_soc = solver_soc)

## Kalai-Smorodinsky bargaining (ks) ###########################################
# Calculate impedance
# Zh = calculate_pos_seq_harmonic_impedance(data, collect(1.0:0.2:50.0), collect(1:9))

# # solve HHC problem -- NLP
# hdata_nlp_ks = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "Kalai-Smorodinsky bargaining")
# HPM.calculate_maximum_harmonic_source_current_injection!(hdata_nlp_ks, Zh)
# results_hhc_nlp_ks = HPM.solve_hhc(hdata_nlp_ks, HarmonicPowerModel, solver_nlp)

# # solve HHC problem -- SOC 
# hdata_soc_ks = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, Zh)
# results_hhc_soc_ks = HPM.solve_hhc(hdata_soc_ks, dHHCPowerModel, solver_nlp)


header  = (
            ["obj.", "mod.", "obj. value", "solve time [s]"]
          );
table_data    = vcat( ["abs. eq." "nl" results_ae_nlp["objective"] results_ae_nlp["solve_time"]],
                ["abs. eq." "soc" results_ae_soc["objective"] results_ae_soc["solve_time"]],
                ["max. eff." "nl" results_me_nlp["objective"] results_me_nlp["solve_time"]],
                ["max. eff." "soc" results_me_soc["objective"] results_me_soc["solve_time"]],
                ["maximin" "nl" results_mm_nlp["objective"] results_mm_nlp["solve_time"]],
                ["maximin" "soc" results_mm_soc["objective"] results_mm_soc["solve_time"]],
                #["KS barg." "nl" results_hhc_nlp_ks["objective"] results_hhc_nlp_ks["solve_time"]],
                #["KS barg." "soc" results_hhc_soc_ks["objective"] results_hhc_soc_ks["solve_time"]]
                )

pretty_table(table_data, header=header)




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

# header  = (
#             ["obj.", "mod.", "obj. value", "solve time [s]"]
#           );
# table_data    = vcat( ["abs. eq." "nl" results_ae_nlp["objective"] results_ae_nlp["solve_time"]],
#                 ["abs. eq." "soc" results_ae_soc["objective"] results_ae_soc["solve_time"]],
#                 ["max. eff." "nl" results_me_nlp["objective"] results_me_nlp["solve_time"]],
#                 ["max. eff." "soc" results_me_soc["objective"] results_me_soc["solve_time"]],
#                 ["maximin" "nl" results_mm_nlp["objective"] results_mm_nlp["solve_time"]],
#                 ["maximin" "soc" results_mm_soc["objective"] results_mm_soc["solve_time"]],
#                 #["KS barg." "nl" results_hhc_nlp_ks["objective"] results_hhc_nlp_ks["solve_time"]],
#                 #["KS barg." "soc" results_hhc_soc_ks["objective"] results_hhc_soc_ks["solve_time"]]
#                 )

# pretty_table(table_data, header=header)




hdata = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc, hhc_principle = "absolute equality")

hpf_data = deepcopy(hdata)
for n in keys(hpf_data["nw"])
    if n ≠ "1"
        delete!(hpf_data["nw"], n)
    end
end


hpf= HPM.solve_hpf(hpf_data, HarmonicPowerModel, solver_nlp)["solution"]["nw"]["1"]
hopf= HPM.solve_hopf(hpf_data, HarmonicPowerModel, solver_nlp)["solution"]["nw"]["1"]
hpf1= HPM.solve_hhc_init(hpf_data, HarmonicPowerModel, solver_nlp)["solution"]["nw"]["1"]
hhc= HPM.solve_hhc(hdata, HarmonicPowerModel, solver_nlp)["solution"]["nw"]["1"]
# hopf_opf = HPM.solve_hopf(hpf_data, HarmonicPowerModel, solver_nlp)["solution"]["nw"]["1"]

for (n, bus) in hpf["bus"]
   println("hhc - hpf, bus: ",n ," ," , sqrt(hhc["bus"][n]["vbr"]^2 + hhc["bus"][n]["vbi"]^2) - sqrt(hpf1["bus"][n]["vbr"]^2 + hpf1["bus"][n]["vbi"]^2))
end


for (s, src) in hhc["hsrc"]
   println("hhc - hpf, src: ",s ," ," , sqrt(hhc["hsrc"][s]["crr"]^2 + hhc["hsrc"][s]["cri"]^2) - sqrt(hpf1["hsrc"][s]["crr"]^2 + hpf1["hsrc"][s]["cri"]^2))
end


for (g, gen) in hhc["gen"]
   println("hhc - hpf, gen: ",g ," ," , hhc["gen"][g]["cgm"] - hpf1["gen"][g]["cgm"])
end






sum([src["crr"] for (s, src) in results_me_nlp["solution"]["nw"]["1"]["hsrc"]])
sum([gen["cgr"] for (g, gen) in results_me_nlp["solution"]["nw"]["1"]["gen"]])

sum([src["crr"] for (s, src) in results_me_nlp["solution"]["nw"]["1"]["hsrc"]]) / sum([gen["cgr"] for (g, gen) in results_me_nlp["solution"]["nw"]["1"]["gen"]])


sum([gen["pg"] for (g, gen) in results_me_nlp["solution"]["nw"]["1"]["gen"]])
sum([src["p_fund"] for (s, src) in hdata_me["nw"]["1"]["hsrc"]])


ib = [max(branch["cbm_fr"], branch["cbm_to"]) / hpf_data["nw"]["1"]["branch"][b]["i_rms_max"] for  (b, branch) in hpf["branch"]]
Plots.plot(ib)

ub = [bus["vbm"] for (b, bus) in hpf["bus"]]
Plots.plot(ub)

ig = [abs(gen["cgm"]) / hpf_data["nw"]["1"]["gen"][g]["i_rms_max"] for (g, gen) in hpf["gen"]]
Plots.plot(ig)


sum([gen["pg"] for (g, gen) in hpf["gen"]])
sum([src["p_fund"] for (s, src) in hpf_data["nw"]["1"]["hsrc"]])