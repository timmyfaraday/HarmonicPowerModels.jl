################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Example considering harmonic power flow for a two-bus example network taken  # 
# from: Harmonic Optimal Power Flow with Transformer Excitation by F. Geth and #
# T. Van Acker, pg. 7, § IV.A.                                                 #
# ---------------------------------------------------------------------------- #   
# Note that Table II contains a mistake, Qˡⁱʲ₃ ≠ 0.006, rather Qˡⁱʲ₃ ≠ 0.0,    #
# respecting the reactive power balance for the third harmonic at bus 1.       # 
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
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

# pkg const
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver = Ipopt.Optimizer

# read-in data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/two_bus_example_hopf.m")
data = PMs.parse_file(path)

# set the ref bus to a clean bus
data["bus"]["1"]["bus_type"] = 4

# define the set of considered harmonics
H = [1, 2, 3]

# COMPUTATION ##################################################################
# solve HPF problem
hdata       = build_hdata_from_matpower_file(data, H=H, prob=:hpf)
results     = solve_hpf(hdata, HarmonicPowerModel, solver)

# RESULTS PROCESSING ###########################################################
## TABLE V: Case 1 - Bus voltages [pu/°], 3rd and 9th harmonic are zero.
Um(nh,nb)   = round(abs(results_wo["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"] +
                        results_wo["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"] * im), digits=3);
Ua(nh,nb)   = round(atand(results_wo["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"],
                          results_wo["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"]), digits=3);
THD(nb)     = round(sqrt(sum(Um(nh,nb)^2 for nh in H if nh ≠ 1) / Um(1,nb)^2), digits=3);
header  = (
            ["i", "|Uᵢ₁|", "∠Uᵢ₁", "|Uᵢ₅|", "∠Uᵢ₅", "|Uᵢ₇|", "∠Uᵢ₇", "|Uᵢ₁₃|", "∠Uᵢ₁₃", "THDᵢ"]
          );
data    = vcat( hcat("1", Um(1,1), Ua(1,1), Um(5,1), Ua(5,1), Um(7,1), Ua(7,1), Um(13,1), Ua(13,1), THD(1)),
                hcat("2", Um(1,2), Ua(1,2), Um(5,2), Ua(5,2), Um(7,2), Ua(7,2), Um(13,2), Ua(13,2), THD(2)),
                hcat("3", Um(1,3), Ua(1,3), Um(5,3), Ua(5,3), Um(7,3), Ua(7,3), Um(13,3), Ua(13,3), THD(3)),
                hcat("4", Um(1,4), Ua(1,4), Um(5,4), Ua(5,4), Um(7,4), Ua(7,4), Um(13,4), Ua(13,4), THD(4)),
                hcat("5", Um(1,5), Ua(1,5), Um(5,5), Ua(5,5), Um(7,5), Ua(7,5), Um(13,5), Ua(13,5), THD(5)),
                hcat("6", Um(1,6), Ua(1,6), Um(5,6), Ua(5,6), Um(7,6), Ua(7,6), Um(13,6), Ua(13,6), THD(6)),
                hcat("7", Um(1,7), Ua(1,7), Um(5,7), Ua(5,7), Um(7,7), Ua(7,7), Um(13,7), Ua(13,7), THD(7)),
                hcat("8", Um(1,8), Ua(1,8), Um(5,8), Ua(5,8), Um(7,8), Ua(7,8), Um(13,8), Ua(13,8), THD(8)));
pretty_table(data, header=header)


# TABLE II: Bus results for two bus line case (per unit)
Um(nh,nb)   = round(abs(results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"] +
                        results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"] * im), digits=3);
Ua(nh,nb)   = round(atand(results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"],
                          results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"]), digits=3);
RMS(nb)     = round(sqrt(sum(Um(nh,nb)^2 for nh in H)), digits=3);
THD(nb)     = round(sqrt(sum(Um(nh,nb)^2 for nh in H if nh ≠ 1) / Um(1,nb)^2), digits=3);

header_b    = (["bus", "harmonic", "|Uᵢₕ|", "∠Uᵢₕ", "RMSᵢ", "THDᵢ"]);
data_b      = vcat( hcat("i",   "h=1",  Um(1,1),    Ua(1,1),    RMS(1),     THD(1)),
                    hcat("",    "h=3",  Um(3,1),    Ua(3,1),    "",         ""),
                    hcat("j",   "h=1",  Um(1,2),    Ua(1,2),    RMS(2),     THD(2)),
                    hcat("",    "h=3",  Um(3,2),    Ua(3,2),    "",         ""));

## load results 
Id(nh,nl)   = round(results["solution"]["nw"]["$nh"]["hload"]["$nl"]["clr"] +
                    results["solution"]["nw"]["$nh"]["hload"]["$nl"]["cli"] * im, digits=3);
Sd(nh,nl)   = (results["solution"]["nw"]["$nh"]["bus"]["2"]["vbr"] +
               results["solution"]["nw"]["$nh"]["bus"]["2"]["vbi"] * im) * conj(Id(nh,nl));
Pd(nh,nl)   = round(real(Sd(nh,nl)), digits=3);
Qd(nh,nl)   = round(imag(Sd(nh,nl)), digits=3);
Idm(nh,nl)  = round(abs(Id(nh,nl)), digits=3);

header_l    = (["load", "harmonic", "Pᵈₕ", "Qᵈₕ", "Iᵈₕ", "|Iᵈₕ|"]);
data_l      = vcat( hcat("d",   "h=1",  Pd(1,1),    Qd(1,1),    Id(1,1),    Idm(1,1)),
                    hcat("",    "h=3",  Pd(3,1),    Qd(3,1),    Id(3,1),    Idm(3,1)));

## generator results
Ig(nh,ng)   = round(results["solution"]["nw"]["$nh"]["gen"]["$ng"]["cgr"] +
                    results["solution"]["nw"]["$nh"]["gen"]["$ng"]["cgi"] * im, digits=3);
Sg(nh,ng)   = (results["solution"]["nw"]["$nh"]["bus"]["1"]["vbr"] +
               results["solution"]["nw"]["$nh"]["bus"]["1"]["vbi"] * im) * conj(Ig(nh,ng));
Pg(nh,ng)   = round(real(Sg(nh,ng)), digits=3);
Qg(nh,ng)   = round(imag(Sg(nh,ng)), digits=3);
Igm(nh,ng)  = round(abs(Ig(nh,ng)), digits=3);

header_g    = (["gen.", "harmonic", "Pᵍₕ", "Qᵍₕ", "Iᵍₕ", "|Iᵍₕ|"]);
data_g      = vcat( hcat("g",   "h=1",  Pg(1,1),    Qg(1,1),    Ig(1,1),    Igm(1,1)),
                    hcat("",    "h=3",  Pg(3,1),    Qg(3,1),    Ig(3,1),    Igm(3,1)));

## branch results
Ib_fr(nh,nb)  = round(results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbr_fr"] +
                      results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbi_fr"] * im, digits=3);
Sb_fr(nh,nb)  = (results["solution"]["nw"]["$nh"]["bus"]["1"]["vbr"] +
                 results["solution"]["nw"]["$nh"]["bus"]["1"]["vbi"] * im) * conj(Ib_fr(nh,nb));
Pb_fr(nh,nb)  = round(real(Sb_fr(nh,nb)), digits=3);
Qb_fr(nh,nb)  = round(imag(Sb_fr(nh,nb)), digits=3);
Ib_to(nh,nb)  = round(results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbr_to"] +
                      results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbi_to"] * im, digits=3);
Sb_to(nh,nb)  = (results["solution"]["nw"]["$nh"]["bus"]["2"]["vbr"] +
                  results["solution"]["nw"]["$nh"]["bus"]["2"]["vbi"] * im) * conj(Ib_to(nh,nb));
Pb_to(nh,nb)  = round(real(Sb_to(nh,nb)), digits=3);
Qb_to(nh,nb)  = round(imag(Sb_to(nh,nb)), digits=3);
Ploss(nh,nb)  = Pb_fr(nh,nb) + Pb_to(nh,nb)

header_br     = (["line", "harmonic", "Pˡⁱʲₕ", "Qˡⁱʲₕ", "Iˡⁱʲₕ", "Pˡᵒˢˢₕ"]);
data_br       = vcat( hcat("g",   "h=1",  Pb_fr(1,1), Qb_fr(1,1), Ib_fr(1,1), Ploss(1,1)),
                      hcat("",    "h=3",  Pb_fr(3,1), Qb_fr(3,1), Ib_fr(3,1), Ploss(3,1)));

## PRINT #######################################################################
println("TABLE II: Bus results for two bus line case (per unit)")
pretty_table(data_b, header=header_b)
pretty_table(data_l, header=header_l)
pretty_table(data_g, header=header_g)
pretty_table(data_br, header=header_br)

# COMPUTATION ##################################################################
# solve HPF problem
hdata       = build_hdata_from_matpower_file(data, H=H, prob=:hpf)
results_wo  = solve_hpf(hdata, HarmonicPowerModel, solver)

