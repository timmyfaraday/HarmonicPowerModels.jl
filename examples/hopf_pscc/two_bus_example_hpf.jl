################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Example considering harmonic power flow for a two-bus example network taken  # 
# from: Harmonic Optimal Power Flow with Transformer Excitation by F. Geth and #
# T. Van Acker, pg. 7, § IV.A.                                                 #
################################################################################
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
path = joinpath(HPM.BASE_DIR,"test/data/matpower/two_bus_example_hpf.m")
data = PMs.parse_file(path)

# set the ref bus to a clean bus
data["bus"]["1"]["bus_type"] = 4

# solve PF problem
results_fund = PMs.solve_pf_iv(data, PMs.IVRPowerModel, solver)

# define the set of considered harmonics
H = [1, 3]

# COMPUTATION ##################################################################
# solve HPF problem
hdata   = build_hdata_from_matpower_file(data, H=H, prob=:hpf)
results = solve_hpf(hdata, HarmonicPowerModel, solver)

# RESULTS ######################################################################
# TABLE: Bus results for two bus line case (per unit)

## bus results
Um(nh,nb)   = round(abs(results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"] +
                        results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"] * im), digits=3);
Ua(nh,nb)   = round(atand(results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"],
                          results["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"]), digits=3);
RMS(nb)     = round(sqrt(sum(Um(nh,nb)^2 for nh in H)), digits=3);
THD(nb)     = round(sqrt(sum(Um(nh,nb)^2 for nh in H if nh ≠ 1) / Um(1,nb)^2), digits=3);

header  = (
            ["bus", "harmonic", "|Uᵢₕ|", "∠Uᵢₕ", "RMSᵢ", "THDᵢ"]
          );
data    = vcat( hcat("i",   "h=1",  Um(1,1),    Ua(1,1),    RMS(1),     THD(1)),
                hcat("",    "h=3",  Um(3,1),    Ua(3,1),    "",         ""),
                hcat("j",   "h=1",  Um(1,2),    Ua(1,2),    RMS(2),     THD(2)),
                hcat("",    "h=3",  Um(3,2),    Ua(3,2),    "",         ""));
pretty_table(data, header=header)

## load results 
Id(nh,nl)   = round(results["solution"]["nw"]["$nh"]["hload"]["$nl"]["clr"] +
                    results["solution"]["nw"]["$nh"]["hload"]["$nl"]["cli"] * im, digits=3);
Sd(nh,nl)   = (results["solution"]["nw"]["$nh"]["bus"]["2"]["vbr"] +
               results["solution"]["nw"]["$nh"]["bus"]["2"]["vbi"] * im) * conj(Id(nh,nl));
Pd(nh,nl)   = round(real(Sd(nh,nl)), digits=3);
Qd(nh,nl)   = round(imag(Sd(nh,nl)), digits=3);
Idm(nh,nl)  = round(abs(Id(nh,nl)), digits=3);

header  = (
            ["load", "harmonic", "Pᵈₕ", "Qᵈₕ", "Iᵈₕ", "|Iᵈₕ|"]
          );
data    = vcat( hcat("d",   "h=1",  Pd(1,1),    Qd(1,1),    Id(1,1),    Idm(1,1)),
                hcat("",    "h=3",  Pd(3,1),    Qd(3,1),    Id(3,1),    Idm(3,1)));
pretty_table(data, header=header)

## generator results
Ig(nh,ng)   = round(results["solution"]["nw"]["$nh"]["gen"]["$ng"]["cgr"] +
                    results["solution"]["nw"]["$nh"]["gen"]["$ng"]["cgi"] * im, digits=3);
Sg(nh,ng)   = (results["solution"]["nw"]["$nh"]["bus"]["1"]["vbr"] +
               results["solution"]["nw"]["$nh"]["bus"]["1"]["vbi"] * im) * conj(Ig(nh,ng));
Pg(nh,ng)   = round(real(Sg(nh,ng)), digits=3);
Qg(nh,ng)   = round(imag(Sg(nh,ng)), digits=3);
Igm(nh,ng)  = round(abs(Ig(nh,ng)), digits=3);

header  = (
            ["gen.", "harmonic", "Pᵍₕ", "Qᵍₕ", "Iᵍₕ", "|Iᵍₕ|"]
          );
data    = vcat( hcat("g",   "h=1",  Pg(1,1),    Qg(1,1),    Ig(1,1),    Igm(1,1)),
                hcat("",    "h=3",  Pg(3,1),    Qg(3,1),    Ig(3,1),    Igm(3,1)));
pretty_table(data, header=header)

## branch results
Ib_fr(nh,nb)   = round(results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbr_fr"] +
                       results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbi_fr"] * im, digits=3);
Sb_fr(nh,nb)   = (results["solution"]["nw"]["$nh"]["bus"]["1"]["vbr"] +
                  results["solution"]["nw"]["$nh"]["bus"]["1"]["vbi"] * im) * conj(Ib_fr(nh,nb));
Pb_fr(nh,nb)    = round(real(Sb_fr(nh,nb)), digits=3);
Qb_fr(nh,nb)    = round(imag(Sb_fr(nh,nb)), digits=3);
Ib_to(nh,nb)   = round(results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbr_to"] +
                       results["solution"]["nw"]["$nh"]["branch"]["$nb"]["cbi_to"] * im, digits=3);
Sb_to(nh,nb)   = (results["solution"]["nw"]["$nh"]["bus"]["2"]["vbr"] +
                  results["solution"]["nw"]["$nh"]["bus"]["2"]["vbi"] * im) * conj(Ib_to(nh,nb));
Pb_to(nh,nb)    = round(real(Sb_to(nh,nb)), digits=3);
Qb_to(nh,nb)    = round(imag(Sb_to(nh,nb)), digits=3);
Ploss(nh,nb)    = Pb_fr(nh,nb) + Pb_to(nh,nb)

header  = (
            ["line", "harmonic", "Pˡⁱʲₕ", "Qˡⁱʲₕ", "Iˡⁱʲₕ", "Pˡᵒˢˢₕ"]
          );
data    = vcat( hcat("g",   "h=1",  Pb_fr(1,1), Qb_fr(1,1), Ib_fr(1,1), Ploss(1,1)),
                hcat("",    "h=3",  Pb_fr(3,1), Qb_fr(3,1), Ib_fr(3,1), Ploss(3,1)));
pretty_table(data, header=header)