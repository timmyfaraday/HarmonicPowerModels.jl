################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Example considering optimal harmonic power flow for an industrial power      #
# system taken from: Harmonic Optimal Power Flow with Transformer Excitation   #
# by F. Geth and T. Van Acker, pg. 7, § IV.B.                                  #
# ---------------------------------------------------------------------------- #   
# Please note that the current implementation is no longer exactly             #
# corresponding to the one presented in the paper.                             #
# - ihd and thd limits relative to standard, rather than hard coded            #
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
using Dierckx
using PrettyTables

# pkg const
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver = Ipopt.Optimizer

# read-in data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hopf.m")
data = PowerModels.parse_file(path)

# set the ref bus to a clean bus
data["bus"]["1"]["bus_type"] = 4

# define the set of considered harmonics
H = [1, 3, 5, 7, 9, 13]

# build xfmr magnetization data
B⁺  = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
H⁺  = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
Bᵗ  = vcat(reverse(-B⁺),0.0,B⁺)
Hᵗ  = vcat(reverse(-H⁺),0.0,H⁺) 
BH_powercore_h100_23 = Dierckx.Spline1D(Bᵗ, Hᵗ; k=3, bc="nearest")
magn = Dict("Hᴱ"    => [1, 5], 
            "Hᴵ"    => [1, 3, 5, 7, 9, 13],
            "Emax"  => 1.1,
            "IHD"   => [1.0, 0.06],
            "pcs"   => [21, 11],
            "xfmr"  => Dict(1 => Dict(  "l"     => 11.4,
                                        "A"     => 0.5,
                                        "N"     => 500,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 150000),
                            2 => Dict(  "l"     => 8.0,
                                        "A"     => 0.2,
                                        "N"     => 300,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 36000),
                            3 => Dict(  "l"     => 3.1,
                                        "A"     => 0.07,
                                        "N"     => 240,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 10000)))

# COMPUTATION ##################################################################
## hpf w. magnetizing current
hdata_w         = build_hdata_from_matpower_file(data, H=H, xfmr_magn=magn, prob=:hopf, bus_id=6)
results_hpf     = solve_hpf(hdata_w, HarmonicPowerModel, solver)

## hopf w. magnetizing current
results_hopf_w  = solve_hopf(hdata_w, HarmonicPowerModel, solver)

## hopf w/o magnetizing current
hdata_wo        = build_hdata_from_matpower_file(data, H=H, prob=:hopf, bus_id=6)
results_hopf_wo = solve_hopf(hdata_wo, HarmonicPowerModel, solver)

# RESULTS ######################################################################
## TABLE V: Case 1 - Bus voltages [pu/°], 3rd and 9th harmonic are zero.
Um(nh,nb)     = round(abs(results_hpf["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"] +
                          results_hpf["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"] * im), digits=3);
Ua(nh,nb)     = round(atand(results_hpf["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbi"],
                            results_hpf["solution"]["nw"]["$nh"]["bus"]["$nb"]["vbr"]), digits=2);
THD(nb)       = round(sqrt(sum(Um(nh,nb)^2 for nh in H if nh ≠ 1) / Um(1,nb)^2), digits=3);
header_hpf_v  = (["i", "|Uᵢ₁|", "∠Uᵢ₁", "|Uᵢ₅|", "∠Uᵢ₅", "|Uᵢ₇|", "∠Uᵢ₇", "|Uᵢ₁₃|", "∠Uᵢ₁₃", "THDᵢ"]);
data_hpf_v    = vcat( hcat("1", Um(1,1), Ua(1,1), Um(5,1), Ua(5,1), Um(7,1), Ua(7,1), Um(13,1), Ua(13,1), THD(1)),
                      hcat("2", Um(1,2), Ua(1,2), Um(5,2), Ua(5,2), Um(7,2), Ua(7,2), Um(13,2), Ua(13,2), THD(2)),
                      hcat("3", Um(1,3), Ua(1,3), Um(5,3), Ua(5,3), Um(7,3), Ua(7,3), Um(13,3), Ua(13,3), THD(3)),
                      hcat("4", Um(1,4), Ua(1,4), Um(5,4), Ua(5,4), Um(7,4), Ua(7,4), Um(13,4), Ua(13,4), THD(4)),
                      hcat("5", Um(1,5), Ua(1,5), Um(5,5), Ua(5,5), Um(7,5), Ua(7,5), Um(13,5), Ua(13,5), THD(5)),
                      hcat("6", Um(1,6), Ua(1,6), Um(5,6), Ua(5,6), Um(7,6), Ua(7,6), Um(13,6), Ua(13,6), THD(6)),
                      hcat("7", Um(1,7), Ua(1,7), Um(5,7), Ua(5,7), Um(7,7), Ua(7,7), Um(13,7), Ua(13,7), THD(7)),
                      hcat("8", Um(1,8), Ua(1,8), Um(5,8), Ua(5,8), Um(7,8), Ua(7,8), Um(13,8), Ua(13,8), THD(8)));

## TABLE VI: Excitation voltages and magnitizing currents [pu/°].
Em(nh,nx)     = round(abs(results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["exr"] +
                          results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["exi"] * im), digits=3);
Ea(nh,nx)     = round(atand(results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["exi"],
                            results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["exr"]), digits=1);
Im(nh,nx)     = round(abs(results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cxmr"] +
                          results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cxmi"] * im), digits=10);
Ia(nh,nx)     = round(atand(results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cxmi"],
                            results_hpf["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cxmr"]), digits=1);
header_hpf_e  = (["x", "|Eₓ₁|", "∠Eₓ₁", "|Eₓ₅|", "|Iᵐₓ₁|", "∠Iᵐₓ₁", "|Iᵐₓ₁|", "|Iᵐₓ₁|", "|Iᵐₓ₁|", "|Iᵐₓ₁|", "|Iᵐₓ₁|"]);
data_hpf_e    = vcat( hcat("1", Em(1,1), Ea(1,1), Em(5,1), Im(1,1), Ia(1,1), Im(3,1), Im(5,1), Im(7,1), Im(9,1), Im(13,1)),
                      hcat("2", Em(1,2), Ea(1,2), Em(5,2), Im(1,2), Ia(1,2), Im(3,2), Im(5,2), Im(7,2), Im(9,2), Im(13,2)),
                      hcat("3", Em(1,3), Ea(1,3), Em(5,3), Im(1,3), Ia(1,3), Im(3,3), Im(5,3), Im(7,3), Im(9,3), Im(13,3)));

## PRINT #######################################################################
println("TABLE V: Case 1 - Bus voltages [pu/°], 3rd and 9th harmonic are zero.")
pretty_table(data_hpf_v, header=header_hpf_v)
println("TABLE VI: Excitation voltages and magnitizing currents [pu/°].")
pretty_table(data_hpf_e, header=header_hpf_e)