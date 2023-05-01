################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
"""
Example considering optimal harmonic power flow for an industrial power system
taken from:
>   Harmonic Optimal Power Flow with Transformer Excitation by F. Geth and T. 
    Van Acker, pg. 7, § IV.B.
"""

# load pkgs
using HarmonicPowerModels
using Dierckx, Ipopt, JuMP, PowerModels

# set the solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

# read-in data
path = joinpath(HarmonicPowerModels.BASE_DIR,"test/data/matpower/industrial_network_hopf.m")
data = PowerModels.parse_file(path)

# set the formulation
form = PowerModels.IVRPowerModel

# BH-curve
B⁺ = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
H⁺ = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
B = vcat(reverse(-B⁺),0.0,B⁺)
H = vcat(reverse(-H⁺),0.0,H⁺) 
BH_powercore_h100_23 = Dierckx.Spline1D(B, H; k=3, bc="nearest")

# xfmr magnetizing data
magn = Dict("Hᴱ"    => [1, 5], 
            "Hᴵ"    => [1, 3, 5, 7, 9, 13],
            "Emax"  => 1.1,
            "IDH"   => [1.0, 0.06],
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
                                        "Vbase" => 10000),
                            4 => Dict(  "l"     => 1.0,
                                        "A"     => 0.001,
                                        "N"     => 240,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 690)
                            )
            )

# build the harmonic data, dictionary with tuple key where,
# key[1] : xfmr_magn - 0 = without xfmr_magn / 1 = with xfmr_magn
# key[2] : active filter - 0 = without active filter / 1 = with active filter
hdata = Dict()
hdata[(0,1)] = HarmonicPowerModels.replicate(data)
hdata[(1,1)] = HarmonicPowerModels.replicate(data, xfmr_magn=magn)

# build the harmonic data without active filter
for (key,val) in hdata
    hdata[(key[1],0)] = deepcopy(val)
    for (nw,ntw) in hdata[(key[1],0)]["nw"], (ng,gen) in ntw["gen"] 
        if gen["isfilter"] == 1
            delete!(ntw["gen"],ng)
        end
    end 
end 


# add the individual harmonic distortion limits                                 # @F: this should preferably be integrated in the date file
ihdmax = Dict("1" => 1.10, "3" => 0.05, "5" => 0.06, "7" => 0.05, "9" => 0.015, "13" => 0.03)
for (key,val) in hdata
    for (nw,ntw) in val["nw"], (nb,bus) in ntw["bus"]
        bus["ihdmax"] = ihdmax[nw]
    end 
end

# harmonic power flow - with xfmr_magn, without active filter → key = (1,0)
results_hpf = HarmonicPowerModels.solve_hpf(hdata[(1,0)], form, solver)
println("Results for the harmonic power flow")
println("Fundamental harmonic:")
print_summary(results_hpf["solution"]["nw"]["1"])
println("Third harmonic:")
print_summary(results_hpf["solution"]["nw"]["3"])
println("Fifth harmonic:")
print_summary(results_hpf["solution"]["nw"]["5"])
println("Seventh harmonic:")
print_summary(results_hpf["solution"]["nw"]["7"])
println("Nineth harmonic:")
print_summary(results_hpf["solution"]["nw"]["9"])
println("Thirteen harmonic:")
print_summary(results_hpf["solution"]["nw"]["13"])

# harmonic optimal power flow - without xfmr_magn, with active filter → key = (0,1)
results_hopf_wo_magn = HarmonicPowerModels.solve_hopf(hdata[(0,1)], form, solver)
println("Results for the harmonic optimal power flow without xfmr magnetization")
println("Fundamental harmonic:")
print_summary(results_hopf_wo_magn["solution"]["nw"]["1"])
println("Third harmonic:")
print_summary(results_hopf_wo_magn["solution"]["nw"]["3"])
println("Fifth harmonic:")
print_summary(results_hopf_wo_magn["solution"]["nw"]["5"])
println("Seventh harmonic:")
print_summary(results_hopf_wo_magn["solution"]["nw"]["7"])
println("Nineth harmonic:")
print_summary(results_hopf_wo_magn["solution"]["nw"]["9"])
println("Thirteen harmonic:")
print_summary(results_hopf_wo_magn["solution"]["nw"]["13"])

# harmonic optimal power flow - with xfmr_magn, with active filter → key = (1,1)
results_hopf_w_magn = HarmonicPowerModels.solve_hopf(hdata[(1,1)], form, solver)
println("Fundamental harmonic:")
print_summary(results_hopf_w_magn["solution"]["nw"]["1"])
println("Third harmonic:")
print_summary(results_hopf_w_magn["solution"]["nw"]["3"])
println("Fifth harmonic:")
print_summary(results_hopf_w_magn["solution"]["nw"]["5"])
println("Seventh harmonic:")
print_summary(results_hopf_w_magn["solution"]["nw"]["7"])
println("Nineth harmonic:")
print_summary(results_hopf_w_magn["solution"]["nw"]["9"])
println("Thirteen harmonic:")
print_summary(results_hopf_w_magn["solution"]["nw"]["13"])