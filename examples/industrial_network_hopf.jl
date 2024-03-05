################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

"""
Example considering optimal harmonic power flow for an industrial power system
taken from:
>   Harmonic Optimal Power Flow with Transformer Excitation by F. Geth and T. 
    Van Acker, pg. 7, § IV.B.
Please note that the current implementation is no longer exactly corresponding 
to the one presented in the paper.
- short-circuit impedance of the inf. grid added
- ihd and thd limits relative to standard, rather than hard coded
"""

# load pkgs
using HarmonicPowerModels, PowerModels
using Dierckx
using Ipopt

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver = Ipopt.Optimizer

# read-in data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hopf.m")
data = PowerModels.parse_file(path)

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

# solve HOPF problem w/o xfmr magnitization
hdata_wo = HPM.replicate(data, H=H, bus_id=6)
results_wo = HPM.solve_hopf(hdata_wo, PMs.IVRPowerModel, solver)

# solve HOPF problem w. xfmr magnitization
hdata_w = HPM.replicate(data, H=H, xfmr_magn=magn, bus_id=6)
results_w = HPM.solve_hopf(hdata_wo, PMs.IVRPowerModel, solver)