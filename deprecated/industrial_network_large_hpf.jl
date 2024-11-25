################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

# using pkgs
using HarmonicPowerModels, PowerModels
using Ipopt 
using Plots
using Revise

# pkg cte
const PMs = PowerModels
const HPM = HarmonicPowerModels

# set the solver
solver_nlp = Ipopt.Optimizer

# INPUT ########################################################################
# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_medium_hhc.m")
data = PMs.parse_file(path)

# add the NTs, DTs and necessary nodes


# define the set of considered harmonics
H = collect(1:50)

# replicate data
hdata   = HPM.replicate(data, H=H)

# add shunt impedance at 35kV
Ybase   = 100.0e6 / (35.4e3)^2 
imp     = Dict( "100"   => [0.6396, 0.0085, 2.8057e-6],
                "200"   => [0.2970, 0.0064, 4.1204e-6],
                "500"   => [0.6436, 0.0089, 3.8392e-6],
                "700"   => [0.6100, 0.0081, 3.3892e-6],
                "800"   => [0.7380, 0.0095, 4.1512e-6],
                "1000"  => [Inf, Inf, 2.8692e-06])

for nh in H, (ns,shunt) in hdata["nw"]["$nh"]["shunt"]
    R, L, C = imp[string(shunt["source_id"][2])]
    y = im * 2 * pi * 50 * nh * C / Ybase
    if nh == 1
        shunt["gs"] = 0.0
        shunt["bs"] = 0.0
    else
        shunt["gs"] = real(y)
        shunt["bs"] = imag(y)
    end
end

# RESULTS ######################################################################