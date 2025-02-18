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
using CSV
using DataFrames
using Ipopt 
using Plots, StatsPlots
using PrettyTables
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

# define the set of considered harmonics
H = collect(1:50)

# replicate data
hdata   = HPM.replicate(data, H=H)

# add shunt impedance [(R,L) // C]
Ybase   = 100.0e6 / (35.4e3)^2 
imp     = Dict( "100"   => [0.6396, 0.0085, 2.8057e-6],
                "200"   => [0.2970, 0.0064, 4.1204e-6],
                "500"   => [0.6436, 0.0089, 3.8392e-6],
                "700"   => [0.6100, 0.0081, 3.3892e-6],
                "800"   => [0.7380, 0.0095, 4.1512e-6],
                "1000"  => [Inf, Inf, 2.8692e-06])

for nh in H, (ns,shunt) in hdata["nw"]["$nh"]["shunt"]
    R, L, C = imp[string(shunt["source_id"][2])]
    y = (1 / (sqrt(nh) * R + im * 2 * pi * 50 * nh * L) + (im * 2 * pi * 50 * nh * C)) / Ybase
    if nh == 1
        shunt["gs"] = 0.0
        shunt["bs"] = 0.0
    else
        shunt["gs"] = real(y)
        shunt["bs"] = imag(y)
    end
end

# RESULTS ######################################################################
# solve HHC problem
results = HPM.solve_hhc(hdata, dHHC_SOC, solver_nlp, solver_nlp)

# IEC61000-3-6 RESULTS ######################################################### 
# read-in data
path = joinpath(HPM.BASE_DIR,"test/data/excel/exposia_post_cs_base_harmonic_impedance_2023_07_20_init.csv")
df   = DataFrame(CSV.File(path))
df.Frequency .= df.Frequency ./ 50.0

dU   = 35.4e3 .* HPM.ihd_limits["IEC61000-2-4:2002, Cl. 2"]

Ilim = [dU[nh] / first(filter(x -> x.Frequency == nh, df).BASF_R100) for nh in setdiff(H,1)]

# VISUALISATION ################################################################
# FIGURE: 
Ibase = 100.0e6 / sqrt(3) / 35.4e3
Im(nh,nl) = Ibase * round(results["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"], digits=7)

bar(string.(50.0 .* setdiff(H,1)), Ilim ./ [Im(h,"1") for h in setdiff(H,1)], yaxis=:log, ylim=(0.1,1e3))