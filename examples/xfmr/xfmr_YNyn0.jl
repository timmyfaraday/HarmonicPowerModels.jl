# load pkgs
using Dierckx
using Ipopt
using HarmonicPowerModels
using Plots
using PowerModels

# pkg const
const PMs = PowerModels
const HPM = HarmonicPowerModels

# path to the data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/xfmr/case_xfmr_YNyn0.m")

# load data
data = PMs.parse_file(path)

# BH-curve
B⁺ = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
H⁺ = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
B = vcat(reverse(-B⁺),0.0,B⁺)
H = vcat(reverse(-H⁺),0.0,H⁺) 
BH_powercore_h100_23 = Dierckx.Spline1D(B, H; k=3, bc="nearest")

# xfmr magnetizing data
magn = Dict("Hᴱ"    => [1,3], 
            "Hᴵ"    => collect(1:2:25),
            "Fᴱ"    => :rectangular,
            "Fᴵ"    => :rectangular,
            "xfmr"  => Dict(1 => Dict(  "l"  => 10.0,
                                        "A"  => 0.5,
                                        "N"  => 75,
                                        "BH" => BH_powercore_h100_23)
                            )
            )

# harmonic data
hdata = HPM.replicate(data, xfmr_magn=magn)

# set the solver
solver = Ipopt.Optimizer

# solve the hopf
model  = PMs.instantiate_model(hdata, PMs.IVRPowerModel, HPM.build_hopf_iv; ref_extensions=[HPM.ref_add_xfmr!]);
result = optimize_model!(model, optimizer=solver, solution_processors=[HPM.sol_data_model!])

HPM.append_indicators!(result, hdata)

for (n,nw) in result["solution"]["nw"]
    for (i,bus) in nw["bus"]
            bus["vm"] =  abs(bus["vr"] +im*  bus["vi"])
            bus["va"] =  angle(bus["vr"] +im*  bus["vi"])*180/pi
    end
end
println("Harmonic 3")
_PMs.print_summary(result["solution"]["nw"]["2"])
println("Harmonic 1")
_PMs.print_summary(result["solution"]["nw"]["1"])
result["objective"]
result["termination_status"]

##
# Yy -> third harmonic of the excitation cannot flow -> current is inherently forced to 0
# this means the third harmonic voltage is 0
# 