# load pkgs
using Ipopt, HarmonicPowerModels, PowerModels

# pkg const
const _PMs = PowerModels
const _HPM = HarmonicPowerModels

# path to the data
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_xfmr_YNyn0_simplified.m")

# transformer excitation data
xfmr = Dict("voltage_harmonics" => [1,3],
            "current_harmonics" => [1,3],
            "N" => 50,
            "current_type" => :rectangular,
            "excitation_type" => :sigmoid,
            "inom" => 0.4,
            "ψmax" => 0.5,
            "voltage_type" => :rectangular,
            "dv" => [0.1,0.1],
            "vmin" => [0.0,0.0],
            "vmax" => [1.1,1.1],
            "dθ" => [π/5,π/5],
            "θmin" => [0.0,0.0],
            "θmax" => [2π,2π])

# load data
data  = _PMs.parse_file(path)

vmmin = 0.8
for (b,branch) in data["branch"]
    branch["c_rating"] = branch["rate_a"]/vmmin
end
for (d,load) in data["load"]
    load["c_rating"] = abs(load["pd"] + im* load["qd"])/vmmin
end
for (g,gen) in data["gen"]
    gen["c_rating"] = abs(gen["pmax"] + im* gen["qmax"])/vmmin
end

hdata = _HPM.replicate(data)
# hdata = _HPM.replicate(data, xfmr_exc=xfmr)

# set the solver
solver = Ipopt.Optimizer

# solve the hopf
# run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
pm = _PMs.instantiate_model(hdata, _PMs.IVRPowerModel, _HPM.build_hopf_iv; ref_extensions=[_HPM.ref_add_xfmr!]);
result = optimize_model!(pm, optimizer=solver, solution_processors=[ _HPM.sol_data_model!])
_HPM.append_indicators!(result, hdata)


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