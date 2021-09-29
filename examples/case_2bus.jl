using Pkg
Pkg.activate("./")
# load pkgs
using Ipopt, HarmonicPowerModels, PowerModels
using JuMP #adding avoids problems with Revise

# pkg const
const _PMs = PowerModels
const _HPM = HarmonicPowerModels

# set the solver
solver = Ipopt.Optimizer

# path to the data
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_2bus.m")

# load data
data  = _PMs.parse_file(path)
# sanity check
pm_fundamental = _PMs.instantiate_model(data, _PMs.IVRPowerModel, _PMs.build_opf_iv);
result_fundamental = optimize_model!(pm_fundamental, optimizer=solver)
result_fundamental["termination_status"]

# set up current ratings consistently
#assume minimum voltage
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

##
# solve the hopf
# result = run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
pm = _PMs.instantiate_model(hdata, _PMs.IVRPowerModel, _HPM.build_hopf_iv; ref_extensions=[_HPM.ref_add_xfmr!]);
result = optimize_model!(pm, optimizer=solver)

##
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