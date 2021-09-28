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
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_3bus2gen.m")

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
result["termination_status"]
sol = result["solution"]
pd1 = hdata["nw"]["1"]["load"]["1"]["pd"]
qd1 = hdata["nw"]["1"]["load"]["1"]["qd"]

pd1 = sol["nw"]["1"]["load"]["1"]["pd"]
qd1 = sol["nw"]["1"]["load"]["1"]["qd"]

vmsb1   = abs(sol["nw"]["1"]["bus"]["1"]["vr"] + im* sol["nw"]["1"]["bus"]["1"]["vi"])
vmload1 = abs(sol["nw"]["1"]["bus"]["2"]["vr"] + im* sol["nw"]["1"]["bus"]["2"]["vi"])
vmgen1 = abs(sol["nw"]["1"]["bus"]["3"]["vr"] + im* sol["nw"]["1"]["bus"]["3"]["vi"])

cmload1 = hypot(sol["nw"]["1"]["load"]["1"]["crd"], sol["nw"]["1"]["load"]["1"]["cid"])
cmgen1 = hypot(sol["nw"]["1"]["gen"]["1"]["crg"], sol["nw"]["1"]["gen"]["1"]["cig"])
cmgen2 = hypot(sol["nw"]["1"]["gen"]["2"]["crg"], sol["nw"]["1"]["gen"]["2"]["cig"])

vmsb2 = abs(sol["nw"]["2"]["bus"]["1"]["vr"] + im*sol["nw"]["2"]["bus"]["1"]["vi"] )
vmload2 = abs(sol["nw"]["2"]["bus"]["2"]["vr"] + im* sol["nw"]["2"]["bus"]["2"]["vi"])
cmgen1_2 = hypot(sol["nw"]["2"]["gen"]["1"]["crg"], sol["nw"]["2"]["gen"]["1"]["cig"])
cmgen2_2 = hypot(sol["nw"]["2"]["gen"]["2"]["crg"], sol["nw"]["2"]["gen"]["2"]["cig"])

cmload2 = hypot(sol["nw"]["2"]["load"]["1"]["crd"], sol["nw"]["2"]["load"]["1"]["cid"])

multiplier = hdata["nw"]["2"]["load"]["1"]["multiplier"]
cmload2/cmload1

##
# print(pm.model)

# result