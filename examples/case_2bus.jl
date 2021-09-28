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
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_line.m")

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
ccmd1 = sol["nw"]["1"]["load"]["1"]["ccmd"]

vmsb1 = sol["nw"]["1"]["bus"]["1"]["vm"]
vmload1 = sol["nw"]["1"]["bus"]["2"]["vm"]

zbranch1 = hdata["nw"]["1"]["branch"]["1"]["br_r"] + im* hdata["nw"]["1"]["branch"]["1"]["br_x"]
cmbranch1 = hypot(sol["nw"]["1"]["branch"]["1"]["cr_fr"], sol["nw"]["1"]["branch"]["1"]["ci_fr"])


vmsb2 = sol["nw"]["2"]["bus"]["1"]["vm"]
vmload2 = sol["nw"]["2"]["bus"]["2"]["vm"]

zbranch2 = hdata["nw"]["2"]["branch"]["1"]["br_r"] + im* hdata["nw"]["2"]["branch"]["1"]["br_x"]
cmbranch2 = hypot(sol["nw"]["2"]["branch"]["1"]["cr_fr"], sol["nw"]["2"]["branch"]["1"]["ci_fr"])


multiplier = hdata["nw"]["2"]["load"]["1"]["multiplier"]
cmbranch2/cmbranch1

pd2 = sol["nw"]["2"]["load"]["1"]["pd"]
qd2 = sol["nw"]["2"]["load"]["1"]["qd"]
ccmd2 = sol["nw"]["2"]["load"]["1"]["ccmd"]

##
# print(pm.model)

result