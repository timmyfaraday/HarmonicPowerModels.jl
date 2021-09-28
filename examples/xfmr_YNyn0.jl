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

hdata = _HPM.replicate(data, xfmr_exc=xfmr)

# set the solver
solver = Ipopt.Optimizer

# solve the hopf
# run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
pm = _PMs.instantiate_model(hdata, _PMs.IVRPowerModel, _HPM.build_hopf_iv; ref_extensions=[_HPM.ref_add_xfmr!]);
result = optimize_model!(pm, optimizer=solver)


##
print(pm.model)

result["termination_status"]
sol = result["solution"]
pd1 = hdata["nw"]["1"]["load"]["1"]["pd"]
qd1 = hdata["nw"]["1"]["load"]["1"]["qd"]

pd1 = sol["nw"]["1"]["load"]["1"]["pd"]
qd1 = sol["nw"]["1"]["load"]["1"]["qd"]
ccmd1 = sol["nw"]["1"]["load"]["1"]["ccmd"]

vmsb1 = sol["nw"]["1"]["bus"]["1"]["vm"]
vmload1 = sol["nw"]["1"]["bus"]["2"]["vm"]

cmxfmr1 = hypot(sol["nw"]["1"]["xfmr"]["1"]["crt_fr"], sol["nw"]["1"]["xfmr"]["1"]["cit_fr"])


vmsb2 = sol["nw"]["2"]["bus"]["1"]["vm"]
vmload2 = sol["nw"]["2"]["bus"]["2"]["vm"]

cmbranch2 = hypot(sol["nw"]["2"]["xfmr"]["1"]["crt_fr"], sol["nw"]["2"]["xfmr"]["1"]["cit_fr"])


multiplier = hdata["nw"]["2"]["load"]["1"]["multiplier"]
cmbranch2/cmbranch1

pd2 = sol["nw"]["2"]["load"]["1"]["pd"]
qd2 = sol["nw"]["2"]["load"]["1"]["qd"]
ccmd2 = sol["nw"]["2"]["load"]["1"]["ccmd"]