using Pkg
Pkg.activate("./")
# load pkgs
using Ipopt, HarmonicPowerModels, PowerModels
using JuMP #avoids problems with Revise

# pkg const
const _PMs = PowerModels
const _HPM = HarmonicPowerModels

# path to the data
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_basf.m")

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


            # xfmr = Dict("voltage_harmonics" => [1,3],
            # "current_harmonics" => [1,3,5,7,9,11,13],
            # "N" => 50,
            # "current_type" => :rectangular,
            # "excitation_type" => :sigmoid,
            # "inom" => 0.4,
            # "ψmax" => 0.5,
            # "voltage_type" => :rectangular,
            # "dv" => [0.1,0.1],
            # "vmin" => [0.0,0.0],
            # "vmax" => [1.1,1.1],
            # "dθ" => [π/5,π/5],
            # "θmin" => [0.0,0.0],
            # "θmax" => [2π,2π])

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
result = run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)