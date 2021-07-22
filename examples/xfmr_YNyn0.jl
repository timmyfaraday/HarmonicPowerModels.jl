# load pkgs
using Ipopt, HarmonicPowerModels, PowerModels

# pkg const
const _PMs = PowerModels
const _HPM = HarmonicPowerModels

# path to the data
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_xfmr_YNyn0.m")

# transformer excitation data
xfmr = Dict("voltage_harmonics" => [1,3],
            "current_harmonics" => [1,3,5,7,9,11,13],
            "N" => 50,
            "current_type" => :rectangular,
            "excitation_type" => :sigmoid,
            "inom" => 0.4,
            "ψmax" => 0.5,
            "voltage_type" => :rectangular,
            "dv" => [0.1,0.05],
            "vmin" => [-0.1,-0.1],
            "vmax" => [1.1,0.3],
            "dθ" => [π/5,π/5],
            "θmin" => [0.0,0.0],
            "θmax" => [2π,2π])

# load data
data  = _PMs.parse_file(path)
hdata = _HPM.replicate(data, xfmr_exc=xfmr)

# set the solver
solver = Ipopt.Optimizer

# solve the hopf
run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)