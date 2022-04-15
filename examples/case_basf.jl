# using Pkg
# Pkg.activate("./")
# load pkgs
using Ipopt, HarmonicPowerModels, PowerModels
using JuMP #avoids problems with Revise
using Dierckx

# pkg const
const _PMs = PowerModels
const _HPM = HarmonicPowerModels

# path to the data
# path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_basf_simplified_no_filter.m")
path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_basf_simplified_with_filter.m")
# path = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_basf.m")

# transformer excitation data
# exc_1 = Dict("voltage_harmonics" => [1,5],
#             "current_harmonics" => [1,3,5,7,9,13],
#             "N" => 50,
#             "current_type" => :rectangular,
#             "excitation_type" => :sigmoid,
#             "inom" => 0.13,
#             "ψmax" => 1,
#             "voltage_type" => :rectangular,
#             "dv" => [0.1,0.1],
#             "vmin" => [-1.1,-1.1],
#             "vmax" => [1.1,1.1],
#             "dθ" => [π/5,π/5],
#             "θmin" => [0.0,0.0],
#             "θmax" => [2π,2π])

# exc_2 = Dict("voltage_harmonics" => [1,3],
#             "current_harmonics" => [1,3],
#             "N" => 50,
#             "current_type" => :rectangular,
#             "excitation_type" => :sigmoid,
#             "inom" => 0.4,
#             "ψmax" => 0.5,
#             "voltage_type" => :rectangular,
#             "dv" => [0.1,0.1],
#             "vmin" => [-1.1,-1.1],
#             "vmax" => [1.1,1.1],
#             "dθ" => [π/5,π/5],
#             "θmin" => [0.0,0.0],
#             "θmax" => [2π,2π])

# xfmr_exc = Dict("1" => exc_1, "2" => exc_2)

# BH-curve
B⁺ = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
H⁺ = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
B = vcat(reverse(-B⁺),0.0,B⁺)
H = vcat(reverse(-H⁺),0.0,H⁺) 
BH_powercore_h100_23 = Dierckx.Spline1D(B, H; k=3, bc="nearest")

# xfmr magnetizing data
magn = Dict("Hᴱ"    => [1, 5], 
            "Hᴵ"    => collect(1:2:19),
            "Fᴱ"    => :rectangular,
            "Fᴵ"    => :rectangular,
            "Emax"  => 1.1,
            "IDH"   => [1.0, 0.06],
            "pcs"   => [6, 6],
            "xfmr"  => Dict(1 => Dict(  "l"     => 11.4,
                                        "A"     => 0.5,
                                        "N"     => 500,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 150000),
                                        2 => Dict(  "l"     => 11.4,
                                        "A"     => 0.5,
                                        "N"     => 500,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 150000),
                                        3 => Dict(  "l"     => 11.4,
                                        "A"     => 0.5,
                                        "N"     => 500,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 150000),
                                        4 => Dict(  "l"     => 11.4,
                                        "A"     => 0.5,
                                        "N"     => 500,
                                        "BH"    => BH_powercore_h100_23,
                                        "Vbase" => 150000)
                            )
            )

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

# hdata = _HPM.replicate(data, xfmr_exc=exc_1)
hdata = _HPM.replicate(data, xfmr_magn=magn)

for (n, nw) in hdata["nw"]
    for (t, xfmr) in nw["xfmr"]
        if t != "1"
            # @show xfmr
            xfmr["voltage_harmonics"] = String[]
            xfmr["current_harmonics"] = String[]
            xfmr["voltage_harmonics_ntws"] = []
            xfmr["current_harmonics_ntws"] = []
        end
    end
end

# set the solver
solver = Ipopt.Optimizer


#solve power flow
resultpf = run_hpf_iv(hdata, _PMs.IVRPowerModel, solver)
# @assert resultpf["termination_status"] == LOCALLY_SOLVED
# _HPM.append_indicators!(resultpf, hdata)

# println("Harmonic 13")
# _PMs.print_summary(resultpf["solution"]["nw"]["6"])
# println("Harmonic 9")
# _PMs.print_summary(resultpf["solution"]["nw"]["5"])
# println("Harmonic 7")
# _PMs.print_summary(resultpf["solution"]["nw"]["4"])
# println("Harmonic 5")
# _PMs.print_summary(resultpf["solution"]["nw"]["3"])
# vm = resultpf["solution"]["nw"]["3"]["bus"]["1"]

# println("Harmonic 3")
# _PMs.print_summary(resultpf["solution"]["nw"]["2"])
# vm = resultpf["solution"]["nw"]["2"]["bus"]["6"]
# println("Harmonic 1")
# _PMs.print_summary(resultpf["solution"]["nw"]["1"])

# pg = resultpf["solution"]["nw"]["1"]["gen"]["1"]


##
# solve the hopf
result = run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
pm = _PMs.instantiate_model(hdata, _PMs.IVRPowerModel, _HPM.build_hopf_iv; ref_extensions=[_HPM.ref_add_xfmr!]);
result = optimize_model!(pm, optimizer=solver, solution_processors=[ _HPM.sol_data_model!])
@assert result["termination_status"] == LOCALLY_SOLVED
_HPM.append_indicators!(result, hdata)


pg = result["solution"]["nw"]["1"]["gen"]["1"]
println("Harmonic 13")
_PMs.print_summary(result["solution"]["nw"]["6"])
println("Harmonic 9")
_PMs.print_summary(result["solution"]["nw"]["5"])
println("Harmonic 7")
_PMs.print_summary(result["solution"]["nw"]["4"])
println("Harmonic 5")
_PMs.print_summary(result["solution"]["nw"]["3"])
println("Harmonic 3")
_PMs.print_summary(result["solution"]["nw"]["2"])
println("Harmonic 1")
_PMs.print_summary(result["solution"]["nw"]["1"])
result["objective"]
result["termination_status"]


Dict(n=>(nw["gen"]["2"]["crg"]+ im*nw["gen"]["2"]["cig"] ) for (n,nw) in result["solution"]["nw"])
result["solution"]["nw"]["1"]["gen"]["1"]["pg"] + im* result["solution"]["nw"]["1"]["gen"]["1"]["qg"]
# @show result["solution"]["nw"]["1"]["xfmr"]["1"]["pexc"]


#
hdata["nw"]["1"]["xfmr"]["1"]["NWᴱ"]
hdata["nw"]["2"]["xfmr"]["1"]["NWᴱ"]
hdata["nw"]["3"]["xfmr"]["1"]["NWᴱ"]
hdata["nw"]["4"]["xfmr"]["1"]["NWᴱ"]
hdata["nw"]["5"]["xfmr"]["1"]["NWᴱ"]
