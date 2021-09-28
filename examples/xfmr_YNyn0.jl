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

# hdata = _HPM.replicate(data)
hdata = _HPM.replicate(data, xfmr_exc=xfmr)

# set the solver
solver = Ipopt.Optimizer

# solve the hopf
# run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
pm = _PMs.instantiate_model(hdata, _PMs.IVRPowerModel, _HPM.build_hopf_iv; ref_extensions=[_HPM.ref_add_xfmr!]);
result = optimize_model!(pm, optimizer=solver)


##
print(pm.model)

# Checking Kirchhoff

# Remarks for F:
# 1) scaling of the harmonic load seems to be inconsistent 

cd = [  result["solution"]["nw"]["1"]["load"]["1"]["crd"],
        result["solution"]["nw"]["1"]["load"]["1"]["cid"],
        result["solution"]["nw"]["2"]["load"]["1"]["crd"],
        result["solution"]["nw"]["2"]["load"]["1"]["cid"]]
ct_to =[result["solution"]["nw"]["1"]["xfmr"]["1"]["crt_to"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["cit_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["crt_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["cit_to"]]
cst_to =[result["solution"]["nw"]["1"]["xfmr"]["1"]["csrt_to"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["csit_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csrt_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csit_to"]]
cet =[  result["solution"]["nw"]["1"]["xfmr"]["1"]["cert"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["ceit"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["cert"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["ceit"]]
cst_fr =[result["solution"]["nw"]["1"]["xfmr"]["1"]["csrt_fr"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["csit_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csrt_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csit_fr"]]
ct_fr =[result["solution"]["nw"]["1"]["xfmr"]["1"]["crt_fr"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["cit_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["crt_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["cit_fr"]]
cg =[   result["solution"]["nw"]["1"]["gen"]["1"]["crg"],
        result["solution"]["nw"]["1"]["gen"]["1"]["cig"],
        result["solution"]["nw"]["2"]["gen"]["1"]["crg"],
        result["solution"]["nw"]["2"]["gen"]["1"]["cig"]]


# enforced by constraint_current_balance
# 1_crt[(1, 2, 1)] + 1_crd[1] == 0.0
# 1_cit[(1, 2, 1)] + 1_cid[1] == 0.0
# 2_crt[(1, 2, 1)] + 2_crd[1] == 0.0
# 2_cit[(1, 2, 1)] + 2_cid[1] == 0.0

# TESTS
@assert isapprox(cd[3], 0.092 * cd[1], rtol=1e-6)
@assert isapprox(cd[4], 0.092 * cd[2], rtol=1e-6)
@assert all(isapprox.(cd, -ct_to, rtol=1e-6))
@assert all(isapprox.(ct_to, cst_to, rtol=1e-6))
@assert all(isapprox.(cst_fr, -cst_to, rtol=1e-6))
@assert all(isapprox.(ct_fr, cst_fr, rtol=1e-6))
@assert all(isapprox.(ct_fr, cg, rtol=1e-6))

### Kirchhoff without shunts works and without excitation

vd = [  result["solution"]["nw"]["1"]["bus"]["2"]["vr"],
        result["solution"]["nw"]["1"]["bus"]["2"]["vi"],
        result["solution"]["nw"]["2"]["bus"]["2"]["vr"],
        result["solution"]["nw"]["2"]["bus"]["2"]["vi"]]
vt_to =[result["solution"]["nw"]["1"]["xfmr"]["1"]["vrt_to"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["vit_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["vrt_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["vit_to"]]
cst_to =[result["solution"]["nw"]["1"]["xfmr"]["1"]["csrt_to"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["csit_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csrt_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csit_to"]]
cet =[  result["solution"]["nw"]["1"]["xfmr"]["1"]["cert"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["ceit"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["cert"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["ceit"]]
cst_fr =[result["solution"]["nw"]["1"]["xfmr"]["1"]["csrt_fr"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["csit_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csrt_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["csit_fr"]]
ct_fr =[result["solution"]["nw"]["1"]["xfmr"]["1"]["crt_fr"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["cit_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["crt_fr"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["cit_fr"]]
cg =[   result["solution"]["nw"]["1"]["gen"]["1"]["crg"],
        result["solution"]["nw"]["1"]["gen"]["1"]["cig"],
        result["solution"]["nw"]["2"]["gen"]["1"]["crg"],
        result["solution"]["nw"]["2"]["gen"]["1"]["cig"]]

        vd = [  result["solution"]["nw"]["1"]["bus"]["2"]["vr"],
        result["solution"]["nw"]["1"]["bus"]["2"]["vi"],
        result["solution"]["nw"]["2"]["bus"]["2"]["vr"],
        result["solution"]["nw"]["2"]["bus"]["2"]["vi"]]
vt_to =[result["solution"]["nw"]["1"]["xfmr"]["1"]["vrt_to"],
        result["solution"]["nw"]["1"]["xfmr"]["1"]["vit_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["vrt_to"],
        result["solution"]["nw"]["2"]["xfmr"]["1"]["vit_to"]]


@assert all(isapprox.(vd .- vt_to, 0.01 .* ct_to, rtol=1e-6))