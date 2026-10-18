################################################################################
# Star test (3- or 5-node): solve the deterministic dHHC and print the bus       #
# voltage phasors + per-bus IHD utilisation for ALL buses.  DATA env selects the #
# case file (default 5node_test.m).                                              #
################################################################################
using HarmonicPowerModels, PowerModels, MosekTools, Ipopt, HSL_jll, JuMP, Printf
const HPM = HarmonicPowerModels; const PMs = PowerModels; PMs.silence()
include(joinpath(@__DIR__, "pce_params.jl"))

DATA = get(ENV, "DATA", "5node_test.m")
mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "ma27",
            "hsllib" => HSL_jll.libhsl_path, "print_level" => 0)

H = [1, 5, 7, 11, 13]; harm = [5, 7, 11, 13]
ihd = HPM.ihd_limits["IEC61000-3-6:2008"]
Eω = PCE_PARAMS["alpha_hat"] / (PCE_PARAMS["alpha_hat"] + PCE_PARAMS["beta_hat"])
Eφ(h) = PCE_PARAMS["harmonics"][h].v_min +
        Eω * (PCE_PARAMS["harmonics"][h].v_max - PCE_PARAMS["harmonics"][h].v_min)

d  = HPM.parse_file(joinpath(@__DIR__, "data", DATA))
hd = HPM.replicate(d, H = H)
for h in harm, (_, l) in hd["nw"]["$h"]["load"]
    l["ref_angle"] = rad2deg(Eφ(h))
end
sol = HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt)
vm(i) = hd["nw"]["1"]["bus"][i]["vm"]
bids = sort(parse.(Int, collect(keys(hd["nw"]["1"]["bus"]))))

@printf("DATA=%s   term=%s   objective (Σ hosting)=%.6e\n", DATA, sol["termination_status"], sol["objective"])
ph(b) = (hypot(b["vr"], b["vi"]), rad2deg(atan(b["vi"], b["vr"])))

println("\nVOLTAGE PHASORS  |V| ∠ deg   [IHD utilisation |V|/(ihd·V₁)]")
println(@sprintf("%-3s", "h") * join((@sprintf(" | %-24s", "bus $i"*(i==first(bids) ? " (gen)" : "")) for i in bids)))
for h in harm
    snw = sol["solution"]["nw"]["$h"]
    cells = map(bids) do i
        m, a = ph(snw["bus"]["$i"]); @sprintf(" | %9.3e∠%+7.2f[%.3f]", m, a, m/(ihd[h]*vm("$i")))
    end
    println(@sprintf("%-3d", h) * join(cells))
end

println("\nHARMONIC HOSTING cmd per unit (h5):")
for k in sort(collect(keys(sol["solution"]["nw"]["5"]["load"])))
    @printf("  load %s: cmd=%.6e\n", k, sol["solution"]["nw"]["5"]["load"][k]["cmd"])
end
