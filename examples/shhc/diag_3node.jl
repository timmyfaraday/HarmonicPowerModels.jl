################################################################################
# 3-node symmetric test: solve the deterministic dHHC and print the bus voltage  #
# phasors + per-bus IHD utilisation, to see WHICH bus binds as a function of the  #
# branch / generator impedances.                                                 #
#                                                                               #
# Set the impedances below:                                                       #
#   CFG = :unit  → large branches, stiff gen  ⇒ unit buses (2,3) bind            #
#   CFG = :gen   → small branches, weak  gen  ⇒ generator bus (1) binds          #
################################################################################
using HarmonicPowerModels, PowerModels, MosekTools, Ipopt, HSL_jll, JuMP, Printf
const HPM = HarmonicPowerModels; const PMs = PowerModels; PMs.silence()
include(joinpath(@__DIR__, "pce_params.jl"))

CFG = Symbol(get(ENV, "CFG", "gen"))
# (branch r, branch x, gen rsc, gen xsc)
br_r, br_x, g_r, g_x = CFG === :gen ? (0.0002, 0.002, 0.10, 2.0) :   # gen-bus binding
                                       (0.01,   0.10,  0.001, 0.02)   # unit-bus binding

mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "ma27",
            "hsllib" => HSL_jll.libhsl_path, "print_level" => 0)

H = [1, 5, 7, 11, 13]; harm = [5, 7, 11, 13]
ihd = HPM.ihd_limits["IEC61000-3-6:2008"]              # per-harmonic IHD limit
Eω = PCE_PARAMS["alpha_hat"] / (PCE_PARAMS["alpha_hat"] + PCE_PARAMS["beta_hat"])
Eφ(h) = PCE_PARAMS["harmonics"][h].v_min +
        Eω * (PCE_PARAMS["harmonics"][h].v_max - PCE_PARAMS["harmonics"][h].v_min)

d = HPM.parse_file(joinpath(@__DIR__, "data", "3node_test.m"))
for (_, br) in d["branch"]; br["br_r"] = br_r; br["br_x"] = br_x; end
for (_, g)  in d["gen"];    g["rsc"]   = g_r;  g["xsc"]   = g_x;  end

hd = HPM.replicate(d, H = H)
for h in harm, (_, l) in hd["nw"]["$h"]["load"]
    l["ref_angle"] = rad2deg(Eφ(h))
end
sol = HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt)
vm(i) = hd["nw"]["1"]["bus"][i]["vm"]                  # fundamental voltage (HPF result)

@printf("CFG=%s  (branch %.4f+j%.4f, gen_sc %.3f+j%.3f)\n", CFG, br_r, br_x, g_r, g_x)
@printf("term=%s   objective (Σ hosting)=%.6e\n", sol["termination_status"], sol["objective"])
@printf("fundamental |V|: bus1=%.4f bus2=%.4f bus3=%.4f\n", vm("1"), vm("2"), vm("3"))
ph(b) = (hypot(b["vr"], b["vi"]), rad2deg(atan(b["vi"], b["vr"])))

println("\nVOLTAGE PHASORS  |V| ∠ deg     and    IHD utilisation |V|/(ihd·V₁) in [ ]")
@printf("%-3s | %-26s %-26s %-26s\n", "h", "bus 1 (gen)", "bus 2 (unit)", "bus 3 (unit)")
for h in harm
    nw = "$h"; snw = sol["solution"]["nw"][nw]
    cell(i) = (m = ph(snw["bus"][i])[1]; a = ph(snw["bus"][i])[2];
               u = m / (ihd[h] * vm(i)); @sprintf("%9.3e∠%+7.2f [%.3f]", m, a, u))
    @printf("%-3d | %-26s %-26s %-26s\n", h, cell("1"), cell("2"), cell("3"))
end
