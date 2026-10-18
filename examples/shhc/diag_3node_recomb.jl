################################################################################
# 3-node, 2-unit recombination test.  Independent random angles per unit (MC).   #
# Compares two impedance regimes:                                                 #
#   :unit → unit buses bind (self-dominated)   ⇒ expect ~no recombination         #
#   :gen  → generator bus binds (aggregating)  ⇒ expect recombination             #
# Reports the deterministic dHHC@E[φ] (coherent) baseline and the MC distribution #
# of per-unit hosting; recombination gain = (q05 − dHHC@E[φ]) / dHHC@E[φ].        #
################################################################################
using HarmonicPowerModels, PowerModels, MosekTools, Ipopt, HSL_jll
using JuMP, Printf, Statistics, Distributions, Random
const HPM = HarmonicPowerModels; const PMs = PowerModels; PMs.silence()
include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "shhc_equivalence.jl"))

mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "ma27",
            "hsllib" => HSL_jll.libhsl_path, "print_level" => 0)
H = [1, 5, 7, 11, 13]; harm = [5, 7, 11, 13]

function mc_test(br_r, br_x, g_r, g_x; nmc = 2000)
    d = HPM.parse_file(joinpath(@__DIR__, "data", "3node_test.m"))
    for (_, br) in d["branch"]; br["br_r"] = br_r; br["br_x"] = br_x; end
    for (_, g)  in d["gen"];    g["rsc"]   = g_r;  g["xsc"]   = g_x;  end
    cfg = (; data = d, H, harm, units = sort(parse.(Int, collect(keys(d["load"])))),
           base_ranges = Dict(h => (PCE_PARAMS["harmonics"][h].v_min,
                                    PCE_PARAMS["harmonics"][h].v_max) for h in harm),
           α̂ = PCE_PARAMS["alpha_hat"], β̂ = PCE_PARAMS["beta_hat"], ε = PCE_PARAMS["epsilon"],
           pce_params = PCE_PARAMS, ipopt = ipopt, mosek = mosek,
           name = "t", outdir = tempdir(), restrict_units = false, qcqp_ok = false,
           nmc = nmc, nsamp_lambda = 100)
    HE = sum(get_I(solve_dhhc_pm(cfg)[1], h) for h in harm)         # coherent dHHC@E[φ]
    Random.seed!(1); Bd = Beta(cfg.α̂, cfg.β̂); o = Float64[]
    for m in 1:nmc
        ω   = Dict(dd => rand(Bd) for dd in cfg.units)              # independent germ per unit
        rng = Dict(h => Dict(dd => (cfg.base_ranges[h][1] +
                    (cfg.base_ranges[h][2]-cfg.base_ranges[h][1])*ω[dd],) for dd in cfg.units) for h in harm)
        s, _ = solve_dhhc(cfg, rng); push!(o, sum(get_I(s, h) for h in harm))
    end
    (HE = HE, q05 = quantile(o,.05), med = quantile(o,.5), q95 = quantile(o,.95),
     spread = quantile(o,.95)/quantile(o,.05), gain = 100*(quantile(o,.05)-HE)/HE)
end

@printf("%-26s %-11s %-11s %-11s %-11s %-9s %-9s\n",
        "config", "dHHC@E[φ]", "q05", "median", "q95", "gain@q05", "gain@med")
for (lbl, p) in (("unit-binding (self)", (0.01, 0.10, 0.001, 0.02)),
                 ("gen-binding (aggregate)", (0.0002, 0.002, 0.10, 2.0)))
    r = mc_test(p...)
    @printf("%-26s %-11.4e %-11.4e %-11.4e %-11.4e %+7.2f%% %+7.2f%%\n",
            lbl, r.HE, r.q05, r.med, r.q95, r.gain, 100*(r.med-r.HE)/r.HE)
end
