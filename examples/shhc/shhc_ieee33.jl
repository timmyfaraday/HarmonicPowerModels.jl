################################################################################
# HarmonicPowerModels.jl — SHHC recombination test: IEEE 33-bus radial feeder   #
#                                                                               #
# Tests whether a RADIAL distribution feeder (single substation bottleneck,     #
# many loads aggregating upstream) shows more harmonic recombination than the    #
# meshed IEEE-30 transmission case.                                             #
#                                                                               #
# With all 32 loads as independent-germ units, the degree-2 PCE basis is 561    #
# terms (the dimensionality wall the germ-clustering step will address), so this #
# script uses the PCE-FREE path: the Monte-Carlo dHHC distribution (independent  #
# per-unit angles, re-solved deterministically) + dHHC@E[φ].  The MC spread      #
# (q95/q05) and the recombination gain (q05 vs dHHC@E[φ]) quantify recombination.#
# Output: examples/shhc/output/ieee33_stoch_*.{txt,csv}.                         #
# Author: Tom Van Acker                                                          #
################################################################################

using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Statistics

include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "shhc_equivalence.jl"))

const TLIM = 360.0
ipopt_ma27 = optimizer_with_attributes(Ipopt.Optimizer,
                "linear_solver" => "ma27", "hsllib" => HSL_jll.libhsl_path,
                "print_level" => 0, "max_wall_time" => TLIM,
                "tol" => 1e-10, "constr_viol_tol" => 1e-10,
                "acceptable_tol" => 1e-10, "acceptable_constr_viol_tol" => 1e-10)
mosek_soc  = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true,
                "MSK_DPAR_OPTIMIZER_MAX_TIME" => TLIM)

data = HarmonicPowerModels.parse_file(joinpath(@__DIR__, "data", "33bus_network_shhc.m"))
H    = [1, 5, 7, 11, 13]
harm = [5, 7, 11, 13]
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)

cfg = (; data, H, harm,
       units         = sort(parse.(Int, collect(keys(data["load"])))),
       base_ranges   = Dict(h => (PCE_PARAMS["harmonics"][h].v_min,
                                  PCE_PARAMS["harmonics"][h].v_max) for h in harm),
       α̂ = PCE_PARAMS["alpha_hat"], β̂ = PCE_PARAMS["beta_hat"], ε = PCE_PARAMS["epsilon"],
       pce_params    = PCE_PARAMS,
       ipopt = ipopt_ma27, mosek = mosek_soc,
       name = "ieee33", outdir = outdir,
       restrict_units = false, qcqp_ok = false,
       nmc = 2000, nsamp_lambda = 20000)

println("IEEE-33 radial recombination test: $(length(cfg.units)) units (PCE-free MC)")
r = stochastic_equivalence(cfg; bounds = String[])

o = r.objmc; q(p) = quantile(o, p)
HE = sum(get_I(solve_dhhc_pm(cfg)[1], h) for h in harm)
@printf("\n%s\nIEEE-33 SUMMARY (N=%d radial)\n%s\n", repeat("=",72), length(cfg.units), repeat("=",72))
@printf("dHHC@E[φ] = %.5e   (MC pct %.1f%%)\n", HE, 100*count(<(HE),o)/length(o))
@printf("MC  q05=%.5e  median=%.5e  q95=%.5e\n", q(.05), q(.5), q(.95))
@printf("MC spread q95/q05         = %.4f\n", q(.95)/q(.05))
@printf("recombination gain (q05 vs dHHC@E[φ]) = %+.2f%%\n", 100*(q(.05)-HE)/HE)
@printf("%s\n", repeat("=",72))
