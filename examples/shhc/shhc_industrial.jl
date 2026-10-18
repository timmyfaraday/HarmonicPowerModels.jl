################################################################################
# HarmonicPowerModels.jl — SHHC paper test case: Industrial network (9-bus)    #
#                                                                               #
# Demonstrates, for the 4 LCC converter units, the two equivalences (shared      #
# machinery in shhc_equivalence.jl):                                            #
#   1. DETERMINISTIC — at E[φ], dHHC ≡ sHHC_SOC ≡ sHHC_QCQP in objective,        #
#      injected current phasors and bus voltage phasors.                        #
#   2. STOCHASTIC — sHHC_SOC ≡ sHHC_QCQP for the empirical / Gaussian / Cantelli #
#      chance bounds, cross-checked against a Monte-Carlo dHHC distribution.     #
#                                                                               #
# Run from the MAIN Julia environment (HarmonicPowerModels `Pkg.develop`-ed).    #
# Solvers: Ipopt+MA27 (HPF, sHHC_QCQP), Mosek (dHHC_SOC, sHHC_SOC).              #
# Output: examples/shhc/output/industrial_*.{txt,csv}.  Plots: plot_phasors.jl, #
#         plot_mc_overlay.jl.                                                    #
# Author: Tom Van Acker                                                          #
################################################################################

using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf

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

data = HarmonicPowerModels.parse_file(joinpath(@__DIR__, "data", "industrial_network_shhc.m"))
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
       name = "industrial", outdir = outdir,
       restrict_units = false, qcqp_ok = true,
       nmc = 2000, nsamp_lambda = 20000)

deterministic_equivalence(cfg)
stochastic_equivalence(cfg)

println("\n✓ industrial equivalence done — see $(outdir)/industrial_*.{txt,csv}")
