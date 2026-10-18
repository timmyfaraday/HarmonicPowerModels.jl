################################################################################
# HarmonicPowerModels.jl — SHHC paper test case: IEEE 30-bus system            #
#                                                                               #
# Same two equivalences as the industrial case (shared machinery in             #
# shhc_equivalence.jl), at a representative number of harmonic units N where the #
# QCQP is solvable within the time limit:                                       #
#   1. DETERMINISTIC — dHHC ≡ sHHC_SOC ≡ sHHC_QCQP at E[φ] (obj + phasors).      #
#   2. STOCHASTIC — sHHC_SOC ≡ sHHC_QCQP for empirical/Gaussian/Cantelli + MC.   #
#                                                                               #
# Run from the MAIN Julia environment. Set SHHC_N to change the unit count.      #
# Output: examples/shhc/output/ieee30_*.{txt,csv}.                              #
# Author: Tom Van Acker                                                          #
################################################################################

using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf

include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "shhc_equivalence.jl"))

const TLIM = 360.0
const NUNIT = parse(Int, get(ENV, "SHHC_N", "4"))   # # harmonic units (QCQP solvable at N≈4)
ipopt_ma27 = optimizer_with_attributes(Ipopt.Optimizer,
                "linear_solver" => "ma27", "hsllib" => HSL_jll.libhsl_path,
                "print_level" => 0, "max_wall_time" => TLIM,
                "tol" => 1e-10, "constr_viol_tol" => 1e-10,
                "acceptable_tol" => 1e-10, "acceptable_constr_viol_tol" => 1e-10)
mosek_soc  = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true,
                "MSK_DPAR_OPTIMIZER_MAX_TIME" => TLIM)

data = PowerModels.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc.m"))
H    = [1, 5, 7, 11, 13]
harm = [5, 7, 11, 13]
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)
all_units = sort(parse.(Int, collect(keys(data["load"]))))

cfg = (; data, H, harm,
       units         = all_units[1:NUNIT],
       base_ranges   = Dict(h => (PCE_PARAMS["harmonics"][h].v_min,
                                  PCE_PARAMS["harmonics"][h].v_max) for h in harm),
       α̂ = PCE_PARAMS["alpha_hat"], β̂ = PCE_PARAMS["beta_hat"], ε = PCE_PARAMS["epsilon"],
       pce_params    = PCE_PARAMS,
       ipopt = ipopt_ma27, mosek = mosek_soc,
       name = "ieee30", outdir = outdir,
       restrict_units = true, qcqp_ok = NUNIT ≤ 6,   # QCQP within time limit only at small N
       nmc = 2000, nsamp_lambda = 20000)

println("IEEE-30 equivalence at N=$NUNIT units: $(cfg.units)")
deterministic_equivalence(cfg)
stochastic_equivalence(cfg)

println("\n✓ ieee30 equivalence done — see $(outdir)/ieee30_*.{txt,csv}")
