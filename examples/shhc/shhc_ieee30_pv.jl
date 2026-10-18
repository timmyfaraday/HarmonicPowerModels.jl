################################################################################
# HarmonicPowerModels.jl — SHHC recombination study: IEEE 30-bus, single slack  #
#                                                                               #
# Baseline for the recombination levers (feature/shhc-recombination, Phase 0)   #
# on data/30bus_network_shhc_pv.m: bus 1 is the only fundamental reference, the #
# other generator buses are PV, no bus is a clean (ideal harmonic ground) bus,  #
# and the generator short-circuit data follow the underground-cables paper.     #
#                                                                               #
# Environment variables:                                                         #
#   SHHC_N         number of harmonic units (first N loads), default 8           #
#   SHHC_ANGLE_REF "global" (default) or "local": harmonic angles relative to   #
#                  the network reference, or to each unit's own fundamental bus  #
#                  voltage (rotated by h·θ_i)                                    #
#   SHHC_DET       "1" to also run the deterministic equivalence (needs QCQP)    #
#   SHHC_NMC       Monte-Carlo samples, default 1000                             #
#   SHHC_CONIC     conic solver: "mosek" (default) or "gurobi"                   #
#                                                                               #
# Solvers: Mosek (default) or Gurobi for the conic hosting-capacity models,     #
# Ipopt (MUMPS) for the fundamental power flow and the QCQP — Ipopt with ma27   #
# fails on this power flow. Gurobi returns LOCALLY_SOLVED points that overstate #
# sHHC_SOC by 1.5-6.7% against Mosek at N=21, so Mosek is the reference.        #
# Output: examples/shhc/output/ieee30pv_N<N>_<ref>_*.{txt,csv}.                 #
# Author: Tom Van Acker                                                          #
################################################################################

using HarmonicPowerModels, PowerModels
using Ipopt, Gurobi, MosekTools
using JuMP, Printf

include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "shhc_equivalence.jl"))

const TLIM  = 360.0
const NUNIT = parse(Int, get(ENV, "SHHC_N", "8"))
const AREF  = get(ENV, "SHHC_ANGLE_REF", "global")
const DET   = get(ENV, "SHHC_DET", "0") == "1"
const NMC   = parse(Int, get(ENV, "SHHC_NMC", "1000"))
const CONIC = get(ENV, "SHHC_CONIC", "mosek")

ipopt  = optimizer_with_attributes(Ipopt.Optimizer,
            "print_level" => 0, "max_wall_time" => TLIM,
            "tol" => 1e-10, "constr_viol_tol" => 1e-10,
            "acceptable_tol" => 1e-10, "acceptable_constr_viol_tol" => 1e-10)
const GRB_ENV = Gurobi.Env()
# Gurobi's default barrier hits numerical trouble on dHHC_SOC for this case and
# its default cone tolerance leaves sHHC_SOC ~0.5% off Mosek; these settings
# solve both and agree with Mosek to ~1e-4.
gurobi = optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => 0, "TimeLimit" => TLIM,
            "Presolve" => 0, "BarHomogeneous" => 1, "NumericFocus" => 3,
            "BarQCPConvTol" => 1e-10)
mosek  = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true,
            "MSK_DPAR_OPTIMIZER_MAX_TIME" => TLIM)
conic  = CONIC == "gurobi" ? gurobi : mosek

PowerModels.silence()
data = HPM.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc_pv.m"))
H    = [1, 5, 7, 11, 13]
harm = [5, 7, 11, 13]
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)
all_units = sort(parse.(Int, collect(keys(data["load"]))))

cfg = (; data, H, harm,
       units           = all_units[1:NUNIT],
       base_ranges     = Dict(h => (PCE_PARAMS["harmonics"][h].v_min,
                                    PCE_PARAMS["harmonics"][h].v_max) for h in harm),
       α̂ = PCE_PARAMS["alpha_hat"], β̂ = PCE_PARAMS["beta_hat"], ε = PCE_PARAMS["epsilon"],
       pce_params      = PCE_PARAMS,
       ipopt = ipopt, mosek = conic,           # `mosek` is the conic-solver slot
       name = "ieee30pv_N$(NUNIT)_$(AREF)" * (CONIC == "mosek" ? "" : "_$(CONIC)"), outdir = outdir,
       angle_reference = AREF,
       restrict_units = true, qcqp_ok = NUNIT ≤ 6,
       nmc = NMC, nsamp_lambda = 20000)

println("IEEE-30 single slack, N=$NUNIT units $(cfg.units), angle reference: $AREF")
DET && deterministic_equivalence(cfg)
stochastic_equivalence(cfg)

println("\n✓ done — see $(outdir)/$(cfg.name)_*")
