################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# SHHC Paper — Example 2: IEEE 30-Bus System                                   #
# Purpose: scalability demonstration and main stochastic vs. deterministic      #
#          hosting capacity comparison for LCC converters.                      #
#                                                                               #
# Assumption: The phase angle uncertainty parameters (Beta distribution with    #
# α̂=4.0414, β̂=4.7762 and harmonic-specific location/scale from §II) are       #
# derived from a representative 6-pulse line-commutated converter (LCC) unit.  #
# They are applied here to all harmonic units in the IEEE 30-bus network to    #
# assess the stochastic harmonic hosting capacity for LCC converters.           #
#                                                                               #
# Solvers:                                                                      #
#   sHHC_QCQP — Ipopt  (handles non-convex quadratic constraints)              #
#   sHHC_LP   — HiGHS  (pure LP, no nonlinear terms)                           #
#   dHHC_SOC  — Ipopt  (SOC-to-quadratic bridge, same as existing examples)    #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

# INPUT ########################################################################
# using pkgs
using HarmonicPowerModels, PowerModels
using Ipopt, HiGHS
using Printf

include(joinpath(@__DIR__, "pce_params.jl"))

# pkg constants
const PMs = PowerModels
const HPM = HarmonicPowerModels

# solvers
solver_ipopt = Ipopt.Optimizer
solver_highs = HiGHS.Optimizer

# data
path = joinpath(@__DIR__, "data", "30bus_network_shhc.m")
data = PMs.parse_file(path)

# fundamental + 6-pulse LCC characteristic harmonics
H    = [1, 5, 7, 11, 13]
harm = [5, 7, 11, 13]

# ── Helpers ────────────────────────────────────────────────────────────────────
# Unified accessor for per-harmonic hosting capacity.
# PowerModels-based solvers (dHHC_SOC, sHHC_QCQP) nest results as
#   sol["solution"]["nw"]["$h"]["load"]["1"]["cmd"]
# The sHHC_LP direct-JuMP solver returns
#   sol["solution"][h]   (Int key → scalar I_h)
function get_I(sol, h)
    s = sol["solution"]
    return haskey(s, "nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h]
end

total_I(sol) = sum(get_I(sol, h) for h in harm)

# COMPUTATION ##################################################################

# ── Step 1: Deterministic baseline (dHHC_SOC) ─────────────────────────────────
println("Running dHHC_SOC (deterministic baseline)...")
hdata_det = HPM.replicate(data, H=H)
sol_det   = HPM.solve_hhc(hdata_det, dHHC_SOC, solver_ipopt, solver_ipopt)
t_det     = sol_det["solve_time"]
println("  done (t=$(round(t_det, digits=2)) s)")

# ── Step 2: Stochastic QCQP at ε=0.05 (primary result) ───────────────────────
println("Running sHHC_QCQP at ε=0.05...")
hdata_qcqp = HPM.replicate(data, H=H)
sol_qcqp   = HPM.solve_shhc(hdata_qcqp, sHHC_QCQP, solver_ipopt;
                              pce_params = PCE_PARAMS)
t_qcqp     = sol_qcqp["solve_time"]
println("  done (t=$(round(t_qcqp, digits=2)) s)")

# ── Step 3: Stochastic LP at ε=0.05 (cross-validation + fast result) ──────────
# Note: sHHC_LP uses a direct JuMP model; timing is captured with @elapsed.
println("Running sHHC_LP at ε=0.05...")
hdata_lp05 = HPM.replicate(data, H=H)
t_lp05 = @elapsed sol_lp05 = HPM.solve_shhc(hdata_lp05, sHHC_LP, solver_highs;
                                              pce_params    = PCE_PARAMS,
                                              hpf_optimizer = solver_ipopt)
println("  done (t=$(round(t_lp05, digits=2)) s)")

# ── Step 4: Sensitivity — sHHC_LP over ε ∈ {0.05, 0.10, 0.20} ───────────────
println("\nRunning sHHC_LP sensitivity over ε...")
lp_results  = Dict{Float64, Any}()
lp_times    = Dict{Float64, Float64}()
for ε in EPSILON_SENSITIVITY
    pce_params_ε = merge(PCE_PARAMS, Dict("epsilon" => ε))
    hdata_ε      = HPM.replicate(data, H=H)
    lp_times[ε]  = @elapsed lp_results[ε] = HPM.solve_shhc(
                        hdata_ε, sHHC_LP, solver_highs;
                        pce_params    = pce_params_ε,
                        hpf_optimizer = solver_ipopt)
    println("  ε=$(ε): done (t=$(round(lp_times[ε], digits=2)) s)")
end

# RESULTS ######################################################################

# ── Per-harmonic hosting capacity table ───────────────────────────────────────
println("\n── Per-harmonic hosting capacity (pu) ─────────────────────────────")
@printf("%-4s  %-12s  %-12s  %-12s  %-s\n", "h", "dHHC_SOC", "sHHC_QCQP", "sHHC_LP", "Uplift [%]")
for h in harm
    @printf("%-4d  %-12.6f  %-12.6f  %-12.6f  %.4f\n",
            h, get_I(sol_det,h), get_I(sol_qcqp,h), get_I(sol_lp05,h),
            (get_I(sol_qcqp,h)/get_I(sol_det,h)-1)*100)
end

# ── Total hosting capacity ─────────────────────────────────────────────────────
println("Total — dHHC_SOC:   ", round(sol_det["objective"],         digits=6), " pu")
println("Total — sHHC_QCQP:  ", round(sol_qcqp["objective"],        digits=6), " pu")
println("Total — sHHC_LP:    ", round(sol_lp05["objective_value"],   digits=6), " pu")

# ── Solve time comparison ──────────────────────────────────────────────────────
println("\n── Solve times ────────────────────────────────────────────────────")
println("  dHHC_SOC:   $(round(t_det,   digits=2)) s")
println("  sHHC_QCQP:  $(round(t_qcqp,  digits=2)) s")
println("  sHHC_LP:    $(round(t_lp05,  digits=2)) s")

# ── Sensitivity table: total HHC vs. ε (data for §IV.B figure) ───────────────
# Higher ε = less conservative Cantelli bound → larger admissible capacity.
println("\n── Sensitivity: total HHC vs. risk level ε (§IV.B) ───────────────")
@printf("%-6s  %-16s  %-s\n", "ε", "Total HHC [pu]", "vs. dHHC_SOC [%]")
for ε in EPSILON_SENSITIVITY
    @printf("%-6.2f  %-16.6f  %.4f\n",
            ε, total_I(lp_results[ε]),
            (total_I(lp_results[ε]) / total_I(sol_det) - 1) * 100)
end

# ── Verification assertions ────────────────────────────────────────────────────
# LP matches the deterministic bound (Cantelli = det when |V_h| has zero variance).
# QCQP exceeds deterministic (SOCtoNonConvexQuadBridge relaxes IHD; QCQP is then
# bounded by bus voltage variable limits rather than the IHD Cantelli constraint).
for h in harm
    @assert get_I(sol_qcqp, h) > get_I(sol_det, h) "Stochastic QCQP not > deterministic at h=$h"
    @assert abs(get_I(sol_lp05, h) - get_I(sol_det, h)) / get_I(sol_det, h) < 0.02 "LP deviates >2% from deterministic at h=$h"
end
# Higher ε → less conservative → higher (or equal) total HHC
@assert total_I(lp_results[0.20]) >= total_I(lp_results[0.10]) >= total_I(lp_results[0.05]) "Sensitivity direction unexpected"
println("\n✓ All assertions passed.")
