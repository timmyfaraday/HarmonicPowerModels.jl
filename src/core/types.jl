################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
################################################################################

"""
Deterministic Harmonic Hosting Capacity (NLP)
"""
mutable struct dHHC_NLP <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
Deterministic Harmonic Hosting Capacity (SOC)
"""
mutable struct dHHC_SOC <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
    sHHC_QCQP <: _PMs.AbstractIVRModel

Stochastic Harmonic Hosting Capacity solved as a non-convex Quadratically
Constrained Quadratic Program (QCQP).

Each harmonic unit's phase angle is modelled as a Beta-distributed random
variable expanded in a degree-2 Jacobi polynomial chaos expansion (PCE).
The resulting network equations couple PCE coefficients through the
multiplication tensor, producing quadratic constraints:

- Lifting (product) constraints for the k=0 PCE index are second-order cone
  (SOC) inequalities that bound the 0-th moment of squared current magnitudes.
- Lifting constraints for k≠0 PCE indices are quadratic equalities linking
  higher-order coefficients.
- Chance constraints on bus-voltage distortion and line-current magnitude are
  enforced as SOC inequalities derived via the one-sided Cantelli bound:
      P(X > μ + λ·σ) ≤ ε  ⟺  μ + λ(ε)·σ ≤ limit,  λ(ε) = √((1-ε)/ε)

Recommended solvers: Ipopt (non-convex NLP mode) or Gurobi (non-convex QCQP).
"""
mutable struct sHHC_QCQP <: _PMs.AbstractIVRModel _PMs.@pm_fields end

"""
    sHHC_SOC <: _PMs.AbstractIVRModel

Stochastic Harmonic Hosting Capacity solved as a Second-Order Cone Program (SOCP).

Applicable under the absolute-equality fairness assumption: all harmonic units
are assigned an equal scalar hosting capacity I_h. The squared bus voltages and
branch currents are then quadratic in I_h (computed analytically via the harmonic
admittance matrix Y_h and the PCE multiplication tensor), so the chance
constraints reduce to a small program in the per-harmonic variables I_h:

- The individual-harmonic-distortion (IHD) chance constraint is per harmonic and
  factorises into a LINEAR bound on I_h².
- The total-harmonic-distortion (THD), RMS-voltage (per bus) and RMS-current
  (per branch) chance constraints couple the harmonics quadratically and become
  SECOND-ORDER CONE constraints in I_h².

Each I_h is additionally capped by the smallest `load["c_rating"]` among the
harmonic units, mirroring the `cmd ≤ c_rating` bound of the dHHC and QCQP models:
under absolute equality the common allocation cannot exceed what the weakest unit
can carry.

Because the single scalar I_h is what makes the formulation conic, absolute
equality is the ONLY fairness principle this model can represent. `solve_shhc`
rejects any other `hdata["principle"]` rather than ignoring it; use `sHHC_QCQP`,
which carries a per-unit current variable.

Both the Cantelli (one-sided Chebyshev) and the Gaussian chance-constraint
bounds are supported via `pce_params["bound"]` (see `update_hdata_with_pce_data!`).
Requires a conic (SOCP) solver — e.g. Mosek or Clarabel (HiGHS cannot handle the
second-order cones).
"""
mutable struct sHHC_SOC <: _PMs.AbstractIVRModel _PMs.@pm_fields end