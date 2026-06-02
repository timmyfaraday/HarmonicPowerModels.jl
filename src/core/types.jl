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
    sHHC_LP <: _PMs.AbstractIVRModel

Stochastic Harmonic Hosting Capacity solved as a Linear Program (LP).

Applicable under the absolute-equality fairness assumption: all harmonic
units are assigned an equal scalar hosting capacity I_h. Under this
assumption the harmonic current injections are linear in I_h, and the PCE
coefficients of every bus voltage and branch current are linear in I_h
(computed analytically via the harmonic admittance matrix Y_h). Consequently:

- PCE coefficients of all voltages and currents are linear in I_h.
- The Cantelli-based chance constraints on voltage distortion and current
  magnitude reduce to linear inequality constraints in I_h.
- No lifting (product) constraints are required.

The resulting problem is a pure LP and can be solved with any LP solver
(e.g., HiGHS, GLPK, Gurobi).
"""
mutable struct sHHC_LP <: _PMs.AbstractIVRModel _PMs.@pm_fields end