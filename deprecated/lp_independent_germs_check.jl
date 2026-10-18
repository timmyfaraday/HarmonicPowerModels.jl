################################################################################
# examples/shhc/lp_independent_germs_check.jl                                 #
#                                                                               #
# Numerical verification for Section IV.D of the SHHC paper.                  #
# Run as:  julia --project examples/shhc/lp_independent_germs_check.jl        #
#                                                                               #
# ── What is being checked ───────────────────────────────────────────────────  #
# Under absolute-equality fairness, the stochastic HHC reduces to the closed-  #
# form LP:                                                                      #
#                                                                               #
#   I_h* = min_n  sqrt( V_lim²_{n,h} / C_{n,h} )                             #
#                                                                               #
# where                                                                         #
#   C_{n,h} = α_{n,h,0} + λ sqrt( Σ_{k>0}  ‖ψ_k‖² · α_{n,h,k}² )           #
#                                                                               #
# and the α_{n,h,k} are the PCE coefficients of the squared bus voltage,       #
# obtained by convolving the (complex) voltage-response PCE vectors C_re, C_im #
# via the normalised M tensor.                                                  #
#                                                                               #
# ── Why independent germs matter ────────────────────────────────────────────  #
# When each of the N load units is driven by its OWN independent Beta-          #
# distributed germ, the multivariate PCE lives in a higher-dimensional space:  #
#   independent germs  →  P = C(N+δ, δ) basis functions                       #
#   shared germ        →  P = δ + 1      basis functions                       #
#                                                                               #
# Higher germ dimension introduces more k>0 multi-indices, each contributing   #
# a positive term to the Cantelli σ.  The additional variance terms typically  #
# lower I_h* compared with the shared-germ model, where all units respond to   #
# the same random angle and variance terms partially cancel.                    #
#                                                                               #
# ── Three tables printed ────────────────────────────────────────────────────  #
# TABLE 1  C_{n,h} per non-slack bus for h ∈ {5, 7}.  Assert > 0.             #
# TABLE 2  Optimal I_h* for shared germ vs. independent germs, per harmonic.   #
# TABLE 3  Tightness: α[n*,1]·I*² = (C[n*] − λσ[n*])·I*²  (rel. gap < 1e-10)#
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

using HarmonicPowerModels
using PowerModels
using PrettyTables

const HPM = HarmonicPowerModels
const PMs = PowerModels

# ── parameters ────────────────────────────────────────────────────────────────
const ALPHA_SHAPE = 4.0414    # Beta shape parameter α̂  (Table II)
const BETA_SHAPE  = 4.7762    # Beta shape parameter β̂  (Table II)
const DELTA       = 2         # PCE truncation degree
const EPSILON     = 0.05      # Cantelli risk level ε

# Characteristic 6-pulse harmonics; h=3 skipped (zero-sequence — delta windings
# block current on most buses, making the LP ill-defined without modifications).
const HARMONICS = [5, 7]

# Angle bounds in degrees per harmonic (Table II).
# All four loads share the same marginal Beta distribution.
const ANGLE_DEG = Dict(
    5 => ( 4.75,  40.05),
    7 => (-22.77, 27.51),
)

beta_params = (alpha_shape = ALPHA_SHAPE, beta_shape = BETA_SHAPE, delta = DELTA)

# ── load and parse network data ───────────────────────────────────────────────
data_path = joinpath(HPM.BASE_DIR, "test", "data", "matpower", "industrial_network_hhc.m")
data      = PMs.parse_file(data_path)

# ── extract sorted load info ──────────────────────────────────────────────────
load_entries = sort(
    [(id = parse(Int, k), bus = v["load_bus"]) for (k, v) in data["load"]],
    by = x -> x.id,
)
unit_ids     = [l.id  for l in load_entries]   # load IDs (1,2,3,4)
load_bus_ids = [l.bus for l in load_entries]   # corresponding bus IDs
n_loads      = length(unit_ids)

# ── multivariate PCE workspace: N = n_loads independent germs ────────────────
# (monic Beta01 probability-measure basis, see src/core/pce_basis.jl)
pb_ind = HPM.PCEBasis(n_loads, ALPHA_SHAPE, BETA_SHAPE, DELTA)
P_ind  = pb_ind.dim

lambda = sqrt((1.0 - EPSILON) / EPSILON)

# ─────────────────────────────────────────────────────────────────────────────
# Main loop: one pass per harmonic
# ─────────────────────────────────────────────────────────────────────────────
results = Dict{Int, NamedTuple}()

for h in HARMONICS
    v_min_deg, v_max_deg = ANGLE_DEG[h]
    angle_params_h       = Dict(uid => (v_min_deg, v_max_deg) for uid in unit_ids)

    # ── active bus IDs after Kron reduction (slack + delta-blocked removed) ───
    active_ids = HPM._active_bus_ids_lp(data, h)

    # ── harmonic transfer impedance Z_h ───────────────────────────────────────
    # Z_h[:,u] = voltage at all active buses when load u injects 1+0j A.
    Z_h = HPM.harmonic_transfer_impedance(data, h, load_bus_ids)

    # ── angle PCE: one germ per load (independent germs) ─────────────────────
    # Unit at position u in the vector occupies germ dimension u.
    units_ind = [(id      = unit_ids[u],
                  h_angles= Dict(h => (v_min_deg, v_max_deg)))
                 for u in 1:n_loads]
    cos_c_ind, sin_c_ind = HPM.angle_pce_coefficients(pb_ind, units_ind, [h])

    # ── voltage-response PCE coefficient matrices ─────────────────────────────
    # C_re[n,k] = Σ_u [ Re(Z[n,u])·cos_k[u] − Im(Z[n,u])·sin_k[u] ]
    # C_im[n,k] = Σ_u [ Im(Z[n,u])·cos_k[u] + Re(Z[n,u])·sin_k[u] ]
    C_re, C_im = HPM.compute_C_coefficients(
        Z_h, cos_c_ind, sin_c_ind, unit_ids, pb_ind, h,
    )

    # ── PCE coefficients of the squared voltage magnitude ─────────────────────
    # alpha[n,k] = Σ_{κ1,κ2} M[κ1,κ2,k] (C_re[n,κ1]C_re[n,κ2] + C_im[n,κ1]C_im[n,κ2])
    alpha = HPM.compute_alpha(C_re, C_im, pb_ind)

    # ── Cantelli scalar and LP solution ───────────────────────────────────────
    C_scal             = HPM.compute_C_scalar(alpha, pb_ind, EPSILON)
    V_lim_sq           = HPM.bus_voltage_limits(data, h)
    I_ind, binding_ind = HPM.lp_solution(C_scal, V_lim_sq, active_ids)

    # ── shared-germ case for comparison ───────────────────────────────────────
    comp = HPM.compare_shared_vs_independent(
        data, h, EPSILON, beta_params, angle_params_h,
    )
    I_sh       = comp.I_star_shared
    binding_sh = comp.binding_bus_shared

    # ── tightness check at the binding bus ────────────────────────────────────
    # By construction: C_scal[n*] = alpha[n*,1] + lambda * sigma[n*]
    # Hence at the LP optimum I*:
    #   lhs = alpha[n*,1] * I*²   (k=0 term of |V|², = E[|V|²/I²]·I*²)
    #   rhs = (C_scal[n*] − lambda·σ) * I*²
    # The relative gap quantifies floating-point round-off in the pipeline.
    # Probability-measure basis: σ = √(Σ_{k≥2} ⟨Ψ_k,Ψ_k⟩·alpha[k]²), no
    # total-measure correction.
    bidx  = findfirst(==(binding_ind), active_ids)
    sigma = sqrt(sum(pb_ind.sp2[k] * alpha[bidx, k]^2 for k in 2:P_ind; init=0.0))
    lhs_t = alpha[bidx, 1] * I_ind^2              # E[|V|²] = alpha[1] * I²
    rhs_t = (C_scal[bidx] - lambda * sigma) * I_ind^2
    gap_t = abs(lhs_t - rhs_t) / max(abs(lhs_t), 1e-16)

    results[h] = (
        active_ids  = active_ids,
        C_scal      = C_scal,
        I_ind       = I_ind,
        I_sh        = I_sh,
        binding_ind = binding_ind,
        binding_sh  = binding_sh,
        lhs_t       = lhs_t,
        rhs_t       = rhs_t,
        gap_t       = gap_t,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# TABLE 1 — Cantelli C scalars per non-slack bus
# ─────────────────────────────────────────────────────────────────────────────
println()
act5 = results[5].active_ids
act7 = results[7].active_ids

# For positive-sequence (h=5) and negative-sequence (h=7) there is no zero-seq
# delta blocking, so the active bus sets are identical.
n_active = length(act5)
t1       = Matrix{Any}(undef, n_active, 3)
for n in 1:n_active
    t1[n, 1] = act5[n]
    t1[n, 2] = results[5].C_scal[n]
    t1[n, 3] = n <= length(act7) ? results[7].C_scal[n] : NaN
end
pretty_table(t1;
    title      = "TABLE 1  Cantelli C scalars  (independent germs, N=$(n_loads), ε=$(EPSILON))",
    column_labels = ["Bus ID", "C  h=5", "C  h=7"],
    formatters = [fmt__printf("%.6g", [2, 3])],
)

all_pos = all(all(results[h].C_scal .> 0) for h in HARMONICS)
println("All C scalars > 0: $all_pos\n")

# ─────────────────────────────────────────────────────────────────────────────
# TABLE 2 — LP solution: shared germ vs. independent germs
# ─────────────────────────────────────────────────────────────────────────────
t2 = Matrix{Any}(undef, length(HARMONICS), 4)
for (i, h) in enumerate(HARMONICS)
    t2[i, 1] = h
    t2[i, 2] = results[h].I_sh
    t2[i, 3] = results[h].I_ind
    t2[i, 4] = results[h].binding_ind
end
pretty_table(t2;
    title      = "TABLE 2  LP solution  (P_ind = $(P_ind) = C($(n_loads)+$(DELTA),$(DELTA))  vs.  P_sh = $(DELTA+1))",
    column_labels = ["h", "I*_shared [pu]", "I*_independent [pu]", "binding bus (indep.)"],
    formatters = [fmt__printf("%.8g", [2, 3])],
)

# ─────────────────────────────────────────────────────────────────────────────
# TABLE 3 — k=0 lifting tightness at the LP optimum
# ─────────────────────────────────────────────────────────────────────────────
t3 = Matrix{Any}(undef, length(HARMONICS), 5)
for (i, h) in enumerate(HARMONICS)
    r      = results[h]
    t3[i, 1] = h
    t3[i, 2] = r.binding_ind
    t3[i, 3] = r.lhs_t
    t3[i, 4] = r.rhs_t
    t3[i, 5] = r.gap_t
end
pretty_table(t3;
    title      = "TABLE 3  k=0 lifting tightness at LP optimum",
    column_labels = ["h", "binding bus", "α[n*,1]·I*²", "(C[n*]−λσ)·I*²", "rel. gap"],
    formatters = [fmt__printf("%.4e", [3, 4, 5])],
)

# ─────────────────────────────────────────────────────────────────────────────
# Assertions
# ─────────────────────────────────────────────────────────────────────────────
println("── Assertions ──────────────────────────────────────────────────────────")
for h in HARMONICS
    r = results[h]
    @assert all(r.C_scal .> 0)   "C_scalar not positive at h=$h"
    @assert isfinite(r.I_ind)     "I*_independent not finite at h=$h"
    @assert isfinite(r.I_sh)      "I*_shared not finite at h=$h"
    @assert r.gap_t < 1e-10       "k=0 tightness gap $(r.gap_t) ≥ 1e-10 at h=$h"
end
println("All assertions passed.")
