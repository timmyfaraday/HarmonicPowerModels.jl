################################################################################
# HarmonicPowerModels.jl                                                       #
# PCE constraint functions for sHHC_QCQP (stochastic HHC, non-convex QCQP).  #
#                                                                               #
# Every deterministic IV constraint becomes one linear or quadratic constraint  #
# per PCE index k ∈ K = 0:δ. Lifted squared-magnitude variables (w_pce,        #
# j_pce) are linked to the voltage/current PCE variables via the degree-δ       #
# multiplication tensor M[j,l,k] = ⟨ψ_j ψ_l, ψ_k⟩/⟨ψ_k,ψ_k⟩ (1-indexed).   #
# Cantelli SOC chance constraints enforce distortion limits with risk ε.        #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

# ── Internal PCE data accessors ────────────────────────────────────────────────

_pce_M(pm)           = pm.data["pce"]["M"]
_pce_lambda(pm)      = pm.data["pce"]["lambda"]
_pce_norms(pm)       = pm.data["pce"]["psi_norms"]
_pce_cos_phi(pm,h,d) = pm.data["pce"]["cos_phi"][h][d]
_pce_sin_phi(pm,h,d) = pm.data["pce"]["sin_phi"][h][d]

# ── helpers for constraint_current_balance_pce ─────────────────────────────────
# Returns transformer arc currents summed over `arcs` for PCE index k, or 0
# when no transformer variables have been created for this network.
function _xfmr_cr_pce(pm, nw, arcs, k)
    isempty(arcs) && return 0.0
    crx_pce = _PMs.var(pm, nw, :crx_pce)
    return sum(crx_pce[t, k] for t in arcs)
end
function _xfmr_ci_pce(pm, nw, arcs, k)
    isempty(arcs) && return 0.0
    cix_pce = _PMs.var(pm, nw, :cix_pce)
    return sum(cix_pce[t, k] for t in arcs)
end

# ══════════════════════════════════════════════════════════════════════════════
# 1. Reference bus
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_ref_bus_voltage_pce(pm::sHHC_QCQP, i; nw, K)

Fix all PCE coefficients of the reference bus voltage.

For harmonic network `nw` (nw ≠ fundamental), the grid injects no harmonic
voltage at the reference bus (deterministic zero source). Therefore:

    vr_pce[i, k] = 0,  vi_pce[i, k] = 0   ∀ k ∈ K

If called for the fundamental network, the k=0 coefficient is set to 1.0
(per-unit nominal) and all k≥1 coefficients remain zero, since the fundamental
voltage is treated as deterministic.
"""
function constraint_ref_bus_voltage_pce(pm::sHHC_QCQP, i::Int;
                                         nw::Int=fundamental(pm), K=0:2)
    vr_pce = _PMs.var(pm, nw, :vr_pce)
    vi_pce = _PMs.var(pm, nw, :vi_pce)

    v_ref_re = (nw == fundamental(pm)) ? 1.0 : 0.0

    JuMP.@constraint(pm.model, vr_pce[i, 0] == v_ref_re)
    JuMP.@constraint(pm.model, vi_pce[i, 0] == 0.0)
    for k in K
        k == 0 && continue
        JuMP.@constraint(pm.model, vr_pce[i, k] == 0.0)
        JuMP.@constraint(pm.model, vi_pce[i, k] == 0.0)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 2. Harmonic unit injection current
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_load_current_pce(pm::sHHC_QCQP, d; nw, K)

PCE expansion of the harmonic unit injection current (one pair per k ∈ K).

Harmonic unit d injects I_d = cmd[d] · exp(j φ_{d,h}) where φ_{d,h} is
Beta-distributed. The PCE coefficients of cos(φ) and sin(φ) are precomputed
parameters stored in `pm.data["pce"]["cos_phi"][nw][d]`. Taking k-th index:

    crd_pce[d, k] = cos_phi[k+1] · cmd[d]
    cid_pce[d, k] = sin_phi[k+1] · cmd[d]

These are linear constraints: cmd[d] is the deterministic scalar hosting
capacity variable and the trig PCE coefficients are constant parameters.
"""
function constraint_load_current_pce(pm::sHHC_QCQP, d::Int;
                                      nw::Int=fundamental(pm), K=0:2)
    crd_pce = _PMs.var(pm, nw, :crd_pce)
    cid_pce = _PMs.var(pm, nw, :cid_pce)
    cmd     = _PMs.var(pm, nw, :cmd, d)

    cos_phi = _pce_cos_phi(pm, nw, d)   # length-|K| vector, 1-indexed
    sin_phi = _pce_sin_phi(pm, nw, d)

    for k in K
        JuMP.@constraint(pm.model, crd_pce[d, k] == cos_phi[k+1] * cmd)
        JuMP.@constraint(pm.model, cid_pce[d, k] == sin_phi[k+1] * cmd)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 3. Kirchhoff current law
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_current_balance_pce(pm::sHHC_QCQP, i; nw, K)

PCE Kirchhoff current law at bus i (one linear constraint pair per k ∈ K).

Because network parameters (shunt conductance gs, susceptance bs) are
deterministic, the PCE coefficients decouple across k and satisfy independent
linear balance equations:

    Σ_a cr_pce[a,k] = −Σ_d crd_pce[d,k] − gs_total·vr_pce[i,k] + bs_total·vi_pce[i,k]
    Σ_a ci_pce[a,k] = −Σ_d cid_pce[d,k] − gs_total·vi_pce[i,k] − bs_total·vr_pce[i,k]

Transformer and filter currents are omitted in the sHHC_QCQP model; add
crx_pce / crf_pce terms here when those components are included.
"""
function constraint_current_balance_pce(pm::sHHC_QCQP, i::Int;
                                         nw::Int=fundamental(pm), K=0:2)
    bus_arcs      = _PMs.ref(pm, nw, :bus_arcs,      i)
    bus_arcs_xfmr = _PMs.ref(pm, nw, :bus_arcs_xfmr, i)
    bus_loads     = _PMs.ref(pm, nw, :bus_loads,      i)
    bus_shunts    = _PMs.ref(pm, nw, :bus_shunts,     i)

    bus_gs = Dict(s => _PMs.ref(pm, nw, :shunt, s, "gs") for s in bus_shunts)
    bus_bs = Dict(s => _PMs.ref(pm, nw, :shunt, s, "bs") for s in bus_shunts)

    vr_pce  = _PMs.var(pm, nw, :vr_pce)
    vi_pce  = _PMs.var(pm, nw, :vi_pce)
    cr_pce  = _PMs.var(pm, nw, :cr_pce)
    ci_pce  = _PMs.var(pm, nw, :ci_pce)
    crd_pce = _PMs.var(pm, nw, :crd_pce)
    cid_pce = _PMs.var(pm, nw, :cid_pce)

    gs_total = isempty(bus_shunts) ? 0.0 : sum(values(bus_gs))
    bs_total = isempty(bus_shunts) ? 0.0 : sum(values(bus_bs))

    for k in K
        JuMP.@constraint(pm.model,
            sum(cr_pce[a, k] for a in bus_arcs)
            + _xfmr_cr_pce(pm, nw, bus_arcs_xfmr, k)
            ==
            - sum(crd_pce[d, k] for d in bus_loads)
            - gs_total * vr_pce[i, k]
            + bs_total * vi_pce[i, k]
        )
        JuMP.@constraint(pm.model,
            sum(ci_pce[a, k] for a in bus_arcs)
            + _xfmr_ci_pce(pm, nw, bus_arcs_xfmr, k)
            ==
            - sum(cid_pce[d, k] for d in bus_loads)
            - gs_total * vi_pce[i, k]
            - bs_total * vr_pce[i, k]
        )
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 4. Ohm's law / voltage drop
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_voltage_drop_pce(pm::sHHC_QCQP, b; nw, K)

PCE voltage drop across branch b (one linear constraint pair per k ∈ K).

Because branch resistance r and reactance x are deterministic parameters, the
PCE coefficient of index k satisfies the same linear Ohm's law structure as
the deterministic model:

    vr_pce[t_bus, k] = vr_pce[f_bus, k] − r·cr_pce[f_idx, k] + x·ci_pce[f_idx, k]
    vi_pce[t_bus, k] = vi_pce[f_bus, k] − r·ci_pce[f_idx, k] − x·cr_pce[f_idx, k]

Simplified π-model (tap ratio tm = 1, no branch shunt admittance), consistent
with the HPM harmonic network convention where branch shunts are zero.
"""
function constraint_voltage_drop_pce(pm::sHHC_QCQP, b::Int;
                                      nw::Int=fundamental(pm), K=0:2)
    branch = _PMs.ref(pm, nw, :branch, b)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (b, f_bus, t_bus)
    r      = branch["br_r"]
    x      = branch["br_x"]

    vr_pce = _PMs.var(pm, nw, :vr_pce)
    vi_pce = _PMs.var(pm, nw, :vi_pce)
    cr_pce = _PMs.var(pm, nw, :cr_pce)
    ci_pce = _PMs.var(pm, nw, :ci_pce)

    for k in K
        JuMP.@constraint(pm.model,
            vr_pce[t_bus, k] == vr_pce[f_bus, k]
                                 - r * cr_pce[f_idx, k]
                                 + x * ci_pce[f_idx, k]
        )
        JuMP.@constraint(pm.model,
            vi_pce[t_bus, k] == vi_pce[f_bus, k]
                                 - r * ci_pce[f_idx, k]
                                 - x * cr_pce[f_idx, k]
        )
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 5. Lifted voltage auxiliary constraints
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_lifted_voltage_pce(pm::sHHC_QCQP, i; nw, K)

Link the lifted variable w_pce[i,k] (PCE coefficient of |V_{h,i}|²) to the
voltage PCE coefficients via the multiplication tensor M.

The k-th PCE coefficient of |V|² = Vr² + Vi² is:
    (|V|²)_k = Σ_{j,l ∈ K} M[j+1,l+1,k+1] · (vr_pce[i,j]·vr_pce[i,l]
                                               + vi_pce[i,j]·vi_pce[i,l])

k=0 (convex quadratic inequality, equivalent to SOC):
    Σ_{j,l} M[j+1,l+1,1] · (vr_j·vr_l + vi_j·vi_l) ≤ w_pce[i,0]
    Since M[:,:,0] is diagonal, this simplifies to Σ_j M[j,j,0]·(vr_j²+vi_j²) ≤ w_pce[i,0].
    The inequality is tight at the optimal solution of a maximisation problem.

k≥1 (non-convex quadratic equality, responsible for PCE variance accuracy):
    Σ_{j,l} M[j+1,l+1,k+1] · (vr_j·vr_l + vi_j·vi_l) = w_pce[i,k]
    Off-diagonal M entries produce bilinear terms, making this non-convex.
"""
function constraint_lifted_voltage_pce(pm::sHHC_QCQP, i::Int;
                                        nw::Int=fundamental(pm), K=0:2)
    M      = _pce_M(pm)
    vr_pce = _PMs.var(pm, nw, :vr_pce)
    vi_pce = _PMs.var(pm, nw, :vi_pce)
    w_pce  = _PMs.var(pm, nw, :w_pce)

    for k in K
        lhs = sum(
            M[j+1, l+1, k+1] * (vr_pce[i, j] * vr_pce[i, l]
                                + vi_pce[i, j] * vi_pce[i, l])
            for j in K, l in K
        )
        if k == 0
            JuMP.@constraint(pm.model, lhs <= w_pce[i, 0])
        else
            JuMP.@constraint(pm.model, lhs == w_pce[i, k])
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 6. Lifted current auxiliary constraints
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_lifted_current_pce(pm::sHHC_QCQP, b; nw, K)

Link the lifted variable j_pce[f_idx,k] (PCE coefficient of |I_{b,h}|²) to
the branch from-arc current PCE coefficients via the multiplication tensor M.

Mirrors `constraint_lifted_voltage_pce` with cr_pce/ci_pce in place of
vr_pce/vi_pce and j_pce in place of w_pce:

    Σ_{j,l ∈ K} M[j+1,l+1,k+1] · (cr_pce[f_idx,j]·cr_pce[f_idx,l]
                                   + ci_pce[f_idx,j]·ci_pce[f_idx,l])
        ≤ j_pce[f_idx, 0]   (k=0, convex inequality used in Cantelli SOC)
        = j_pce[f_idx, k]   (k≥1, quadratic equality for variance correctness)
"""
function constraint_lifted_current_pce(pm::sHHC_QCQP, b::Int;
                                        nw::Int=fundamental(pm), K=0:2)
    M      = _pce_M(pm)
    branch = _PMs.ref(pm, nw, :branch, b)
    f_bus  = branch["f_bus"]
    t_bus  = branch["t_bus"]
    f_idx  = (b, f_bus, t_bus)

    cr_pce = _PMs.var(pm, nw, :cr_pce)
    ci_pce = _PMs.var(pm, nw, :ci_pce)
    j_pce  = _PMs.var(pm, nw, :j_pce)

    for k in K
        lhs = sum(
            M[j+1, l+1, k+1] * (cr_pce[f_idx, j] * cr_pce[f_idx, l]
                                + ci_pce[f_idx, j] * ci_pce[f_idx, l])
            for j in K, l in K
        )
        if k == 0
            JuMP.@constraint(pm.model, lhs <= j_pce[f_idx, 0])
        else
            JuMP.@constraint(pm.model, lhs == j_pce[f_idx, k])
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# 7. IHD voltage chance constraint
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_voltage_ihd_chance(pm::sHHC_QCQP, i; nw)

One-sided Cantelli chance constraint on the individual harmonic distortion
(IHD) of the bus voltage at bus i for harmonic network nw.

P(|V_{h,i}|² > limit) ≤ ε is enforced via the Cantelli bound:
    E[|V_{h,i}|²] + λ · √Var[|V_{h,i}|²] ≤ (ihdmax · v_fund)²

where λ = √((1-ε)/ε) and, in terms of the lifted PCE variable w_pce:
    E[|V|²]   = w_pce[i, 0]
    Var[|V|²] = Σ_{k≥1} psi_norms[k+1] · w_pce[i,k]²

Rearranging to second-order cone (SOC) form:
    [(limit − w_pce[i,0]) / λ ;  √psi_norms[k+1]·w_pce[i,k]  ∀k≥1] ∈ SOC
"""
function constraint_voltage_ihd_chance(pm::sHHC_QCQP, i::Int;
                                        nw::Int=fundamental(pm))
    λ         = _pce_lambda(pm)
    psi_norms = _pce_norms(pm)
    ihdmax    = _PMs.ref(pm, nw,              :bus, i, "ihdmax")
    v_fund    = _PMs.ref(pm, fundamental(pm), :bus, i, "vm")

    limit  = (ihdmax * v_fund)^2
    K_high = 1:length(psi_norms)-1     # k = 1, …, δ

    w_pce = _PMs.var(pm, nw, :w_pce)

    JuMP.@constraint(pm.model,
        [(limit - w_pce[i, 0]) / λ;
         [sqrt(psi_norms[k+1]) * w_pce[i, k] for k in K_high]]
        in JuMP.SecondOrderCone()
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# 8. THD voltage chance constraint
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_voltage_thd_chance(pm::sHHC_QCQP, i; harmonics, K)

One-sided Cantelli chance constraint on the total harmonic distortion (THD)
of the bus voltage at bus i, summed over all harmonic network IDs in `harmonics`.

    Σ_h w_h[i,0] + λ · √(Σ_{k≥1} psi_norms[k+1] · (Σ_h w_h[i,k])²)
        ≤ (thdmax · v_fund)²

Defining W0 = Σ_h w_h[i,0] and Wk = Σ_h w_h[i,k], the SOC form is:
    [(limit − W0) / λ ;  √psi_norms[k+1]·Wk  ∀k≥1] ∈ SOC
"""
function constraint_voltage_thd_chance(pm::sHHC_QCQP, i::Int;
                                        harmonics, K=0:2)
    λ         = _pce_lambda(pm)
    psi_norms = _pce_norms(pm)
    thdmax    = _PMs.ref(pm, fundamental(pm), :bus, i, "thdmax")
    v_fund    = _PMs.ref(pm, fundamental(pm), :bus, i, "vm")

    limit  = (thdmax * v_fund)^2
    K_high = [k for k in K if k != 0]

    W0 = sum(_PMs.var(pm, nw, :w_pce)[i, 0] for nw in harmonics)
    Wk = [sum(_PMs.var(pm, nw, :w_pce)[i, k] for nw in harmonics) for k in K_high]

    JuMP.@constraint(pm.model,
        [(limit - W0) / λ;
         [sqrt(psi_norms[k+1]) * Wk[idx] for (idx, k) in enumerate(K_high)]]
        in JuMP.SecondOrderCone()
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# 9. RMS voltage chance constraint
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_voltage_rms_chance(pm::sHHC_QCQP, i; nw_ids, K)

One-sided Cantelli chance constraint on the RMS bus voltage at bus i,
restricted to the harmonic networks in `nw_ids` (fundamental excluded).

The fundamental voltage V_fund is treated as deterministic, so the harmonic
contribution must satisfy:

    Σ_h w_h[i,0] + λ · √(Σ_{k≥1} psi_norms[k+1] · (Σ_h w_h[i,k])²)
        ≤ vmaxrms² − v_fund²

The right-hand side margin reserves room for the deterministic fundamental
component. SOC form with W0 = Σ_h w_h[i,0] and Wk = Σ_h w_h[i,k]:
    [(limit − W0) / λ ;  √psi_norms[k+1]·Wk  ∀k≥1] ∈ SOC
"""
function constraint_voltage_rms_chance(pm::sHHC_QCQP, i::Int;
                                        nw_ids, K=0:2)
    λ         = _pce_lambda(pm)
    psi_norms = _pce_norms(pm)
    vmaxrms   = _PMs.ref(pm, fundamental(pm), :bus, i, "vmaxrms")
    v_fund    = _PMs.ref(pm, fundamental(pm), :bus, i, "vm")

    limit  = vmaxrms^2 - v_fund^2
    K_high = [k for k in K if k != 0]

    W0 = sum(_PMs.var(pm, nw, :w_pce)[i, 0] for nw in nw_ids)
    Wk = [sum(_PMs.var(pm, nw, :w_pce)[i, k] for nw in nw_ids) for k in K_high]

    JuMP.@constraint(pm.model,
        [(limit - W0) / λ;
         [sqrt(psi_norms[k+1]) * Wk[idx] for (idx, k) in enumerate(K_high)]]
        in JuMP.SecondOrderCone()
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# 10. RMS current chance constraint
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_current_rms_chance(pm::sHHC_QCQP, b; nw_ids, K)

One-sided Cantelli chance constraint on the RMS current through branch b,
restricted to the harmonic networks in `nw_ids` (fundamental excluded).

The harmonic current contribution must satisfy:

    Σ_h j_h[f_idx,0] + λ · √(Σ_{k≥1} psi_norms[k+1] · (Σ_h j_h[f_idx,k])²)
        ≤ i_rated² − i_fund²

where f_idx = (b, f_bus, t_bus) is the from-arc index, i_rated is the branch
thermal rating, and i_fund is the (deterministic) fundamental current magnitude.
The margin i_rated² − i_fund² is the remaining headroom for harmonic currents.

SOC form with J0 = Σ_h j_h[f_idx,0] and Jk = Σ_h j_h[f_idx,k]:
    [(limit − J0) / λ ;  √psi_norms[k+1]·Jk  ∀k≥1] ∈ SOC
"""
function constraint_current_rms_chance(pm::sHHC_QCQP, b::Int;
                                        nw_ids, K=0:2)
    λ         = _pce_lambda(pm)
    psi_norms = _pce_norms(pm)
    branch    = _PMs.ref(pm, fundamental(pm), :branch, b)
    f_bus     = branch["f_bus"]
    t_bus     = branch["t_bus"]
    f_idx     = (b, f_bus, t_bus)
    i_rated   = branch["c_rating"]
    i_fund    = branch["cm_fr"]

    limit  = i_rated^2 - i_fund^2
    K_high = [k for k in K if k != 0]

    J0 = sum(_PMs.var(pm, nw, :j_pce)[f_idx, 0] for nw in nw_ids)
    Jk = [sum(_PMs.var(pm, nw, :j_pce)[f_idx, k] for nw in nw_ids) for k in K_high]

    JuMP.@constraint(pm.model,
        [(limit - J0) / λ;
         [sqrt(psi_norms[k+1]) * Jk[idx] for (idx, k) in enumerate(K_high)]]
        in JuMP.SecondOrderCone()
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# Transformer PCE constraints
# ══════════════════════════════════════════════════════════════════════════════

"""
    constraint_xfmr_core_magnetization_pce(pm::sHHC_QCQP, x; nw, K)

PCE expansion of the transformer core magnetisation (one constraint pair per k ∈ K).

In the SHHC formulation network parameters are deterministic, so the magnetising
current is also deterministic. For the linearised (no-load loss neglected) case,
every PCE coefficient is fixed to zero:

    cmrx_pce[x, k] = 0    ∀ k ∈ K
    cmix_pce[x, k] = 0    ∀ k ∈ K

The nonlinear saturation case (with interpolation functions `Im_A`, `Im_B`) would
require a PCE expansion of the nonlinear function of the excitation voltages —
this is deferred to future work; a warning is emitted if saturation data is found.
"""
function constraint_xfmr_core_magnetization_pce(pm::sHHC_QCQP, x::Int;
                                                  nw::Int=fundamental(pm), K=0:2)
    xfmr = _PMs.ref(pm, nw, :xfmr, x)
    if haskey(xfmr, "Hᴵ") && nw in xfmr["Hᴵ"]
        @warn "constraint_xfmr_core_magnetization_pce: nonlinear saturation data " *
              "found for xfmr $x at nw=$nw but is not yet supported in the PCE " *
              "formulation; magnetising current is set to zero."
    end

    cmrx_pce = _PMs.var(pm, nw, :cmrx_pce)
    cmix_pce = _PMs.var(pm, nw, :cmix_pce)

    for k in K
        JuMP.@constraint(pm.model, cmrx_pce[x, k] == 0.0)
        JuMP.@constraint(pm.model, cmix_pce[x, k] == 0.0)
    end
end

"""
    constraint_xfmr_core_voltage_drop_pce(pm::sHHC_QCQP, x; nw, K)

PCE voltage drop across the transformer leakage reactance (one constraint pair per k ∈ K).

The series reactance xsc is a deterministic parameter, so the constraint is linear:

    vrx_pce[f_idx, k] = erx_pce[x, k] − xsc · csix_pce[f_idx, k]
    vix_pce[f_idx, k] = eix_pce[x, k] + xsc · csrx_pce[f_idx, k]
"""
function constraint_xfmr_core_voltage_drop_pce(pm::sHHC_QCQP, x::Int;
                                                 nw::Int=fundamental(pm), K=0:2)
    xfmr  = _PMs.ref(pm, nw, :xfmr, x)
    f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
    f_idx = (x, f_bus, t_bus)
    xsc   = xfmr["xsc"]

    erx_pce  = _PMs.var(pm, nw, :erx_pce)
    eix_pce  = _PMs.var(pm, nw, :eix_pce)
    vrx_pce  = _PMs.var(pm, nw, :vrx_pce)
    vix_pce  = _PMs.var(pm, nw, :vix_pce)
    csrx_pce = _PMs.var(pm, nw, :csrx_pce)
    csix_pce = _PMs.var(pm, nw, :csix_pce)

    for k in K
        JuMP.@constraint(pm.model,
            vrx_pce[f_idx, k] == erx_pce[x, k] - xsc * csix_pce[f_idx, k])
        JuMP.@constraint(pm.model,
            vix_pce[f_idx, k] == eix_pce[x, k] + xsc * csrx_pce[f_idx, k])
    end
end

"""
    constraint_xfmr_core_voltage_phase_shift_pce(pm::sHHC_QCQP, x; nw, K)

PCE expansion of the transformer ideal-ratio voltage phase shift (one constraint
pair per k ∈ K).

The turn-ratio phasors `tr` and `ti` are deterministic, giving linear constraints:

    erx_pce[x, k] = tr · vrx_pce[t_idx, k] − ti · vix_pce[t_idx, k]
    eix_pce[x, k] = tr · vix_pce[t_idx, k] + ti · vrx_pce[t_idx, k]
"""
function constraint_xfmr_core_voltage_phase_shift_pce(pm::sHHC_QCQP, x::Int;
                                                        nw::Int=fundamental(pm), K=0:2)
    xfmr  = _PMs.ref(pm, nw, :xfmr, x)
    f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
    t_idx = (x, t_bus, f_bus)
    tr    = xfmr["tr"]; ti = xfmr["ti"]

    erx_pce = _PMs.var(pm, nw, :erx_pce)
    eix_pce = _PMs.var(pm, nw, :eix_pce)
    vrx_pce = _PMs.var(pm, nw, :vrx_pce)
    vix_pce = _PMs.var(pm, nw, :vix_pce)

    for k in K
        JuMP.@constraint(pm.model,
            erx_pce[x, k] == tr * vrx_pce[t_idx, k] - ti * vix_pce[t_idx, k])
        JuMP.@constraint(pm.model,
            eix_pce[x, k] == tr * vix_pce[t_idx, k] + ti * vrx_pce[t_idx, k])
    end
end

"""
    constraint_xfmr_core_current_balance_pce(pm::sHHC_QCQP, x; nw, K)

PCE expansion of the transformer core current balance (one constraint pair per k ∈ K).

Derived from conj(t)·(iˢ_fr − iᵐ − gsh·e) + iˢ_to = 0 expanded into real/imaginary
parts; all parameters (tr, ti, gsh) are deterministic so the constraints are linear:

    tr·(csrx_pce[f,k] − cmrx_pce[x,k] − gsh·erx_pce[x,k])
    + ti·(csix_pce[f,k] − cmix_pce[x,k] − gsh·eix_pce[x,k]) + csrx_pce[t,k] = 0

    tr·(csix_pce[f,k] − cmix_pce[x,k] − gsh·eix_pce[x,k])
    − ti·(csrx_pce[f,k] − cmrx_pce[x,k] − gsh·erx_pce[x,k]) + csix_pce[t,k] = 0
"""
function constraint_xfmr_core_current_balance_pce(pm::sHHC_QCQP, x::Int;
                                                    nw::Int=fundamental(pm), K=0:2)
    xfmr  = _PMs.ref(pm, nw, :xfmr, x)
    f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
    f_idx = (x, f_bus, t_bus); t_idx = (x, t_bus, f_bus)
    tr    = xfmr["tr"]; ti = xfmr["ti"]; gsh = xfmr["gsh"]

    erx_pce  = _PMs.var(pm, nw, :erx_pce)
    eix_pce  = _PMs.var(pm, nw, :eix_pce)
    csrx_pce = _PMs.var(pm, nw, :csrx_pce)
    csix_pce = _PMs.var(pm, nw, :csix_pce)
    cmrx_pce = _PMs.var(pm, nw, :cmrx_pce)
    cmix_pce = _PMs.var(pm, nw, :cmix_pce)

    for k in K
        JuMP.@constraint(pm.model,
            tr * (csrx_pce[f_idx,k] - cmrx_pce[x,k] - gsh * erx_pce[x,k])
            + ti * (csix_pce[f_idx,k] - cmix_pce[x,k] - gsh * eix_pce[x,k])
            + csrx_pce[t_idx, k]
            == 0.0
        )
        JuMP.@constraint(pm.model,
            tr * (csix_pce[f_idx,k] - cmix_pce[x,k] - gsh * eix_pce[x,k])
            - ti * (csrx_pce[f_idx,k] - cmrx_pce[x,k] - gsh * erx_pce[x,k])
            + csix_pce[t_idx, k]
            == 0.0
        )
    end
end

"""
    constraint_xfmr_winding_config_pce(pm::sHHC_QCQP, x; nw, K)

PCE expansion of the transformer winding voltage-drop equations (one constraint
pair per winding per k ∈ K).

The winding resistance r (and grounding impedances re, xe for zero-sequence) are
deterministic, so the constraint is linear for each PCE index k:

- Positive/negative sequence (h ∉ 𝓗⁰):
    vrx_pce[idx,k] = vr_pce[i,k] − r · crx_pce[idx,k]
    vix_pce[idx,k] = vi_pce[i,k] − r · cix_pce[idx,k]

- Zero sequence, grounded (gnd == 1):
    vrx_pce[idx,k] = vr_pce[i,k] − (r+3re)·crx_pce[idx,k] + 3xe·cix_pce[idx,k]
    vix_pce[idx,k] = vi_pce[i,k] − (r+3re)·cix_pce[idx,k] − 3xe·crx_pce[idx,k]

- Zero sequence, ungrounded (gnd ≠ 1):
    crx_pce[idx,k] = 0,  cix_pce[idx,k] = 0  (winding disconnected)
"""
function constraint_xfmr_winding_config_pce(pm::sHHC_QCQP, x::Int;
                                              nw::Int=fundamental(pm), K=0:2)
    xfmr  = _PMs.ref(pm, nw, :xfmr, x)
    f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
    w_bus = [f_bus, t_bus]
    w_idx = [(x, f_bus, t_bus), (x, t_bus, f_bus)]

    r   = [_PMs.ref(pm, nw, :xfmr, x, nk) for nk in ["r1", "r2"]]
    re  = [_PMs.ref(pm, nw, :xfmr, x, nk) for nk in ["re1","re2"]]
    xe  = [_PMs.ref(pm, nw, :xfmr, x, nk) for nk in ["xe1","xe2"]]
    gnd = [_PMs.ref(pm, nw, :xfmr, x, nk) for nk in ["gnd1","gnd2"]]

    vr_pce  = _PMs.var(pm, nw, :vr_pce)
    vi_pce  = _PMs.var(pm, nw, :vi_pce)
    vrx_pce = _PMs.var(pm, nw, :vrx_pce)
    vix_pce = _PMs.var(pm, nw, :vix_pce)
    crx_pce = _PMs.var(pm, nw, :crx_pce)
    cix_pce = _PMs.var(pm, nw, :cix_pce)

    for w in 1:2
        i   = w_bus[w]; idx = w_idx[w]
        rw  = r[w]; rew = re[w]; xew = xe[w]; gw = gnd[w]

        for k in K
            if !is_zero_sequence(nw)
                JuMP.@constraint(pm.model,
                    vrx_pce[idx,k] == vr_pce[i,k] - rw * crx_pce[idx,k])
                JuMP.@constraint(pm.model,
                    vix_pce[idx,k] == vi_pce[i,k] - rw * cix_pce[idx,k])
            end

            if is_zero_sequence(nw) && gw == 1
                JuMP.@constraint(pm.model,
                    vrx_pce[idx,k] == vr_pce[i,k]
                                      - (rw + 3rew) * crx_pce[idx,k]
                                      + 3xew * cix_pce[idx,k])
                JuMP.@constraint(pm.model,
                    vix_pce[idx,k] == vi_pce[i,k]
                                      - (rw + 3rew) * cix_pce[idx,k]
                                      - 3xew * crx_pce[idx,k])
            end

            if is_zero_sequence(nw) && gw != 1
                JuMP.@constraint(pm.model, crx_pce[idx,k] == 0.0)
                JuMP.@constraint(pm.model, cix_pce[idx,k] == 0.0)
            end
        end
    end
end

"""
    constraint_xfmr_winding_current_balance_pce(pm::sHHC_QCQP, x; nw, K)

PCE expansion of the transformer winding current-balance equations (one constraint
pair per winding per k ∈ K).

The shunt conductance g_sh and susceptance b_sh are deterministic parameters, so
the constraint is linear for each PCE index k. The winding configuration (`cnf`)
determines which form applies:

- Positive/negative sequence, or zero-sequence Y/Z windings:
    crx_pce[idx,k] = csrx_pce[idx,k] − g_sh·vrx_pce[idx,k] + b_sh·vix_pce[idx,k]
    cix_pce[idx,k] = csix_pce[idx,k] − g_sh·vix_pce[idx,k] − b_sh·vrx_pce[idx,k]

- Zero-sequence delta (D) winding with r ≠ 0: additional delta loop term (−vrx/r).
"""
function constraint_xfmr_winding_current_balance_pce(pm::sHHC_QCQP, x::Int;
                                                       nw::Int=fundamental(pm), K=0:2)
    xfmr  = _PMs.ref(pm, nw, :xfmr, x)
    f_bus = xfmr["f_bus"]; t_bus = xfmr["t_bus"]
    w_idx = [(x, f_bus, t_bus), (x, t_bus, f_bus)]

    r    = [_PMs.ref(pm, nw, :xfmr, x, nk) for nk in ["r1","r2"]]
    cnf  = [_PMs.ref(pm, nw, :xfmr, x, nk) for nk in ["cnf1","cnf2"]]
    b_sh = [0.0, 0.0]
    g_sh = [0.0, 0.0]

    vrx_pce  = _PMs.var(pm, nw, :vrx_pce)
    vix_pce  = _PMs.var(pm, nw, :vix_pce)
    crx_pce  = _PMs.var(pm, nw, :crx_pce)
    cix_pce  = _PMs.var(pm, nw, :cix_pce)
    csrx_pce = _PMs.var(pm, nw, :csrx_pce)
    csix_pce = _PMs.var(pm, nw, :csix_pce)

    for w in 1:2
        idx = w_idx[w]; rw = r[w]; bw = b_sh[w]; gw = g_sh[w]; cw = cnf[w]

        for k in K
            if !is_zero_sequence(nw)
                JuMP.@constraint(pm.model,
                    crx_pce[idx,k] == csrx_pce[idx,k] - gw*vrx_pce[idx,k] + bw*vix_pce[idx,k])
                JuMP.@constraint(pm.model,
                    cix_pce[idx,k] == csix_pce[idx,k] - gw*vix_pce[idx,k] - bw*vrx_pce[idx,k])
            end

            if is_zero_sequence(nw) && cw in ['Y','Z']
                JuMP.@constraint(pm.model,
                    crx_pce[idx,k] == csrx_pce[idx,k] - gw*vrx_pce[idx,k] + bw*vix_pce[idx,k])
                JuMP.@constraint(pm.model,
                    cix_pce[idx,k] == csix_pce[idx,k] - gw*vix_pce[idx,k] - bw*vrx_pce[idx,k])
            end

            if is_zero_sequence(nw) && cw in ['D'] && rw != 0.0
                JuMP.@constraint(pm.model,
                    crx_pce[idx,k] == csrx_pce[idx,k] - gw*vrx_pce[idx,k] + bw*vix_pce[idx,k]
                                      - vrx_pce[idx,k] / rw)
                JuMP.@constraint(pm.model,
                    cix_pce[idx,k] == csix_pce[idx,k] - gw*vix_pce[idx,k] - bw*vrx_pce[idx,k]
                                      - vix_pce[idx,k] / rw)
            end
        end
    end
end
