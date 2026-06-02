################################################################################
# HarmonicPowerModels.jl                                                       #
# Stochastic Harmonic Hosting Capacity (SHHC) problem builders.                #
#                                                                               #
# Part A: sHHC_QCQP — full PCE-expanded QCQP via PowerModels infrastructure.  #
# Part B: sHHC_LP   — absolute-equality LP solved directly with JuMP.          #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

import LinearAlgebra: I as _eye

# ──────────────────────────────────────────────────────────────────────────────
# Helper
# ──────────────────────────────────────────────────────────────────────────────

"""
    harmonic_nw_ids(pm)

Return the sorted list of non-fundamental network ids for `pm`.
These correspond to the harmonic orders h ≠ 1 in the multinetwork model.
"""
harmonic_nw_ids(pm) = sort([n for n in _PMs.nw_ids(pm) if n ≠ fundamental(pm)])

# ══════════════════════════════════════════════════════════════════════════════
# Part A — sHHC_QCQP
# ══════════════════════════════════════════════════════════════════════════════

# ── Solver entry point ─────────────────────────────────────────────────────────

"""
    solve_shhc(hdata, sHHC_QCQP, optimizer; pce_params, kwargs...)

Solve the stochastic harmonic hosting capacity problem as a non-convex QCQP.

The function:
1. Precomputes all PCE quantities via `update_hdata_with_pce_data!`.
2. Initialises fairness-principle auxiliary data.
3. Delegates to PowerModels' multinetwork `solve_model` loop with `build_shhc`.

`pce_params` must contain the keys expected by `update_hdata_with_pce_data!`
(see that function's docstring).
"""
function solve_shhc(hdata::Dict, ::Type{sHHC_QCQP}, optimizer;
                    pce_params::Dict, kwargs...)
    update_hdata_with_pce_data!(hdata, pce_params)
    # Pre-compute fundamental HPF results so branch["cm_fr"] and bus["vm"] are
    # available for constraint_current_rms_chance and constraint_voltage_*_chance.
    update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, optimizer)
    update_hdata_with_fairness_principle_data!(hdata, sHHC_QCQP, optimizer)

    return _PMs.solve_model(hdata, sHHC_QCQP, optimizer, build_shhc;
                            ref_extensions      = [ref_add_filter!, ref_add_xfmr!],
                            solution_processors = [_HPM.sol_data_model!],
                            multinetwork        = true,
                            kwargs...)
end

# ── Problem builder ────────────────────────────────────────────────────────────

"""
    build_shhc(pm::sHHC_QCQP)

Build the stochastic HHC QCQP model within the PowerModels solve loop.

Structure
---------
- Variables: PCE-expanded voltage, current, and lifted squared-magnitude
  variables for every harmonic network (nw ≠ 1).
- Objective: maximise total hosting capacity (as in the deterministic model).
- Linear constraints: KCL, Ohm's law, reference bus, load current PCE.
- Quadratic constraints: lifting equations linking PCE variables to w_pce/j_pce.
- SOC chance constraints: Cantelli-based IHD, THD, RMS voltage and current.
"""
function build_shhc(pm::sHHC_QCQP)
    # Bridge SOC constraints to non-convex quadratic form so that Ipopt (and
    # other NLP solvers) can handle the Cantelli chance constraints.
    JuMP.add_bridge(pm.model, _MOI.Bridges.Constraint.SOCtoNonConvexQuadBridge)

    K      = pm.data["pce"]["K"]
    nw_h   = harmonic_nw_ids(pm)

    # ── Variables ─────────────────────────────────────────────────────────────
    for n in nw_h
        variable_fairness_principle(pm, nw=n, bounded=true)
        variable_load_current_magnitude(pm, nw=n, bounded=true)

        variable_bus_voltage_pce(pm, nw=n, K=K)
        variable_branch_current_pce(pm, nw=n, K=K)
        variable_load_current_pce(pm, nw=n, K=K)
        variable_lifted_voltage_pce(pm, nw=n, K=K)
        variable_lifted_current_pce(pm, nw=n, K=K)

        if !isempty(_PMs.ids(pm, :xfmr, nw=n))
            variable_xfmr_voltage_pce(pm, nw=n, K=K)
            variable_xfmr_current_pce(pm, nw=n, K=K)
        end
    end

    # ── Objective ─────────────────────────────────────────────────────────────
    objective_maximum_hosting_capacity(pm)

    # ── Fairness principle ────────────────────────────────────────────────────
    for n in nw_h
        constraint_fairness_principle(pm, nw=n)
    end

    # ── Reference bus: all PCE coefficients fixed to zero (harmonic source = 0)
    for n in nw_h
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_ref_bus_voltage_pce(pm, i, nw=n, K=K)
        end
    end

    # ── Network equality constraints (linear in PCE variables) ────────────────
    for n in nw_h
        ref_buses_n = Set(_PMs.ids(pm, :ref_buses, nw=n))

        # KCL — skip the reference bus because its voltage is fixed and the
        # generator current (which balances KCL there) is not modelled as a
        # PCE-expanded variable.  Applying KCL at the ref bus would incorrectly
        # force all branch arc currents to zero.
        for i in _PMs.ids(pm, :bus, nw=n)
            i in ref_buses_n && continue
            constraint_current_balance_pce(pm, i, nw=n, K=K)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_voltage_drop_pce(pm, b, nw=n, K=K)

            # Arc current conservation: for branches with zero shunt admittance,
            # the total from-arc current and to-arc current are equal in magnitude
            # and opposite in sign (cr_fr + cr_to = 0).  This links the two arc
            # variables that would otherwise be disconnected when the ref-bus KCL
            # is not enforced.
            branch = _PMs.ref(pm, n, :branch, b)
            f_bus  = branch["f_bus"]
            t_bus  = branch["t_bus"]
            f_idx  = (b, f_bus, t_bus)
            t_idx  = (b, t_bus, f_bus)
            g_fr = get(branch, "g_fr", 0.0)
            b_fr = get(branch, "b_fr", 0.0)
            g_to = get(branch, "g_to", 0.0)
            b_to = get(branch, "b_to", 0.0)
            if iszero(g_fr) && iszero(b_fr) && iszero(g_to) && iszero(b_to)
                cr_pce = _PMs.var(pm, n, :cr_pce)
                ci_pce = _PMs.var(pm, n, :ci_pce)
                for k in K
                    JuMP.@constraint(pm.model, cr_pce[f_idx, k] + cr_pce[t_idx, k] == 0)
                    JuMP.@constraint(pm.model, ci_pce[f_idx, k] + ci_pce[t_idx, k] == 0)
                end
            end
        end

        for d in _PMs.ids(pm, :load, nw=n)
            constraint_load_current_pce(pm, d, nw=n, K=K)
        end

        for x in _PMs.ids(pm, :xfmr, nw=n)
            constraint_xfmr_core_magnetization_pce(pm, x, nw=n, K=K)
            constraint_xfmr_core_voltage_drop_pce(pm, x, nw=n, K=K)
            constraint_xfmr_core_voltage_phase_shift_pce(pm, x, nw=n, K=K)
            constraint_xfmr_core_current_balance_pce(pm, x, nw=n, K=K)
            constraint_xfmr_winding_config_pce(pm, x, nw=n, K=K)
            constraint_xfmr_winding_current_balance_pce(pm, x, nw=n, K=K)
        end
    end

    # ── Lifting constraints (quadratic coupling PCE vars → w_pce / j_pce) ─────
    for n in nw_h
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_lifted_voltage_pce(pm, i, nw=n, K=K)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_lifted_current_pce(pm, b, nw=n, K=K)
        end
    end

    # ── Cantelli chance constraints (global — couple all harmonic networks) ────
    for i in _PMs.ids(pm, :bus, nw=fundamental(pm))
        # IHD: one constraint per bus per harmonic
        for n in nw_h
            constraint_voltage_ihd_chance(pm, i, nw=n)
        end
        # THD and RMS: one constraint per bus across all harmonics
        constraint_voltage_thd_chance(pm, i, harmonics=nw_h, K=K)
        constraint_voltage_rms_chance(pm, i, nw_ids=nw_h, K=K)
    end

    for b in _PMs.ids(pm, :branch, nw=fundamental(pm))
        constraint_current_rms_chance(pm, b, nw_ids=nw_h, K=K)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Part B — sHHC_LP
# ══════════════════════════════════════════════════════════════════════════════

# ── Internal: build numeric harmonic admittance matrix ────────────────────────

"""
    _harmonic_admittance(nw_data)

Build a numeric harmonic admittance matrix from `nw_data`, which already
contains harmonically-scaled branch impedances (br_r *= √h, br_x *= h as
applied by `replicate`).

The matrix is indexed by bus integer id (1-based, assumed contiguous).
"""
function _harmonic_admittance(nw_data::Dict, bus_idx::Dict{Int,Int})
    Nn = length(nw_data["bus"])
    Yh = zeros(ComplexF64, Nn, Nn)

    for (_, branch) in nw_data["branch"]
        f = bus_idx[branch["f_bus"]]
        t = bus_idx[branch["t_bus"]]
        r = branch["br_r"]
        x = branch["br_x"]
        if abs(r) + abs(x) > 0.0
            y_s  = 1 / (r + im * x)
            Yh[f, f] += y_s
            Yh[t, t] += y_s
            Yh[f, t] -= y_s
            Yh[t, f] -= y_s
        end
        # Branch shunt (already scaled by h in replicate)
        b_fr = get(branch, "b_fr", 0.0)
        b_to = get(branch, "b_to", 0.0)
        Yh[f, f] += im * b_fr
        Yh[t, t] += im * b_to
    end

    # Transformer stamp: for non-zero sequence, each transformer contributes a
    # 2-port admittance with complex turns ratio a = tr + im*ti (unit magnitude).
    # I_f = y_s*(V_f - a*V_t),  I_t = -conj(a)*I_f  (from core current balance)
    if haskey(nw_data, "xfmr")
        for (_, xfmr) in nw_data["xfmr"]
            f   = bus_idx[xfmr["f_bus"]]
            t   = bus_idx[xfmr["t_bus"]]
            r_s = get(xfmr, "r1", 0.0) + get(xfmr, "r2", 0.0)
            x_s = get(xfmr, "xsc", 0.0)
            abs(r_s) + abs(x_s) > 0.0 || continue
            y_s = 1 / (r_s + im * x_s)
            a   = get(xfmr, "tr", 1.0) + im * get(xfmr, "ti", 0.0)
            Yh[f, f] += y_s
            Yh[t, t] += y_s
            Yh[f, t] -= a * y_s
            Yh[t, f] -= conj(a) * y_s
        end
    end

    return Yh
end

# ── Solver entry point ─────────────────────────────────────────────────────────

"""
    solve_shhc(hdata, sHHC_LP, optimizer; pce_params, kwargs...)

Solve the stochastic HHC problem as a pure LP under absolute-equality fairness.

Under this principle all harmonic units at harmonic h receive the same scalar
hosting capacity I_h. Because the harmonic admittance matrix is linear,
every bus voltage and branch current PCE coefficient is linear in I_h:
    V_{n,h,k} = α_{n,h,k} · I_h
Consequently the Cantelli chance constraints reduce to linear upper bounds
on I_h (see `_shhc_lp_bounds` below), and the problem is a simple LP.

Steps
-----
1. Precompute PCE data (`update_hdata_with_pce_data!`).
2. Solve the fundamental harmonic power flow.
3. Build harmonic admittance matrices Z_h from scaled network data.
4. Compute voltage-response scalars α_{n,h,k} via Z_h.
5. Compute squared-magnitude PCE scalars C_{n,h,k} via the M tensor.
6. Derive per-harmonic IHD Cantelli bounds on I_h.
7. Build and solve the LP with JuMP.

Returns a result Dict with `"termination_status"`, `"objective_value"`, and
`"solution"` (Dict mapping harmonic nw id → optimal I_h).
"""
function solve_shhc(hdata::Dict, ::Type{sHHC_LP}, optimizer;
                    pce_params::Dict, hpf_optimizer=optimizer, kwargs...)
    update_hdata_with_pce_data!(hdata, pce_params)
    # The HPF requires a nonlinear solver (e.g. Ipopt); hpf_optimizer lets
    # callers pass a different solver than the LP optimizer.
    update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, hpf_optimizer)

    pce      = hdata["pce"]
    M        = pce["M"]
    lambda   = pce["lambda"]
    psi_norms = pce["psi_norms"]
    K        = pce["K"]
    cos_phi  = pce["cos_phi"]
    sin_phi  = pce["sin_phi"]

    fund_nw  = hdata["nw"]["1"]
    Nn       = length(fund_nw["bus"])

    # Map original bus IDs (may be 0-based) to 1-based matrix indices.
    bus_ids_sorted = sort([bus["bus_i"] for (_, bus) in fund_nw["bus"]])
    bus_idx        = Dict(id => i for (i, id) in enumerate(bus_ids_sorted))

    harm_nw_strs = sort([nw for nw in keys(hdata["nw"]) if nw != "1"])
    harm_nw_ints = parse.(Int, harm_nw_strs)

    # ── Step 3: Reduced Z_h via Kron elimination of reference buses ───────────
    # The full nodal admittance matrix Y_h is singular (no shunt grounding in
    # harmonic networks).  Kron reduction removes reference-bus rows/columns
    # (which have V_ref = 0), yielding an invertible reduced matrix Y_red.
    # Z_red = Y_red⁻¹ gives the transfer impedance between non-ref buses.
    ref_bus_ids = Set(
        bus_idx[bus["bus_i"]]
        for (_, bus) in fund_nw["bus"]
        if bus["bus_type"] == 3
    )
    non_ref_ids = [i for i in 1:Nn if i ∉ ref_bus_ids]   # row/col indices

    Z_harm = Dict{Int, Matrix{ComplexF64}}()
    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        nw_data = hdata["nw"][nw_str]
        Yh      = _harmonic_admittance(nw_data, bus_idx)
        Y_red   = Yh[non_ref_ids, non_ref_ids]
        nr      = length(non_ref_ids)
        try
            Z_harm[nw_int] = Y_red \ Matrix{ComplexF64}(_eye, nr, nr)
        catch
            @warn "solve_shhc (LP): reduced admittance matrix singular for " *
                  "harmonic $nw_int; skipping."
        end
    end

    # ── Step 4: α_{n,h,k} — voltage response per unit hosting current ─────────
    # For non-ref bus n and load at non-ref bus d:
    #   Vr_{n,h,k} = Σ_d [ Re(Z_red[n_idx,d_idx])·cos_k − Im(·)·sin_k ] · I_h
    #   Vi_{n,h,k} = Σ_d [ Im(Z_red[n_idx,d_idx])·cos_k + Re(·)·sin_k ] · I_h
    # Reference buses have V=0 → α = 0.
    alpha_re = Dict{Int, Matrix{Float64}}()   # [nw_int][bus_id, k_idx]
    alpha_im = Dict{Int, Matrix{Float64}}()

    non_ref_pos = Dict(id => idx for (idx, id) in enumerate(non_ref_ids))

    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(Z_harm, nw_int) || continue
        haskey(cos_phi, nw_int) || continue

        Zr      = Z_harm[nw_int]
        nw_data = hdata["nw"][nw_str]
        first_d = parse(Int, first(keys(nw_data["load"])))
        cp      = cos_phi[nw_int][first_d]
        sp      = sin_phi[nw_int][first_d]

        load_buses = [bus_idx[nw_data["load"][d]["load_bus"]]
                      for d in keys(nw_data["load"])]

        ar = zeros(Nn, length(K))
        ai = zeros(Nn, length(K))

        for n in 1:Nn
            n in ref_bus_ids && continue           # V_ref = 0 → α = 0
            n_pos = non_ref_pos[n]
            for (k_idx, _) in enumerate(K)
                for bus_d in load_buses
                    bus_d in ref_bus_ids && continue
                    d_pos = non_ref_pos[bus_d]
                    z_nd  = Zr[n_pos, d_pos]
                    rz, iz = real(z_nd), imag(z_nd)
                    ar[n, k_idx] += rz * cp[k_idx] - iz * sp[k_idx]
                    ai[n, k_idx] += iz * cp[k_idx] + rz * sp[k_idx]
                end
            end
        end

        alpha_re[nw_int] = ar
        alpha_im[nw_int] = ai
    end

    # ── Step 5: C_{n,h,k} — squared-magnitude PCE scalars ────────────────────
    # C_{n,h,k} = Σ_{j,l ∈ K} M[j+1,l+1,k+1] · (α^r_{n,h,j}·α^r_{n,h,l} + α^i_{n,h,j}·α^i_{n,h,l})
    # Note: w_{n,h,k} = C_{n,h,k} · I_h²   (squared because V ∝ I_h)
    C = Dict{Int, Matrix{Float64}}()   # [nw_int][bus, k+1]

    for nw_int in harm_nw_ints
        haskey(alpha_re, nw_int) || continue
        ar = alpha_re[nw_int]
        ai = alpha_im[nw_int]

        C_h = zeros(Nn, length(K))
        for n in 1:Nn, (k_idx, _) in enumerate(K)
            for (j_idx, _) in enumerate(K), (l_idx, _) in enumerate(K)
                C_h[n, k_idx] += M[j_idx, l_idx, k_idx] *
                                  (ar[n, j_idx] * ar[n, l_idx] +
                                   ai[n, j_idx] * ai[n, l_idx])
            end
        end
        C[nw_int] = C_h
    end

    # ── Step 6: Cantelli IHD bounds → linear upper bounds on I_h ─────────────
    # IHD Cantelli: (C_{n,h,0} + λ·√(Σ_{k≥1} ψ_k·C_{n,h,k}²)) · I_h² ≤ limit²
    # Rearranged: I_h ≤ limit / √D_{n,h}
    # where D_{n,h} = C_{n,h,0} + λ·√(Σ_{k≥1} ψ_k·C_{n,h,k}²)
    #
    # THD and cross-harmonic RMS constraints are omitted here because they mix
    # different I_h values quadratically and cannot be linearised without fixing
    # the harmonic mix. Apply them as post-solve verification constraints.

    I_max = Dict(nw_int => Inf for nw_int in harm_nw_ints)

    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(C, nw_int) || continue
        C_h     = C[nw_int]
        nw_data = hdata["nw"][nw_str]

        for (bi_str, bus) in nw_data["bus"]
            haskey(bus, "ihdmax") || continue
            bus_i   = bus["bus_i"]
            n       = bus_idx[bus_i]                    # matrix row (1-based)
            fund_bus = fund_nw["bus"][bi_str]
            v_fund  = fund_bus["vm"]
            ihdmax  = bus["ihdmax"]

            C_n0    = C_h[n, 1]                         # k=0 coefficient (1-indexed)
            sigma   = sqrt(sum(psi_norms[k+1] * C_h[n, idx]^2
                               for (idx, k) in enumerate(K) if k != 0;
                               init=0.0))
            D_nh    = C_n0 + lambda * sigma

            if D_nh > 0.0
                limit_sq    = (ihdmax * v_fund)^2
                bound_ihd   = sqrt(limit_sq / D_nh)
                I_max[nw_int] = min(I_max[nw_int], bound_ihd)
            end
        end
    end

    # Cache precomputed scalars for external result reporting
    pce["alpha_re"] = alpha_re
    pce["alpha_im"] = alpha_im
    pce["C"]        = C
    pce["I_max"]    = I_max

    # ── Step 7: Build and solve the LP ────────────────────────────────────────
    # Decision variables: I_h ≥ 0 for each harmonic h.
    # Objective: max Σ_h U_h · I_h  (U_h = number of harmonic units at h).
    # Constraints: I_h ≤ I_max[h]   (derived from IHD Cantelli bounds above).
    model = JuMP.Model(optimizer)

    U_h = Dict(nw_int => length(hdata["nw"]["$nw_int"]["load"])
               for nw_int in harm_nw_ints)

    JuMP.@variable(model, I_h[h in harm_nw_ints] >= 0)

    for nw_int in harm_nw_ints
        if isfinite(I_max[nw_int])
            JuMP.@constraint(model, I_h[nw_int] <= I_max[nw_int])
        end
    end

    JuMP.@objective(model, Max,
        sum(U_h[h] * I_h[h] for h in harm_nw_ints; init=0.0))

    JuMP.optimize!(model)

    return Dict(
        "termination_status" => JuMP.termination_status(model),
        "objective_value"    => JuMP.has_values(model) ?
                                    JuMP.objective_value(model) : NaN,
        "solution"           => JuMP.has_values(model) ?
                                    Dict(h => JuMP.value(I_h[h])
                                         for h in harm_nw_ints) :
                                    Dict{Int, Float64}(),
    )
end
