################################################################################
# HarmonicPowerModels.jl                                                       #
# Stochastic Harmonic Hosting Capacity (SHHC) problem builders.                #
#                                                                               #
# Part A: sHHC_QCQP — full PCE-expanded QCQP via PowerModels infrastructure.  #
# Part B: sHHC_SOC   — absolute-equality SOCP solved directly with JuMP.        #
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

    K      = pm.data["pce"]["K"]        # full basis (lifted vars, chance constraints)
    Kp     = pm.data["pce"]["K_phys"]   # reduced single-germ set (physical vars)
    nw_h   = harmonic_nw_ids(pm)
    hunits = Set(pm.data["pce"]["harmonic_units"])   # loads that host harmonics

    # ── Variables ─────────────────────────────────────────────────────────────
    for n in nw_h
        variable_fairness_principle(pm, nw=n, bounded=true)
        variable_load_current_magnitude(pm, nw=n, bounded=true)

        # Normal loads (not harmonic units) inject no harmonic current and host
        # no harmonic capacity → pin their hosting magnitude to zero so they
        # leave the objective and the fairness constraint untouched.
        cmd = _PMs.var(pm, n, :cmd)
        for d in _PMs.ids(pm, :load, nw=n)
            d in hunits || JuMP.set_upper_bound(cmd[d], 0.0)
        end

        # physical voltage/current variables live on the reduced set K_phys;
        # lifted |·|² variables span the full basis K
        variable_bus_voltage_pce(pm, nw=n, K=Kp)
        variable_branch_current_pce(pm, nw=n, K=Kp)
        variable_load_current_pce(pm, nw=n, K=Kp)
        variable_lifted_voltage_pce(pm, nw=n, K=K)
        variable_lifted_current_pce(pm, nw=n, K=K)

        if !isempty(_PMs.ids(pm, :xfmr, nw=n))
            variable_xfmr_voltage_pce(pm, nw=n, K=Kp)
            variable_xfmr_current_pce(pm, nw=n, K=Kp)
        end
    end

    # ── Objective ─────────────────────────────────────────────────────────────
    # Maximises Σ cmd over all loads; normal-load cmd are pinned to 0 above, so
    # this is effectively the sum over harmonic units only.
    objective_maximum_hosting_capacity(pm)

    # ── Fairness principle (over harmonic units only) ─────────────────────────
    for n in nw_h
        hload_ids = sort([d for d in _PMs.ids(pm, :load, nw=n) if d in hunits])
        constraint_fairness_principle(pm, n, hload_ids)
    end

    # ── Clean buses are ideal harmonic grounds: fix all their voltage PCE
    #    coefficients to zero.  Every other bus is free, grounded through the
    #    generator subtransient admittance Y_g(h) (folded into the KCL bus-shunt
    #    terms in constraint_current_balance_pce) and the line/bus shunts.
    clean_buses = Dict(n => Set(i for i in _PMs.ids(pm, :bus, nw=n)
                               if get(_PMs.ref(pm, n, :bus, i), "is_clean", false))
                       for n in nw_h)
    for n in nw_h
        for i in clean_buses[n]
            constraint_ref_bus_voltage_pce(pm, i, nw=n, K=Kp)
        end
    end

    # ── Network equality constraints (linear in PCE variables) ────────────────
    for n in nw_h
        # KCL at every bus except the clean (ideal-ground) buses, whose voltage
        # is fixed and whose ideal-source current balances KCL there.
        for i in _PMs.ids(pm, :bus, nw=n)
            i in clean_buses[n] && continue
            constraint_current_balance_pce(pm, i, nw=n, K=Kp)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_voltage_drop_pce(pm, b, nw=n, K=Kp)

            # Series arc current conservation: the series current component
            # flowing into the from-arc equals the series current flowing out of
            # the to-arc (cr_fr + cr_to = 0).  Any pi-model shunt currents
            # (b_fr, b_to) are separate shunt-to-ground flows that enter the KCL
            # at the respective bus nodes, not the arc current balance.
            branch = _PMs.ref(pm, n, :branch, b)
            f_bus  = branch["f_bus"]
            t_bus  = branch["t_bus"]
            f_idx  = (b, f_bus, t_bus)
            t_idx  = (b, t_bus, f_bus)
            cr_pce = _PMs.var(pm, n, :cr_pce)
            ci_pce = _PMs.var(pm, n, :ci_pce)
            for k in Kp
                JuMP.@constraint(pm.model, cr_pce[f_idx, k] + cr_pce[t_idx, k] == 0)
                JuMP.@constraint(pm.model, ci_pce[f_idx, k] + ci_pce[t_idx, k] == 0)
            end
        end

        for d in _PMs.ids(pm, :load, nw=n)
            constraint_load_current_pce(pm, d, nw=n, K=Kp)
        end

        for x in _PMs.ids(pm, :xfmr, nw=n)
            constraint_xfmr_core_magnetization_pce(pm, x, nw=n, K=Kp)
            constraint_xfmr_core_voltage_drop_pce(pm, x, nw=n, K=Kp)
            constraint_xfmr_core_voltage_phase_shift_pce(pm, x, nw=n, K=Kp)
            constraint_xfmr_core_current_balance_pce(pm, x, nw=n, K=Kp)
            constraint_xfmr_winding_config_pce(pm, x, nw=n, K=Kp)
            constraint_xfmr_winding_current_balance_pce(pm, x, nw=n, K=Kp)
        end
    end

    # ── Lifting constraints (quadratic coupling PCE vars → w_pce / j_pce) ─────
    for n in nw_h
        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_lifted_voltage_pce(pm, i, nw=n, K=K, K_phys=Kp)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_lifted_current_pce(pm, b, nw=n, K=K, K_phys=Kp)
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
# Part B — sHHC_SOC
# ══════════════════════════════════════════════════════════════════════════════

# ── Internal: build numeric harmonic admittance matrix ────────────────────────

"""
    _harmonic_admittance(nw_data)

Build a numeric harmonic admittance matrix from `nw_data`, which already
contains harmonically-scaled branch impedances (br_r *= √h, br_x *= h as
applied by `replicate`).

The matrix is indexed by bus integer id (1-based, assumed contiguous).

Every term stamped here has a counterpart in the QCQP's PCE Kirchhoff current
law (`constraint_current_balance_pce`) and in PowerModels' IVR model used by the
dHHC, so that the reduced-network SOCP and the full network models describe the
*same* network: branch π-model with off-nominal complex tap, branch shunt
admittance (g_fr/b_fr, g_to/b_to), bus shunts (gs/bs), generator subtransient
admittance, and the transformer 2-port. Dropping any of them silently breaks the
dHHC ≡ sHHC_SOC equivalence on networks where it is non-zero.
"""
function _harmonic_admittance(nw_data::Dict, bus_idx::Dict{Int,Int})
    Nn = length(nw_data["bus"])
    Yh = zeros(ComplexF64, Nn, Nn)

    for (_, branch) in nw_data["branch"]
        f = bus_idx[branch["f_bus"]]
        t = bus_idx[branch["t_bus"]]
        r = branch["br_r"]
        x = branch["br_x"]
        # Off-nominal tap / phase shift (branches encoding transformers, e.g. the
        # IEEE-30 regulating transformers). PowerModels' IVR pi-model places the
        # complex tap t = tm·e^{jθ} on the from side; omitting it (t=1) introduces
        # a network-model gap vs the dHHC. The turns ratio is frequency-independent.
        tm = get(branch, "tap", 1.0); tm == 0.0 && (tm = 1.0)
        tc = tm * cis(get(branch, "shift", 0.0))
        # Branch shunt admittance y_fr = g_fr + j*b_fr (b already h-scaled in
        # replicate).  The from-side shunt sits behind the tap, exactly as in the
        # QCQP's KCL, where its current enters as (g_fr*vr − b_fr*vi)/tm².
        y_fr = get(branch, "g_fr", 0.0) + im * get(branch, "b_fr", 0.0)
        y_to = get(branch, "g_to", 0.0) + im * get(branch, "b_to", 0.0)
        if abs(r) + abs(x) > 0.0
            y_s  = 1 / (r + im * x)
            Yh[f, f] += (y_s + y_fr) / abs2(tc)
            Yh[t, t] += y_s + y_to
            Yh[f, t] -= y_s / conj(tc)
            Yh[t, f] -= y_s / tc
        else
            Yh[f, f] += y_fr / abs2(tc)
            Yh[t, t] += y_to
        end
    end

    # Bus shunt stamp: a bus shunt (capacitor bank, shunt filter, shunt reactor)
    # is a passive admittance to ground and enters the QCQP's KCL as the gs/bs
    # totals on the right-hand side, so it belongs in Y_h as a diagonal term.
    for (_, shunt) in get(nw_data, "shunt", Dict{String,Any}())
        s = bus_idx[shunt["shunt_bus"]]
        Yh[s, s] += get(shunt, "gs", 0.0) + im * get(shunt, "bs", 0.0)
    end

    # Generator stamp: at a harmonic each generator is a shunt to ground through
    # its subtransient admittance Y_g(h) = gsc + j*bsc (already h-scaled in
    # replicate).  An "infinite" generator (inf=1) is an ideal source → no shunt.
    for (_, gen) in nw_data["gen"]
        iszero(get(gen, "inf", 0)) || continue
        g = bus_idx[gen["gen_bus"]]
        Yh[g, g] += get(gen, "gsc", 0.0) + im * get(gen, "bsc", 0.0)
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
    solve_shhc(hdata, sHHC_SOC, optimizer; pce_params, hpf_optimizer, kwargs...)

Solve the stochastic HHC problem as a Second-Order Cone Program (SOCP) under
absolute-equality fairness.

Under this principle all harmonic units at harmonic h receive the same scalar
hosting capacity I_h, so the squared bus voltages and branch currents are
quadratic in I_h (computed analytically via the harmonic admittance matrix Z_h
and the PCE multiplication tensor). The chance constraints become a small
program in the per-harmonic variables I_h (with g_h = I_h²):

- IHD is per harmonic ⇒ a LINEAR bound on g_h.
- THD and RMS-voltage (per bus) and RMS-current (per branch) couple the
  harmonics quadratically ⇒ SECOND-ORDER CONE constraints in g_h.

Steps
-----
1. Precompute PCE data (`update_hdata_with_pce_data!`); the chance-constraint
   factor λ honours `pce_params["bound"]` (`:cantelli` default, or `:gaussian`).
2. Solve the fundamental harmonic power flow (`hpf_optimizer`, an NLP solver).
3. Build harmonic admittance matrices Z_h from scaled network data.
4. Compute voltage- and branch-current-response PCE coefficients via Z_h.
5. Compute the squared-magnitude PCE scalars α_V, α_I via the M tensor.
6. Build the SOCP (IHD linear bounds + THD/RMS cones) and solve with `optimizer`
   (a conic solver, e.g. Mosek or Clarabel — not HiGHS).

Returns a result Dict with `"termination_status"`, `"objective_value"`, and
`"solution"` (Dict mapping harmonic nw id → optimal I_h).

Only the absolute-equality fairness principle is representable in this
formulation: the single scalar I_h per harmonic is what makes the squared
magnitudes quadratic in one variable and the chance constraints conic. Any other
`hdata["principle"]` is rejected rather than silently ignored; use `sHHC_QCQP`,
which carries a per-unit current variable.
"""
function solve_shhc(hdata::Dict, ::Type{sHHC_SOC}, optimizer;
                    pce_params::Dict, hpf_optimizer=optimizer, kwargs...)
    principle = get(hdata, "principle", "absolute equality")
    principle == "absolute equality" || error(
        "solve_shhc (sHHC_SOC): fairness principle \"$principle\" is not " *
        "representable in the SOCP, which assigns every harmonic unit the same " *
        "scalar hosting capacity I_h. Use sHHC_QCQP for \"$principle\", or set " *
        "hdata[\"principle\"] = \"absolute equality\".")

    update_hdata_with_pce_data!(hdata, pce_params)
    # The HPF requires a nonlinear solver (e.g. Ipopt); hpf_optimizer lets
    # callers pass a different solver than the SOCP optimizer.
    update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, hpf_optimizer)

    pce      = hdata["pce"]
    pb       = pce["basis"]          # multivariate PCEBasis (one germ per harmonic unit)
    lambda   = pce["lambda"]
    lambdas  = get(pce, "lambdas", nothing)   # per-constraint λ* (empirical) or nothing
    # λ for a given constraint: per-constraint empirical value if available, else scalar
    λfor(kind, key) = (lambdas !== nothing && haskey(lambdas[kind], key)) ?
                       lambdas[kind][key] : lambda
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

    # ── Step 3: Z_h via Kron elimination of clean buses ───────────────────────
    # Clean buses (bus_type 4) are ideal harmonic grounds (V=0); they are
    # Kron-eliminated.  Every other bus is grounded through the generator
    # subtransient admittance Y_g(h) (stamped into Y_h) and the line/bus shunts,
    # which keeps the reduced matrix invertible.
    ref_bus_ids = Set(
        bus_idx[bus["bus_i"]]
        for (_, bus) in fund_nw["bus"]
        if get(bus, "is_clean", false)
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
            @warn "solve_shhc (SOCP): reduced admittance matrix singular for " *
                  "harmonic $nw_int; skipping."
        end
    end

    # ── Step 4: C_re/C_im — voltage response PCE per unit hosting current ─────
    # Independent germs: each harmonic unit d carries its OWN trig PCE vector
    # cos_phi[nw][d] / sin_phi[nw][d] (length P, nonzero only on d's germ
    # dimension). For non-ref bus n:
    #   Cr_{n,h,k} = Σ_d [ Re(Z_red[n,d])·cos_d,k − Im(Z_red[n,d])·sin_d,k ]
    #   Ci_{n,h,k} = Σ_d [ Im(Z_red[n,d])·cos_d,k + Re(Z_red[n,d])·sin_d,k ]
    # Reference buses have V=0 → response = 0; normal loads contribute zero
    # (their cos/sin vectors are all-zero).
    alpha_re = Dict{Int, Matrix{Float64}}()   # [nw_int][bus_id, k_idx]
    alpha_im = Dict{Int, Matrix{Float64}}()

    non_ref_pos = Dict(id => idx for (idx, id) in enumerate(non_ref_ids))

    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(Z_harm, nw_int) || continue
        haskey(cos_phi, nw_int) || continue

        Zr      = Z_harm[nw_int]
        nw_data = hdata["nw"][nw_str]
        cph     = cos_phi[nw_int]   # Dict: load id → length-P coefficient vector
        sph     = sin_phi[nw_int]

        ar = zeros(Nn, length(K))
        ai = zeros(Nn, length(K))

        for n in 1:Nn
            n in ref_bus_ids && continue           # V_ref = 0 → response = 0
            n_pos = non_ref_pos[n]
            for (d_str, load) in nw_data["load"]
                d     = parse(Int, d_str)
                bus_d = bus_idx[load["load_bus"]]
                bus_d in ref_bus_ids && continue
                d_pos  = non_ref_pos[bus_d]
                z_nd   = Zr[n_pos, d_pos]
                rz, iz = real(z_nd), imag(z_nd)
                cp_d   = cph[d]
                sp_d   = sph[d]
                for (k_idx, _) in enumerate(K)
                    ar[n, k_idx] += rz * cp_d[k_idx] - iz * sp_d[k_idx]
                    ai[n, k_idx] += iz * cp_d[k_idx] + rz * sp_d[k_idx]
                end
            end
        end

        alpha_re[nw_int] = ar
        alpha_im[nw_int] = ai
    end

    # ── Step 5: |V|² (per bus) and |I_b|² (per branch) PCE per harmonic ───────
    # alpha_V[h][n,:] = PCE of |V_{n,h}|²/I_h²  (sparse Galerkin square of the
    # voltage response).  alpha_I[h][b]   = PCE of |I_{b,h}|²/I_h² for the series
    # branch current I_b = y_s·(V_f − V_t), y_s = 1/(r+jx) (harmonic-scaled r,x).
    alpha_V = Dict{Int, Matrix{Float64}}()              # [nw][bus, k]
    alpha_I = Dict{Int, Dict{Int, Vector{Float64}}}()   # [nw][branch] → length-|K|

    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(alpha_re, nw_int) || continue
        ar, ai  = alpha_re[nw_int], alpha_im[nw_int]
        nw_data = hdata["nw"][nw_str]

        A = zeros(Nn, length(K))
        for n in 1:Nn
            n in ref_bus_ids && continue
            A[n, :] .= galerkin_product(pb, view(ar, n, :), view(ar, n, :)) .+
                       galerkin_product(pb, view(ai, n, :), view(ai, n, :))
        end
        alpha_V[nw_int] = A

        AI = Dict{Int, Vector{Float64}}()
        for (b_str, br) in nw_data["branch"]
            r, x = br["br_r"], br["br_x"]
            (abs(r) + abs(x) > 0.0) || continue
            f, t = bus_idx[br["f_bus"]], bus_idx[br["t_bus"]]
            ys   = 1.0 / (r + im * x); gy, by = real(ys), imag(ys)
            dVr  = view(ar, f, :) .- view(ar, t, :)
            dVi  = view(ai, f, :) .- view(ai, t, :)
            Ibr  = gy .* dVr .- by .* dVi
            Ibi  = gy .* dVi .+ by .* dVr
            AI[parse(Int, b_str)] = galerkin_product(pb, Ibr, Ibr) .+
                                    galerkin_product(pb, Ibi, Ibi)
        end
        alpha_I[nw_int] = AI
    end

    # ── Step 6: build the SOCP ────────────────────────────────────────────────
    # Variables I_h ≥ 0 (hosting capacity per harmonic) and g_h = I_h² ≥ 0.
    #   • IHD (per bus, per harmonic): Cantelli on |V_{n,h}|² ⇒ a LINEAR bound
    #     C_scalar_{n,h}·g_h ≤ (ihdmax·V₁)² (the per-harmonic case factorises).
    #   • THD and RMS-voltage (per bus): Cantelli on Σ_h|V_{n,h}|² and RMS-current
    #     (per branch): Cantelli on Σ_h|I_{b,h}|² couple the harmonics quadratically
    #     ⇒ SECOND-ORDER CONES in g_h.
    # The cross-harmonic variance sums Σ_h α_{·,h,k}·g_h over the SAME PCE basis
    # (each unit's harmonics share its germ), matching the QCQP THD/RMS chance
    # constraints.
    model = JuMP.Model(optimizer)
    # U_h = number of HARMONIC units at harmonic h (only these host capacity);
    # normal loads are present in the network but inject nothing.
    hunit_set = Set(pce["harmonic_units"])
    U_h = Dict(nw_int => count(d -> parse(Int, d) in hunit_set,
                               keys(hdata["nw"]["$nw_int"]["load"]))
               for nw_int in harm_nw_ints)

    JuMP.@variable(model, I_h[h in harm_nw_ints] >= 0)
    JuMP.@variable(model, g_h[h in harm_nw_ints] >= 0)
    for h in harm_nw_ints                      # I_h² ≤ g_h  (rotated cone)
        JuMP.@constraint(model, [g_h[h], 0.5, I_h[h]] in JuMP.RotatedSecondOrderCone())
    end

    # Per-unit harmonic current rating, mirroring the cmd ≤ load["c_rating"]
    # bound of the dHHC and QCQP models (see `variable_load_current_magnitude`).
    # Under absolute equality every harmonic unit is assigned the same I_h, so
    # the binding cap is the SMALLEST rating among the units — the network's
    # weakest unit limits the common allocation. Omitting this bound lets the
    # SOCP report capacities the units could not physically carry, and makes it
    # disagree with the dHHC/QCQP wherever the rating binds.
    for h in harm_nw_ints
        caps = [ld["c_rating"] for (d, ld) in hdata["nw"]["$h"]["load"]
                if parse(Int, d) in hunit_set && haskey(ld, "c_rating")]
        isempty(caps) && continue
        JuMP.set_upper_bound(I_h[h], minimum(caps))
        JuMP.set_upper_bound(g_h[h], minimum(caps)^2)
    end

    Kh = [k for k in K if k != 0]              # k ≥ 1 PCE indices

    # Cantelli bound  E[S] + λ√Var[S] ≤ limit  for an aggregate S whose per-harmonic
    # PCE coefficients are coef[h] (length |K|); added as a second-order cone in g_h.
    function _cantelli_soc!(coef, limit, λ)
        if limit <= 0.0
            @warn "solve_shhc (SOCP): non-positive limit $limit; constraint skipped " *
                  "(fundamental already at/over the limit)."
            return
        end
        ES = JuMP.@expression(model, sum(coef[h][1] * g_h[h] for h in keys(coef)))
        Wk = [JuMP.@expression(model, sum(coef[h][k+1] * g_h[h] for h in keys(coef))) for k in Kh]
        if λ <= 0.0                       # empirical λ*≈0 ⇒ mean bound (no variance term)
            JuMP.@constraint(model, ES <= limit); return
        end
        JuMP.@constraint(model,
            vcat((limit - ES) / λ,
                 [sqrt(psi_norms[k+1]) * Wk[j] for (j, k) in enumerate(Kh)])
            in JuMP.SecondOrderCone())
    end

    # IHD — linear per (bus, harmonic)
    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(alpha_V, nw_int) || continue
        A = alpha_V[nw_int]
        for (bi_str, bus) in hdata["nw"][nw_str]["bus"]
            haskey(bus, "ihdmax") || continue
            n      = bus_idx[bus["bus_i"]]
            v_fund = fund_nw["bus"][bi_str]["vm"]
            λ      = λfor(:ihd, (nw_int, bus["bus_i"]))
            Cscal  = A[n, 1] + λ * sqrt(sum(psi_norms[k+1] * A[n, idx]^2
                                            for (idx, k) in enumerate(K) if k != 0; init=0.0))
            Cscal > 0.0 && JuMP.@constraint(model, Cscal * g_h[nw_int] <= (bus["ihdmax"] * v_fund)^2)
        end
    end

    # THD + RMS-voltage — per bus (cross-harmonic SOC on Σ_h|V_{n,h}|²)
    for (bi_str, fbus) in fund_nw["bus"]
        n = bus_idx[fbus["bus_i"]]
        n in ref_bus_ids && continue
        coef   = Dict(h => view(alpha_V[h], n, :) for h in harm_nw_ints if haskey(alpha_V, h))
        v_fund = fbus["vm"]
        λ      = λfor(:agg, fbus["bus_i"])
        haskey(fbus, "thdmax")  && _cantelli_soc!(coef, (fbus["thdmax"] * v_fund)^2, λ)
        haskey(fbus, "vmaxrms") && _cantelli_soc!(coef, fbus["vmaxrms"]^2 - v_fund^2, λ)
    end

    # RMS-current — per branch (cross-harmonic SOC on Σ_h|I_{b,h}|²)
    for (b_str, fbr) in fund_nw["branch"]
        haskey(fbr, "c_rating") || continue
        b    = parse(Int, b_str)
        coef = Dict(h => alpha_I[h][b] for h in harm_nw_ints
                    if haskey(alpha_I, h) && haskey(alpha_I[h], b))
        isempty(coef) && continue
        _cantelli_soc!(coef, fbr["c_rating"]^2 - get(fbr, "cm_fr", 0.0)^2, λfor(:rmsi, b))
    end

    JuMP.@objective(model, Max, sum(U_h[h] * I_h[h] for h in harm_nw_ints; init=0.0))
    JuMP.optimize!(model)

    # cache precomputed scalars for external result reporting
    pce["alpha_re"] = alpha_re; pce["alpha_im"] = alpha_im
    pce["alpha_V"]  = alpha_V;  pce["alpha_I"]  = alpha_I

    return Dict(
        "termination_status" => JuMP.termination_status(model),
        "objective_value"    => JuMP.has_values(model) ? JuMP.objective_value(model) : NaN,
        "solution"           => JuMP.has_values(model) ?
                                    Dict(h => JuMP.value(I_h[h]) for h in harm_nw_ints) :
                                    Dict{Int, Float64}(),
    )
end

# ══════════════════════════════════════════════════════════════════════════════
# Part C — empirical (quantile-calibrated) iterative solve
# ══════════════════════════════════════════════════════════════════════════════

# hosting decision g (h => I_h²) extracted from a solution (SOC or QCQP form)
function _hosting_g(sol)
    s = sol["solution"]
    if haskey(s, "nw")                              # PowerModels (sHHC_QCQP) form
        g = Dict{Int,Float64}()
        for (nw_str, nwd) in s["nw"]
            nw_str == "1" && continue
            haskey(nwd, "load") || continue
            g[parse(Int, nw_str)] = maximum(ld["cmd"] for (_, ld) in nwd["load"]; init=0.0)^2
        end
        return g
    else                                            # sHHC_SOC form: Dict(h => I_h)
        return Dict(h => v^2 for (h, v) in s)
    end
end
_obj(sol) = get(sol, "objective_value", get(sol, "objective", NaN))

"""
    solve_shhc_calibrated(hdata, model, optimizer; pce_params, nsamp, maxiter, tol,
                          init_bound, hpf_optimizer, verbose) -> result

Solve the SHHC problem with EMPIRICAL, quantile-calibrated chance-constraint
factors λ*(ε) (see [`empirical_lambdas`](@ref)) for `model ∈ {sHHC_SOC, sHHC_QCQP}`.

Because the per-harmonic IHD factors are independent of the hosting decision but
the cross-harmonic THD/RMS factors are not, the routine iterates
    solve → recalibrate λ* at the solution → re-solve
until the objective converges (`|Δobj|/obj ≤ tol`) or `maxiter` is reached. The
first solve uses `init_bound` (default `"gaussian"`) to obtain a starting `g`.

The returned result is the final solution dict, augmented with
`"empirical_lambdas"`, `"iterations"` and `"objective_history"`.
"""
function solve_shhc_calibrated(hdata::Dict, ::Type{M}, optimizer;
                               pce_params::Dict, nsamp::Int=20000, maxiter::Int=6,
                               tol::Real=1e-3, init_bound::String="gaussian",
                               verbose::Bool=false, kwargs...) where {M}
    # workspace with PCE data (angle ranges, ε, α̂β̂) for λ* sampling — no HPF needed
    hd_pce = deepcopy(hdata)
    update_hdata_with_pce_data!(hd_pce, merge(pce_params, Dict("bound" => init_bound)))

    # initial solve to obtain a hosting decision g
    pp = merge(pce_params, Dict("bound" => init_bound))
    sol = solve_shhc(deepcopy(hdata), M, optimizer; pce_params=pp, kwargs...)
    g   = _hosting_g(sol)
    hist = Float64[_obj(sol)]
    λs   = nothing

    for it in 1:maxiter
        λs  = empirical_lambdas(hd_pce, g; nsamp=nsamp)
        pp2 = merge(pce_params, Dict("bound" => "empirical", "empirical_lambdas" => λs))
        sol = solve_shhc(deepcopy(hdata), M, optimizer; pce_params=pp2, kwargs...)
        ob  = _obj(sol); push!(hist, ob)
        verbose && @info "solve_shhc_calibrated" iter=it objective=ob
        if isfinite(ob) && isfinite(hist[end-1]) && ob > 0 &&
           abs(ob - hist[end-1]) / ob ≤ tol
            break
        end
        g = _hosting_g(sol)
    end

    sol["empirical_lambdas"] = λs
    sol["iterations"]        = length(hist) - 1
    sol["objective_history"] = hist
    return sol
end
