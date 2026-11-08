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
    # Pre-compute fundamental HPF results so branch["cm_fr"] and bus["vm"] are
    # available for constraint_current_rms_chance and constraint_voltage_*_chance.
    # It runs first: with a local angle reference the PCE angle supports are
    # rotated by the fundamental bus voltage angles it produces.
    update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, optimizer)
    update_hdata_with_pce_data!(hdata, pce_params)
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

# π-model stamp of one branch (series r + jx, shunts y_fr and y_to, complex tap tc on the
# from side) into Yh; shared by `_harmonic_admittance` and the robust fast evaluator
function _stamp_branch!(Yh, f, t, r, x, y_fr, y_to, tc)
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
    return Yh
end

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
    nw_data = _zero_sequence_view(nw_data)   # the network itself unless zero sequence
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
        # zero sequence: a winding pair that earths one side (`_zero_sequence_view`)
        if get(branch, "earth_to", false)
            Yh[f, f] += (1 / (r + im * x) + y_fr) / abs2(tc)
        elseif get(branch, "earth_fr", false)
            Yh[t, t] += 1 / (r + im * x) + y_to
        else
            _stamp_branch!(Yh, f, t, r, x, y_fr, y_to, tc)
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

# ar[n,k] += Re(z)·cos_k − Im(z)·sin_k and ai[n,k] += Im(z)·cos_k + Re(z)·sin_k; a function
# barrier, so the loop specialises on the coefficient vectors read from untyped dictionaries
function _accumulate_response!(ar, ai, n, nk, rz, iz, cp_d, sp_d)
    @inbounds for k_idx in 1:nk
        ar[n, k_idx] += rz * cp_d[k_idx] - iz * sp_d[k_idx]
        ai[n, k_idx] += iz * cp_d[k_idx] + rz * sp_d[k_idx]
    end
    return nothing
end

# Harmonic units that have a return path, with floating matrix indices excluded.
_harmonic_units(nw, hunits, bus_idx, floating) = Set{Int}(
    parse(Int, d) for (d, ld) in nw["load"]
    if parse(Int, d) in hunits && bus_idx[ld["load_bus"]] ∉ floating)

"""
    shhc_moments(hdata; networks = hdata["nw"]) -> NamedTuple

Steps 3–5 of [`solve_shhc`](@ref) for `sHHC_SOC`: the PCE moments of the squared bus
voltages (`alpha_V`) and series branch currents (`alpha_I`) per unit hosting current.
`hdata` must already carry the PCE data (`update_hdata_with_pce_data!`).

The harmonic admittances, loads and branch impedances are read from `networks` (nw id
string ⇒ network data), which defaults to `hdata["nw"]`. Passing the networks of another
realization of the network parameters gives that realization's moments, with the same
buses, units and PCE data — the oracle of the robust (cutting-set) formulation.
`active_units[h]` records the participating load IDs: floating zero-sequence units
have zero injection and allocation weight at that harmonic; the base fairness
weights `wts` still apply at harmonics where those units participate.
"""
function shhc_moments(hdata::Dict; networks::Dict = hdata["nw"])
    pce      = hdata["pce"]
    pb       = pce["basis"]          # multivariate PCEBasis (one germ per harmonic unit)
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
    # zero-sequence harmonics: the zero-sequence equivalent networks
    networks = Dict(nw => _zero_sequence_view(networks[nw]) for nw in harm_nw_strs)

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
    # Floating units inject nothing at this harmonic and drop out of its rows and weights.
    # Block the shared capacity only when every harmonic unit floats: no injection remains
    # to bound it. Units in earthed parts retain their common hosting capacity.
    blocked = Set{Int}()
    hunits  = Set(get(pce, "harmonic_units", Int[]))
    active_units = Dict{Int,Set{Int}}()
    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        nw_data = networks[nw_str]
        Yh      = _harmonic_admittance(nw_data, bus_idx)
        Y_red   = Yh[non_ref_ids, non_ref_ids]
        nr      = length(non_ref_ids)
        fl      = Set{Int}()
        if get(nw_data, "harmonic_sequence", :positive) === :zero
            Y_red, floating = _earth_floating(Y_red)
            fl = Set(non_ref_ids[i] for idx in floating for i in idx)
        end
        active_units[nw_int] = _harmonic_units(nw_data, hunits, bus_idx, fl)
        all_units = _harmonic_units(nw_data, hunits, bus_idx, Set{Int}())
        !isempty(all_units) && isempty(active_units[nw_int]) && push!(blocked, nw_int)
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
    # fixed allocation weights: unit d injects w_d·I_h (w_d = 1 under absolute equality)
    wts = fairness_weights(hdata, pce["harmonic_units"])
    alpha_re = Dict{Int, Matrix{Float64}}()   # [nw_int][bus_id, k_idx]
    alpha_im = Dict{Int, Matrix{Float64}}()

    non_ref_pos = Dict(id => idx for (idx, id) in enumerate(non_ref_ids))

    # ── No recombination: worst-case alignment of the unit contributions ─────
    if !get(pce, "recombination", true)
        alpha_V, alpha_I = _worst_case_moments(Z_harm, networks, harm_nw_strs, harm_nw_ints,
            bus_idx, ref_bus_ids, non_ref_pos, wts, active_units, length(K))
        return (; pb, K, psi_norms, alpha_re, alpha_im, alpha_V, alpha_I, Nn, bus_idx,
                  ref_bus_ids, harm_nw_strs, harm_nw_ints, wts, active_units, blocked)
    end

    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(Z_harm, nw_int) || continue
        haskey(cos_phi, nw_int) || continue

        Zr      = Z_harm[nw_int]
        nw_data = networks[nw_str]
        cph     = cos_phi[nw_int]   # Dict: load id → length-P coefficient vector
        sph     = sin_phi[nw_int]

        ar = zeros(Nn, length(K))
        ai = zeros(Nn, length(K))

        for n in 1:Nn
            n in ref_bus_ids && continue           # V_ref = 0 → response = 0
            n_pos = non_ref_pos[n]
            for (d_str, load) in nw_data["load"]
                d     = parse(Int, d_str)
                d in hunits && d ∉ active_units[nw_int] && continue
                bus_d = bus_idx[load["load_bus"]]
                bus_d in ref_bus_ids && continue
                d_pos  = non_ref_pos[bus_d]
                z_nd   = Zr[n_pos, d_pos]
                rz, iz = real(z_nd), imag(z_nd)
                cp_d   = cph[d]
                sp_d   = sph[d]
                w_d    = get(wts, d, 0.0)        # normal loads inject nothing anyway
                rz, iz = w_d * rz, w_d * iz
                _accumulate_response!(ar, ai, n, length(K), rz, iz, cp_d, sp_d)
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
        nw_data = networks[nw_str]

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
            # series current y_s·(V_f/t − V_t) with the off-nominal complex tap t
            # on the from side, as in the QCQP's Ohm's law and the dHHC
            tm   = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0)
            ct   = 1.0 / (tm * cis(get(br, "shift", 0.0))); cr, ci = real(ct), imag(ct)
            # zero sequence: an earthed side of a winding pair has zero voltage
            ef   = get(br, "earth_fr", false) ? 0.0 : 1.0
            et   = get(br, "earth_to", false) ? 0.0 : 1.0
            dVr  = ef .* (cr .* view(ar, f, :) .- ci .* view(ai, f, :)) .- et .* view(ar, t, :)
            dVi  = ef .* (cr .* view(ai, f, :) .+ ci .* view(ar, f, :)) .- et .* view(ai, t, :)
            Ibr  = gy .* dVr .- by .* dVi
            Ibi  = gy .* dVi .+ by .* dVr
            AI[parse(Int, b_str)] = galerkin_product(pb, Ibr, Ibr) .+
                                    galerkin_product(pb, Ibi, Ibi)
        end
        alpha_I[nw_int] = AI
    end

    return (; pb, K, psi_norms, alpha_re, alpha_im, alpha_V, alpha_I, Nn, bus_idx,
              ref_bus_ids, harm_nw_strs, harm_nw_ints, wts, active_units, blocked)
end

# Moments without recombination: every unit contribution aligns in the worst case, so a
# squared quantity per unit g_h is the square of the arithmetic sum of the magnitudes,
#   |V_{n,h}|² / g_h = (Σ_d w_d |Z_{n,d}|)²,  |I_{b,h}|² / g_h = (Σ_d w_d |y_s (Z_{f,d}/t − Z_{t,d})|)²,
# stored as the mean (k = 0) with zero variance. Normal loads inject nothing.
function _worst_case_moments(Z_harm, networks, harm_nw_strs, harm_nw_ints, bus_idx,
                             ref_bus_ids, non_ref_pos, wts, active_units, nK)
    alpha_V = Dict{Int, Matrix{Float64}}()
    alpha_I = Dict{Int, Dict{Int, Vector{Float64}}}()
    Nn = length(bus_idx)
    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(Z_harm, nw_int) || continue
        Zr, nw_data = Z_harm[nw_int], networks[nw_str]
        units = Tuple{Int,Float64}[]
        for (d_str, load) in nw_data["load"]
            d = parse(Int, d_str); d in active_units[nw_int] || continue
            bd = bus_idx[load["load_bus"]]; bd in ref_bus_ids && continue
            w  = get(wts, d, 0.0); w == 0.0 || push!(units, (non_ref_pos[bd], w))
        end
        zrow(n) = n in ref_bus_ids ? nothing : non_ref_pos[n]
        A = zeros(Nn, nK)
        for n in 1:Nn
            pn = zrow(n); pn === nothing && continue
            A[n, 1] = sum((w * abs(Zr[pn, pd]) for (pd, w) in units); init = 0.0)^2
        end
        alpha_V[nw_int] = A
        AI = Dict{Int, Vector{Float64}}()
        for (b_str, br) in nw_data["branch"]
            r, x = br["br_r"], br["br_x"]
            (abs(r) + abs(x) > 0.0) || continue
            ys = 1.0 / (r + im * x)
            tm = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0)
            ct = 1.0 / (tm * cis(get(br, "shift", 0.0)))
            pf = get(br, "earth_fr", false) ? nothing : zrow(bus_idx[br["f_bus"]])
            pt = get(br, "earth_to", false) ? nothing : zrow(bus_idx[br["t_bus"]])
            a  = zero(Float64)
            for (pd, w) in units
                zf = pf === nothing ? zero(ComplexF64) : Zr[pf, pd]
                zt = pt === nothing ? zero(ComplexF64) : Zr[pt, pd]
                a += w * abs(ys * (ct * zf - zt))
            end
            v = zeros(nK); v[1] = a^2
            AI[parse(Int, b_str)] = v
        end
        alpha_I[nw_int] = AI
    end
    return alpha_V, alpha_I
end

"""
    shhc_master(hdata, moments, optimizer) -> NamedTuple

Step 6 of [`solve_shhc`](@ref) for `sHHC_SOC`, without the chance constraints: the model,
the scaled variables I'_h = I_h/s_h and g'_h = g_h/s_h² with I_h² ≤ g_h, the unit
ratings and the normalised objective. The scaling s_h and the objective weights are
fixed by `moments` (the nominal network); add the chance constraints of one or more
network realizations with [`add_shhc_constraints!`](@ref).
"""
function shhc_master(hdata::Dict, moments, optimizer)
    (; K, alpha_V, bus_idx, harm_nw_strs, harm_nw_ints, wts) = moments
    pce     = hdata["pce"]
    fund_nw = hdata["nw"]["1"]

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
    # U_h = Σ w_d over the HARMONIC units at harmonic h, i.e., the total hosting
    # per unit of I_h (the number of participating units under absolute equality); normal loads
    # are present in the network but inject nothing.
    hunit_set = Set(pce["harmonic_units"])
    active_units = get(moments, :active_units, Dict(h => hunit_set for h in harm_nw_ints))
    U_h = Dict(nw_int => sum((wts[parse(Int, d)] for d in keys(hdata["nw"]["$nw_int"]["load"])
                              if parse(Int, d) in active_units[nw_int]); init = 0.0)
               for nw_int in harm_nw_ints)

    # Scaling: I_h ≈ 1e-3 and g_h = I_h² ≈ 1e-7 p.u. sit at the edge of the conic
    # solvers' tolerances (Clarabel returned OPTIMAL 2 % above the optimum). The
    # program is solved in I'_h = I_h/s_h with s_h the bound from the individual
    # harmonic voltage limits alone, and every constraint is divided by its limit.
    s_h = Dict{Int,Float64}()
    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        r = Inf
        if haskey(alpha_V, nw_int)
            for (bi_str, bus) in hdata["nw"][nw_str]["bus"]
                haskey(bus, "ihdmax") || continue
                a0 = alpha_V[nw_int][bus_idx[bus["bus_i"]], 1]
                L  = (bus["ihdmax"] * fund_nw["bus"][bi_str]["vm"])^2
                a0 > 0.0 && L > 0.0 && (r = min(r, sqrt(L / a0)))
            end
        end
        s_h[nw_int] = isfinite(r) ? r : 1.0
    end

    JuMP.@variable(model, Ip_h[h in harm_nw_ints] >= 0)    # I_h / s_h
    JuMP.@variable(model, gp_h[h in harm_nw_ints] >= 0)    # g_h / s_h²
    I_h = Dict(h => s_h[h] * Ip_h[h] for h in harm_nw_ints)
    g_h = Dict(h => s_h[h]^2 * gp_h[h] for h in harm_nw_ints)
    for h in harm_nw_ints                      # I_h² ≤ g_h  (rotated cone)
        JuMP.@constraint(model, [gp_h[h], 0.5, Ip_h[h]] in JuMP.RotatedSecondOrderCone())
    end

    # Per-unit harmonic current rating, mirroring the cmd ≤ load["c_rating"]
    # bound of the dHHC and QCQP models (see `variable_load_current_magnitude`).
    # Unit d is assigned w_d·I_h, so the binding cap is the SMALLEST c_rating/w_d
    # among the units (the smallest rating under absolute equality). Omitting this
    # bound lets the SOCP report capacities the units could not physically carry,
    # and makes it disagree with the dHHC/QCQP wherever the rating binds.
    for h in harm_nw_ints
        caps = [ld["c_rating"] / wts[parse(Int, d)] for (d, ld) in hdata["nw"]["$h"]["load"]
                if parse(Int, d) in active_units[h] && haskey(ld, "c_rating")]
        isempty(caps) && continue
        JuMP.set_upper_bound(Ip_h[h], minimum(caps) / s_h[h])
        JuMP.set_upper_bound(gp_h[h], (minimum(caps) / s_h[h])^2)
    end
    # No harmonic unit has a return path: fix the otherwise unbounded capacity to zero.
    for h in get(moments, :blocked, Set{Int}())
        JuMP.set_upper_bound(Ip_h[h], 0.0)
        JuMP.set_upper_bound(gp_h[h], 0.0)
    end

    # objective normalised by S = Σ_h U_h·s_h; the reported value is rescaled
    S = sum(U_h[h] * s_h[h] for h in harm_nw_ints; init = 0.0); S > 0.0 || (S = 1.0)
    JuMP.@objective(model, Max, sum(U_h[h] * s_h[h] / S * Ip_h[h] for h in harm_nw_ints; init=0.0))

    return (; model, Ip_h, gp_h, I_h, g_h, s_h, U_h, S)
end

"""
    ShhcConstraint

One chance constraint of the `sHHC_SOC` program in the hosting variables g_h (= I_h²):
E[S] + λ·√Var[S] ≤ `limit`, where S has the PCE coefficients `coef[h]` (length |K|) per
unit g_h. `kind` is `:ihd` (`key` = (h, bus id), a single harmonic), `:thd` or `:rmsv`
(`key` = bus id) or `:rmsi` (`key` = branch id).
"""
struct ShhcConstraint
    kind::Symbol
    key::Any
    coef::Dict{Int,Vector{Float64}}
    limit::Float64
    λ::Float64
end

"""
    shhc_constraints(hdata, moments) -> Vector{ShhcConstraint}

The chance constraints of step 6 of [`solve_shhc`](@ref) for `sHHC_SOC`, for the network
realization whose moments are `moments`: IHD per bus and harmonic, THD and RMS voltage
per bus, RMS current per branch. The one definition of these constraints: the solver
builds them with [`add_shhc_constraint!`](@ref) and the robust oracle evaluates them with
[`constraint_value`](@ref).
"""
function shhc_constraints(hdata::Dict, moments)
    (; alpha_V, alpha_I, bus_idx, ref_bus_ids, harm_nw_strs, harm_nw_ints) = moments
    pce     = hdata["pce"]
    fund_nw = hdata["nw"]["1"]
    lambda   = pce["lambda"]
    lambdas  = get(pce, "lambdas", nothing)   # per-constraint λ* (empirical) or nothing
    # λ for a given constraint: per-constraint empirical value if available, else scalar
    λfor(kind, key) = (lambdas !== nothing && haskey(lambdas[kind], key)) ?
                       lambdas[kind][key] : lambda
    # With empirical_moments = "pce", λ is recalibrated so that the bound built on the
    # expansion's own mean and standard deviation equals the sampled (1−ε) quantile;
    # the aggregates are evaluated at the hosting decision g of the calibration.
    pce_moments = get(pce, "empirical_moments", "pce") == "pce" &&
                  lambdas !== nothing && haskey(lambdas, :q)
    function λfor(kind, key, coef)
        pce_moments && haskey(lambdas[:q][kind], key) || return λfor(kind, key)
        g  = lambdas[:g]; gv(h) = kind == :ihd ? 1.0 : get(g, h, 0.0)
        E  = sum(coef[h][1] * gv(h) for h in keys(coef))
        σ  = sqrt(sum(psi_norms[k+1] * sum(coef[h][k+1] * gv(h) for h in keys(coef))^2
                      for k in K if k != 0; init = 0.0))
        return σ > 0.0 ? max((lambdas[:q][kind][key] - E) / σ, 0.0) : 0.0
    end
    psi_norms = pce["psi_norms"]
    K = pce["K"]

    cons = ShhcConstraint[]

    # IHD — one harmonic per (bus, harmonic)
    for (nw_str, nw_int) in zip(harm_nw_strs, harm_nw_ints)
        haskey(alpha_V, nw_int) || continue
        A = alpha_V[nw_int]
        for (bi_str, bus) in hdata["nw"][nw_str]["bus"]
            haskey(bus, "ihdmax") || continue
            n      = bus_idx[bus["bus_i"]]
            v_fund = fund_nw["bus"][bi_str]["vm"]
            push!(cons, ShhcConstraint(:ihd, (nw_int, bus["bus_i"]),
                Dict(nw_int => A[n, :]), (bus["ihdmax"] * v_fund)^2,
                λfor(:ihd, (nw_int, bus["bus_i"]), Dict(nw_int => view(A, n, :)))))
        end
    end

    # THD + RMS-voltage — per bus (cross-harmonic on Σ_h|V_{n,h}|²)
    for (bi_str, fbus) in fund_nw["bus"]
        n = bus_idx[fbus["bus_i"]]
        n in ref_bus_ids && continue
        coef   = Dict(h => alpha_V[h][n, :] for h in harm_nw_ints if haskey(alpha_V, h))
        v_fund = fbus["vm"]
        λ      = λfor(:agg, fbus["bus_i"], coef)
        haskey(fbus, "thdmax")  && push!(cons, ShhcConstraint(:thd, fbus["bus_i"], coef,
                                             (fbus["thdmax"] * v_fund)^2, λ))
        haskey(fbus, "vmaxrms") && push!(cons, ShhcConstraint(:rmsv, fbus["bus_i"], coef,
                                             fbus["vmaxrms"]^2 - v_fund^2, λ))
    end

    # RMS-current — per branch (cross-harmonic on Σ_h|I_{b,h}|²)
    for (b_str, fbr) in fund_nw["branch"]
        haskey(fbr, "c_rating") || continue
        b    = parse(Int, b_str)
        coef = Dict(h => alpha_I[h][b] for h in harm_nw_ints
                    if haskey(alpha_I, h) && haskey(alpha_I[h], b))
        isempty(coef) && continue
        push!(cons, ShhcConstraint(:rmsi, b, coef,
                                   fbr["c_rating"]^2 - get(fbr, "cm_fr", 0.0)^2, λfor(:rmsi, b, coef)))
    end
    return cons
end

"""
    add_shhc_constraint!(master, c, moments)

Add the chance constraint `c` ([`ShhcConstraint`](@ref)) to `master` (from
[`shhc_master`](@ref)): a linear bound for `:ihd`, a second-order cone otherwise, both
normalised by the limit.
"""
function add_shhc_constraint!(master, c::ShhcConstraint, moments)
    (; K, psi_norms) = moments
    (; model, g_h) = master
    coef, limit, λ = c.coef, c.limit, c.λ

    if c.kind === :ihd
        # a single harmonic: E[S] + λ√Var[S] = Cscal·g_h for g_h ≥ 0, a linear bound
        h, A  = only(coef)
        Cscal = A[1] + λ * sqrt(sum(psi_norms[k+1] * A[idx]^2
                                    for (idx, k) in enumerate(K) if k != 0; init=0.0))
        Cscal > 0.0 && JuMP.@constraint(model, (Cscal / limit) * g_h[h] <= 1.0)
        return master
    end

    # Cantelli bound  E[S] + λ√Var[S] ≤ limit  for an aggregate S whose per-harmonic
    # PCE coefficients are coef[h] (length |K|); added as a second-order cone in g_h.
    Kh = [k for k in K if k != 0]              # k ≥ 1 PCE indices
    ES = JuMP.@expression(model, sum(coef[h][1] * g_h[h] for h in keys(coef)))
    if limit <= 0.0
        # The fundamental already uses the whole limit: no harmonic headroom.
        # E[S] ≤ limit forces every contributing g_h to zero (infeasible when
        # the limit is negative), exactly as in the dHHC and QCQP models.
        @warn "solve_shhc (SOCP): non-positive limit $limit " *
              "(fundamental already at/over the limit); harmonic headroom is zero."
        JuMP.@constraint(model, ES <= limit)
        return master
    end
    σ = 1.0 / limit                   # normalise the constraint by its limit
    Wk = [JuMP.@expression(model, sum(coef[h][k+1] * g_h[h] for h in keys(coef))) for k in Kh]
    if λ <= 0.0                       # empirical λ*≈0 ⇒ mean bound (no variance term)
        JuMP.@constraint(model, σ * ES <= 1.0); return master
    end
    JuMP.@constraint(model,
        vcat(σ * (limit - ES) / λ,
             [σ * sqrt(psi_norms[k+1]) * Wk[j] for (j, k) in enumerate(Kh)])
        in JuMP.SecondOrderCone())
    return master
end

"""
    add_shhc_constraints!(master, hdata, moments)

Add every chance constraint of [`shhc_constraints`](@ref) for the network realization
whose moments are `moments` to `master`.
"""
function add_shhc_constraints!(master, hdata::Dict, moments)
    for c in shhc_constraints(hdata, moments)
        add_shhc_constraint!(master, c, moments)
    end
    return master
end

"""
    constraint_value(c, g, moments) -> Float64

E[S] + λ·√Var[S] of the chance constraint `c` at the hosting variables `g` (h ⇒ g_h,
with g_h = I_h² the squared hosting current), to compare with `c.limit`.
"""
function constraint_value(c::ShhcConstraint, g, moments)
    (; K, psi_norms) = moments
    ES = sum(c.coef[h][1] * g[h] for h in keys(c.coef))
    c.λ <= 0.0 && return ES
    var = sum(psi_norms[k+1] * sum(c.coef[h][k+1] * g[h] for h in keys(c.coef))^2
              for k in K if k != 0; init = 0.0)
    return ES + c.λ * sqrt(var)
end

"""
    normalized_violation(c, g, moments) -> Float64

`constraint_value(c, g, moments) / c.limit − 1`: positive when `c` is violated at `g`, in
the units of the normalised constraint of the program. For a non-positive limit, where the
program forces the contributing g_h to zero, it returns `constraint_value − limit`.
"""
function normalized_violation(c::ShhcConstraint, g, moments)
    v = constraint_value(c, g, moments)
    return c.limit > 0.0 ? v / c.limit - 1.0 : v - c.limit
end

"""
    solve_shhc(hdata, sHHC_SOC, optimizer; pce_params, hpf_optimizer, kwargs...)

Solve the stochastic HHC problem as a Second-Order Cone Program (SOCP) for a
fairness principle with fixed allocation ratios: `"absolute equality"` (every
unit receives I_h) or `"weighted equality"` (unit d receives w_d·I_h, with w_d
from `load["fairness_weight"]`, see [`fairness_weights`](@ref)).

With fixed ratios, one scalar hosting capacity I_h per harmonic remains, so the squared bus voltages and branch currents are
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
`"solution"` (Dict mapping harmonic nw id → optimal I_h; a participating unit d hosts
w_d·I_h). At a zero-sequence harmonic, floating units inject zero and are omitted
from its moments, ratings and objective weights; only an entirely floating set of
harmonic units forces I_h to zero.

Only principles with fixed allocation ratios are representable in this
formulation: the single scalar I_h per harmonic is what makes the squared
magnitudes quadratic in one variable and the chance constraints conic. Any other
`hdata["principle"]` is rejected rather than silently ignored; use `sHHC_QCQP`,
which carries a per-unit current variable.
"""
function solve_shhc(hdata::Dict, ::Type{sHHC_SOC}, optimizer;
                    pce_params::Dict, hpf_optimizer=optimizer, kwargs...)
    principle = get(hdata, "principle", "absolute equality")
    principle in ("absolute equality", "weighted equality") || error(
        "solve_shhc (sHHC_SOC): fairness principle \"$principle\" is not " *
        "representable in the SOCP, which fixes the ratios between the units' " *
        "allocations (I_{d,h} = w_d·I_h, one scalar I_h per harmonic). Use sHHC_QCQP " *
        "for \"$principle\", or set hdata[\"principle\"] to \"absolute equality\" " *
        "or \"weighted equality\" (weights in load[\"fairness_weight\"]).")

    # The HPF requires a nonlinear solver (e.g. Ipopt); hpf_optimizer lets
    # callers pass a different solver than the SOCP optimizer. It runs before the
    # PCE data so a local angle reference can rotate the angle supports.
    update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, hpf_optimizer)
    update_hdata_with_pce_data!(hdata, pce_params)

    moments = shhc_moments(hdata)
    master  = shhc_master(hdata, moments, optimizer)
    cons = shhc_constraints(hdata, moments)
    used = Dict(:ihd => Dict{Any,Float64}(), :agg => Dict{Any,Float64}(), :rmsi => Dict{Any,Float64}())
    for c in cons
        add_shhc_constraint!(master, c, moments)
        kind = c.kind in (:thd, :rmsv) ? :agg : c.kind
        used[kind][c.key] = c.λ
    end
    (; model, I_h, S) = master
    (; alpha_re, alpha_im, alpha_V, alpha_I, harm_nw_ints) = moments
    pce = hdata["pce"]

    JuMP.optimize!(model)

    # cache precomputed scalars for external result reporting
    pce["alpha_re"] = alpha_re; pce["alpha_im"] = alpha_im
    pce["alpha_V"]  = alpha_V;  pce["alpha_I"]  = alpha_I

    return Dict(
        "termination_status" => JuMP.termination_status(model),
        "objective_value"    => JuMP.has_values(model) ? S * JuMP.objective_value(model) : NaN,
        "solve_time"         => JuMP.solve_time(model),
        "lambdas_used"       => used,     # λ per constraint as applied (:ihd, :agg, :rmsi)
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
    # workspace with PCE data (angle ranges, ε, α̂β̂) for λ* sampling. The
    # fundamental power flow runs first, as in solve_shhc: with the local angle
    # reference the sampled angles must include the offsets h·θ_i it produces.
    hd_pce  = deepcopy(hdata)
    hpf_opt = get(kwargs, :hpf_optimizer, M === sHHC_QCQP ? optimizer : nothing)
    if hpf_opt !== nothing
        update_hdata_with_fundamental_hpf_results!(hd_pce, dHHC_NLP, hpf_opt)
    elseif get(hdata, "angle_reference", "local") == "local"
        @warn "solve_shhc_calibrated: no hpf_optimizer, so λ* is sampled without the " *
              "local angle offsets." maxlog=1
    end
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
