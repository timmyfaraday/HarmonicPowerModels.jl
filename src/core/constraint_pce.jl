################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.1 - added PCE constraint functions for sHHC                             #
################################################################################

# ── 6.2: Angle-to-rectangular PCE constraint ─────────────────────────────────
# Paper Eqs.(30)-(31): crd_k = cos_k * I_budget,  cid_k = sin_k * I_budget
# for all PCE modes k = 0, …, deg.

function constraint_pce_angle_to_rectangular(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_idx       :: Int,
    u_id        :: Any,
    :: Int
)
    nw_0 = nw_id(h_idx, 0, pce.P_size)

    # Read stochastic angle range for this load and harmonic.
    # :load is the harmonic unit component (AUDIT.md H).
    u_sdata    = _PMs.ref(pm, nw_0, :load, u_id, "sdata")
    h_str      = string(h_idx)
    v_min      = u_sdata[h_str]["angle_min_deg"]
    v_max      = u_sdata[h_str]["angle_max_deg"]

    phi_coeffs = pce_angle_coefficients(pce, Float64(v_min), Float64(v_max))
    cos_k      = galerkin_cos(pce, phi_coeffs)
    sin_k      = galerkin_sin(pce, phi_coeffs)

    # :cmd is the hosting-capacity budget (AUDIT.md I); lives only in k=0 network.
    I_budget = _PMs.var(pm, nw_0, :cmd, u_id)

    for k in 0:pce.deg
        nw   = nw_id(h_idx, k, pce.P_size)
        # :crd / :cid are the real/imaginary load current symbols (AUDIT.md I).
        IR_k = _PMs.var(pm, nw, :crd, u_id)
        II_k = _PMs.var(pm, nw, :cid, u_id)
        JuMP.@constraint(pm.model, cos_k[k+1] * I_budget == IR_k)
        JuMP.@constraint(pm.model, sin_k[k+1] * I_budget == II_k)
    end
end

# ── 6.3: SOC voltage magnitude with PCE multiplication ───────────────────────
# Paper Eq.(42): Σ_{k1,k2} T[k1,k2,k] · (VR_{k1}·VR_{k2} + VI_{k1}·VI_{k2}) ≤ ξ_k

function constraint_pce_soc_voltage(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    nw          :: Int,
    h_idx       :: Int,
    :: Int
)
    k = get_pce_mode(nw, pce.P_size)
    for n_id in _PMs.ids(pm, nw, :bus)
        ξ_k  = _PMs.var(pm, nw, :xi)[n_id]
        expr = JuMP.QuadExpr()
        for k1 in 0:pce.deg, k2 in 0:pce.deg
            m_val = pce.T[k1+1, k2+1, k+1]
            abs(m_val) < 1e-14 && continue
            nw1 = nw_id(h_idx, k1, pce.P_size)
            nw2 = nw_id(h_idx, k2, pce.P_size)
            # :vr / :vi are the bus voltage symbols (AUDIT.md I).
            VR1 = _PMs.var(pm, nw1, :vr, n_id)
            VR2 = _PMs.var(pm, nw2, :vr, n_id)
            VI1 = _PMs.var(pm, nw1, :vi, n_id)
            VI2 = _PMs.var(pm, nw2, :vi, n_id)
            JuMP.add_to_expression!(expr, m_val, VR1, VR2)
            JuMP.add_to_expression!(expr, m_val, VI1, VI2)
        end
        JuMP.@constraint(pm.model, expr <= ξ_k)
    end
end

# ── 6.4: IHD chance constraint ────────────────────────────────────────────────
# Cantelli bound: E[ξ_h] + λ · σ_h ≤ (ihdmax · v_nom)²
# SOC encodes σ_h ≥ ||norms .* ξ_higher||

function constraint_chance_ihd(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_idx       :: Int,
    n_id        :: Any,
    :: Int
)
    nw_0      = nw_id(h_idx, 0, pce.P_size)
    ξ_0       = _PMs.var(pm, nw_0, :xi)[n_id]
    σ         = _PMs.var(pm, nw_0, :sigma_ihd)[n_id]

    # IHD limit from AUDIT.md H: ref(pm, nw, :bus, i, "ihdmax").
    bus_data  = _PMs.ref(pm, nw_0, :bus, n_id)
    ihd_limit = bus_data["ihdmax"]
    # Fundamental bus voltage from AUDIT.md H: ref(pm, fundamental(pm), :bus, i, "vm").
    v_nom     = _PMs.ref(pm, fundamental(pm), :bus, n_id)["vm"]
    # λ from global sdata (set by compute_lambda! in data_pce.jl).
    λ         = pm.data["sdata"]["lambda"]

    xi_higher = [_PMs.var(pm, nw_id(h_idx, k, pce.P_size), :xi)[n_id]
                 for k in 1:pce.deg]
    norms_k   = pce.norms[2:end]

    JuMP.@constraint(pm.model,
        [σ; norms_k .* xi_higher] in JuMP.SecondOrderCone())
    JuMP.@constraint(pm.model, ξ_0 + λ * σ <= (ihd_limit * v_nom)^2)
end

# ── 6.5: THD chance constraint ────────────────────────────────────────────────
# Cantelli bound: E[Σ_h ξ_h] + λ · σ_thd ≤ (thdmax · v_nom)²

function constraint_chance_thd(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_indices   :: Vector{Int},
    n_id        :: Any,
    :: Int
)
    nw_0_ref    = nw_id(h_indices[1], 0, pce.P_size)
    σ_thd       = _PMs.var(pm, nw_0_ref, :sigma_thd)[n_id]
    # THD limit from AUDIT.md H: ref(pm, fundamental(pm), :bus, i, "thdmax").
    thd_limit   = _PMs.ref(pm, fundamental(pm), :bus, n_id)["thdmax"]
    # Fundamental bus voltage from AUDIT.md H.
    v_nom       = _PMs.ref(pm, fundamental(pm), :bus, n_id)["vm"]
    # λ from global sdata.
    λ           = pm.data["sdata"]["lambda"]

    xi_mean_sum = sum(
        _PMs.var(pm, nw_id(h, 0, pce.P_size), :xi)[n_id] for h in h_indices)
    xi_hi_sum   = [sum(_PMs.var(pm, nw_id(h, k, pce.P_size), :xi)[n_id]
                       for h in h_indices)
                   for k in 1:pce.deg]
    norms_k     = pce.norms[2:end]

    JuMP.@constraint(pm.model,
        [σ_thd; norms_k .* xi_hi_sum] in JuMP.SecondOrderCone())
    JuMP.@constraint(pm.model, xi_mean_sum + λ * σ_thd <= (thd_limit * v_nom)^2)
end
