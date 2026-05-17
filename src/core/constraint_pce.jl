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

# ── 6.3: Galerkin projection for voltage magnitude squared ────────────────────
# Exact equality (not SOC relaxation): ξ_k == Σ_{k1,k2} T[k1,k2,k]·(VR_{k1}·VR_{k2} + VI_{k1}·VI_{k2})

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
            VR1 = _PMs.var(pm, nw1, :vr, n_id)
            VR2 = _PMs.var(pm, nw2, :vr, n_id)
            VI1 = _PMs.var(pm, nw1, :vi, n_id)
            VI2 = _PMs.var(pm, nw2, :vi, n_id)
            JuMP.add_to_expression!(expr, m_val, VR1, VR2)
            JuMP.add_to_expression!(expr, m_val, VI1, VI2)
        end
        JuMP.@constraint(pm.model, expr == ξ_k)
    end
end

# ── 6.4: IHD chance constraint ────────────────────────────────────────────────
# Squared Cantelli (no sigma variable):
#   E[ξ_h] ≤ limit²  and  λ² · Var[ξ_h] ≤ (limit² - E[ξ_h])²
# This avoids the SOC auxiliary variable whose Jacobian is zero at the
# trivial point, preventing LICQ failure in Ipopt.

function constraint_chance_ihd(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_idx       :: Int,
    n_id        :: Any,
    :: Int
)
    nw_0      = nw_id(h_idx, 0, pce.P_size)
    ξ_0       = _PMs.var(pm, nw_0, :xi)[n_id]

    bus_data  = _PMs.ref(pm, nw_0, :bus, n_id)
    ihd_limit = bus_data["ihdmax"]
    v_nom     = _PMs.ref(pm, fundamental(pm), :bus, n_id)["vm"]
    λ         = pm.data["sdata"]["lambda"]
    limit_sq  = (ihd_limit * v_nom)^2

    xi_higher = [_PMs.var(pm, nw_id(h_idx, k, pce.P_size), :xi)[n_id]
                 for k in 1:pce.deg]
    norms_k   = pce.norms[2:end]

    JuMP.@constraint(pm.model, ξ_0 <= limit_sq)
    JuMP.@constraint(pm.model,
        λ^2 * sum(norms_k[i] * xi_higher[i]^2 for i in eachindex(norms_k))
        <= (limit_sq - ξ_0)^2)
end

# ── 6.5: THD chance constraint ────────────────────────────────────────────────
# Squared Cantelli (no sigma variable):
#   E[Σ_h ξ_h] ≤ limit²  and  λ² · Var[Σ_h ξ_h] ≤ (limit² - E[Σ_h ξ_h])²

function constraint_chance_thd(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_indices   :: Vector{Int},
    n_id        :: Any,
    :: Int
)
    thd_limit   = _PMs.ref(pm, fundamental(pm), :bus, n_id)["thdmax"]
    v_nom       = _PMs.ref(pm, fundamental(pm), :bus, n_id)["vm"]
    λ           = pm.data["sdata"]["lambda"]
    limit_sq    = (thd_limit * v_nom)^2

    xi_mean_sum = sum(_PMs.var(pm, nw_id(h, 0, pce.P_size), :xi)[n_id]
                      for h in h_indices)
    xi_hi_sum   = [sum(_PMs.var(pm, nw_id(h, k, pce.P_size), :xi)[n_id]
                       for h in h_indices)
                   for k in 1:pce.deg]
    norms_k     = pce.norms[2:end]

    JuMP.@constraint(pm.model, xi_mean_sum <= limit_sq)
    JuMP.@constraint(pm.model,
        λ^2 * sum(norms_k[k] * xi_hi_sum[k]^2 for k in eachindex(norms_k))
        <= (limit_sq - xi_mean_sum)^2)
end

# ── 6.6: Galerkin projection for branch current squared ───────────────────────
# Exact equality: J_branch_k == Σ_{k1,k2} T[k1,k2,k]·(CR_{k1}·CR_{k2} + CI_{k1}·CI_{k2})

function constraint_pce_soc_branch_current(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    nw          :: Int,
    h_idx       :: Int,
    :: Int
)
    k = get_pce_mode(nw, pce.P_size)
    for (b, branch) in _PMs.ref(pm, nw, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (b, f_bus, t_bus)
        J_k   = _PMs.var(pm, nw, :J_branch)[b]

        expr = JuMP.QuadExpr()
        for k1 in 0:pce.deg, k2 in 0:pce.deg
            m_val = pce.T[k1+1, k2+1, k+1]
            abs(m_val) < 1e-14 && continue
            nw1 = nw_id(h_idx, k1, pce.P_size)
            nw2 = nw_id(h_idx, k2, pce.P_size)
            CR_1 = _PMs.var(pm, nw1, :cr, f_idx)
            CR_2 = _PMs.var(pm, nw2, :cr, f_idx)
            CI_1 = _PMs.var(pm, nw1, :ci, f_idx)
            CI_2 = _PMs.var(pm, nw2, :ci, f_idx)
            JuMP.add_to_expression!(expr, m_val, CR_1, CR_2)
            JuMP.add_to_expression!(expr, m_val, CI_1, CI_2)
        end
        JuMP.@constraint(pm.model, expr == J_k)
    end
end

# ── 6.8: RMS voltage chance constraint ────────────────────────────────────────
# Squared Cantelli (no sigma variable):
#   vm_1² + E[Σ_h ξ_h] ≤ vmaxrms²  and  λ² · Var[Σ_h ξ_h] ≤ (vmaxrms² - vm_1² - E[...])²

function constraint_chance_rms_voltage(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_indices   :: Vector{Int},
    n_id        :: Any,
    :: Int
)
    bus_data   = _PMs.ref(pm, fundamental(pm), :bus, n_id)
    vmaxrms    = bus_data["vmaxrms"]
    vm_1       = bus_data["vm"]
    λ          = pm.data["sdata"]["lambda"]
    slack_const = vmaxrms^2 - vm_1^2

    xi_mean    = [_PMs.var(pm, nw_id(h, 0, pce.P_size), :xi)[n_id] for h in h_indices]
    xi_hi_sum  = [sum(_PMs.var(pm, nw_id(h, k, pce.P_size), :xi)[n_id] for h in h_indices)
                  for k in 1:pce.deg]
    norms_k    = pce.norms[2:end]
    xi_mean_sum = sum(xi_mean)

    JuMP.@constraint(pm.model, xi_mean_sum <= slack_const)
    JuMP.@constraint(pm.model,
        λ^2 * sum(norms_k[k] * xi_hi_sum[k]^2 for k in eachindex(norms_k))
        <= (slack_const - xi_mean_sum)^2)
end

# ── 6.9: RMS branch current chance constraint ─────────────────────────────────
# Squared Cantelli (no sigma variable):
#   cm_fund² + E[Σ_h J_h] ≤ c_rating²  and  λ² · Var[Σ_h J_h] ≤ (c_rating² - cm_fund² - E[...])²

function constraint_chance_rms_current(
    pm          :: AbstractSHHCModel,
    pce         :: PCEData,
    h_indices   :: Vector{Int},
    b_id        :: Any,
    :: Int
)
    branch     = _PMs.ref(pm, fundamental(pm), :branch, b_id)
    c_rating   = branch["c_rating"]
    cm_fund    = get(branch, "cm_fr", 0.0)
    λ          = pm.data["sdata"]["lambda"]
    slack_const = c_rating^2 - cm_fund^2

    J_mean     = [_PMs.var(pm, nw_id(h, 0, pce.P_size), :J_branch)[b_id] for h in h_indices]
    J_hi_sum   = [sum(_PMs.var(pm, nw_id(h, k, pce.P_size), :J_branch)[b_id] for h in h_indices)
                  for k in 1:pce.deg]
    norms_k    = pce.norms[2:end]
    J_mean_sum = sum(J_mean)

    JuMP.@constraint(pm.model, J_mean_sum <= slack_const)
    JuMP.@constraint(pm.model,
        λ^2 * sum(norms_k[k] * J_hi_sum[k]^2 for k in eachindex(norms_k))
        <= (slack_const - J_mean_sum)^2)
end
