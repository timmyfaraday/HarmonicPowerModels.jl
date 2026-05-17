################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.1 - added sHHC_SOC entry point and build function                       #
################################################################################

# ── Entry point ──────────────────────────────────────────────────────────────

"""
    solve_shhc_soc(data, optimizer; H, bus_id, kwargs...)

Solve the stochastic harmonic hosting capacity problem using a second-order
cone formulation with polynomial chaos expansion.

`data` must be a flat (single-network) PowerModels data dict augmented with:
  - `data["sdata"]` containing `"alpha"`, `"beta"`, `"epsilon"` and optionally `"deg"`, `"H"`
  - `data["load"][id]["sdata"][h_str]` containing `"angle_min_deg"`, `"angle_max_deg"`

The fundamental bus voltage `"vm"` (used by the chance constraints) must be
pre-populated in `data["bus"][id]` before calling this function (e.g., via a
prior `solve_hpf` run).  If it is absent the chance constraints will error.
"""
function solve_shhc_soc(data::Dict, optimizer;
                        H::Vector{Int} = Int[],
                        bus_id::Int    = 1,
                        kwargs...)
    check_shhc_data(data)

    # Default fairness principle to maximum efficiency if not specified.
    haskey(data, "principle") || (data["principle"] = "absolute equality")

    α   = Float64(data["sdata"]["alpha"])
    β   = Float64(data["sdata"]["beta"])
    deg = Int(get(data["sdata"], "deg", 2))
    _H  = isempty(H) ? get(data["sdata"], "H", H) : H

    pce = build_pce_data(α, β; deg = deg)
    mn  = build_mn_pce_data(data, pce; bus_id = bus_id, H = _H)

    return _PMs.solve_model(mn, sHHC_SOC, optimizer, build_shhc_soc!;
                            ref_extensions      = [ref_add_filter!, ref_add_xfmr!],
                            solution_processors = [_HPM.sol_data_model!],
                            multinetwork        = true,
                            kwargs...)
end

# ── Model builder ─────────────────────────────────────────────────────────────

function build_shhc_soc!(pm::AbstractSHHCModel)
    pce    = pm.data["pce"]
    n_harm = pm.data["n_harmonics"]

    # Sorted list of all harmonic numbers present in the mean-mode sub-networks.
    h_indices = sort(unique([
        pm.data["nw"]["$nw"]["harmonic_idx"]
        for nw in _PMs.nw_ids(pm)
        if pm.data["nw"]["$nw"]["is_mean_mode"]
    ]))
    # Non-fundamental harmonics used for angle constraints and chance constraints.
    h_harm = filter(h -> h ≠ fundamental(pm), h_indices)

    # ── Variable declarations ──────────────────────────────────────────────
    for nw in _PMs.nw_ids(pm)
        nw_data = pm.data["nw"]["$nw"]
        h_idx   = nw_data["harmonic_idx"]
        k       = nw_data["pce_mode"]

        # Fundamental is deterministic: skip h=1 higher PCE modes (k>0).
        h_idx == fundamental(pm) && k > 0 && continue

        # Bus voltage with PCE-aware bounding (:vr, :vi bounded only for mean mode).
        variable_bus_voltage_pce(pm, nw)
        # Squared voltage magnitude per PCE mode (:xi, lb=0 only for mean mode).
        variable_voltage_squared_pce(pm, nw)
        # Lifted squared bus injection current per PCE mode (:J_bus).
        variable_bus_injection_current_squared_pce(pm, nw)
        # Lifted squared branch current per PCE mode (:J_branch).
        variable_branch_current_squared_pce(pm, nw)

        # Standard network-element variables (same as dHHC_SOC).
        variable_xfmr_voltage(pm,  nw=nw, bounded=true)
        variable_branch_current(pm, nw=nw, bounded=true)
        variable_xfmr_current(pm,  nw=nw, bounded=true)
        variable_filter_current(pm, nw=nw, bounded=false)
        variable_gen_current(pm,    nw=nw, bounded=true)

        # Load real/imaginary currents declared for every PCE mode.
        variable_load_current_real(pm,      nw=nw, bounded=true)
        variable_load_current_imaginary(pm, nw=nw, bounded=true)

        # Harmonic mean-mode only: hosting-capacity budget and fairness.
        if nw_data["is_mean_mode"] && h_idx ≠ fundamental(pm)
            variable_hosting_capacity_pce(pm, nw)    # :cmd
            variable_fairness_pce(pm, nw)            # :cmh or :fh per principle
        end
    end

    # ── Objective ─────────────────────────────────────────────────────────
    objective_maximum_hosting_capacity(pm)

    # ── Per-network constraints ────────────────────────────────────────────
    for nw in _PMs.nw_ids(pm)
        nw_data = pm.data["nw"]["$nw"]
        h_idx   = nw_data["harmonic_idx"]
        k       = nw_data["pce_mode"]

        # Fundamental is deterministic: skip h=1 higher PCE modes (k>0).
        h_idx == fundamental(pm) && k > 0 && continue

        # Reference bus.
        if !haskey(pm.setting, "fix_refbus_angle") ||
                pm.setting["fix_refbus_angle"] == true
            for i in _PMs.ids(pm, :ref_buses, nw=nw)
                constraint_voltage_ref_bus(pm, i, nw=nw)
            end
        end

        # KCL: current balance at each bus.
        for i in _PMs.ids(pm, :bus, nw=nw)
            constraint_current_balance(pm, i, nw=nw)
        end

        # KVL / Ohm: branch current injections and voltage drop.
        for b in _PMs.ids(pm, :branch, nw=nw)
            _PMs.constraint_current_from(pm, b, nw=nw)
            _PMs.constraint_current_to(pm, b, nw=nw)
            _PMs.constraint_voltage_drop(pm, b, nw=nw)
        end

        # Generator injection (Thevenin equivalent shunt).
        for g in _PMs.ids(pm, :gen, nw=nw)
            constraint_gen_current(pm, g, nw=nw)
        end

        # Transformer core and winding constraints (RF-1 resolved via AbstractSHHCModel dispatch).
        for x in _PMs.ids(pm, :xfmr, nw=nw)
            constraint_xfmr_core_magnetization(pm, x, nw=nw)
            constraint_xfmr_core_voltage_drop(pm, x, nw=nw)
            constraint_xfmr_core_voltage_phase_shift(pm, x, nw=nw)
            constraint_xfmr_core_current_balance(pm, x, nw=nw)
            constraint_xfmr_winding_config(pm, x, nw=nw)
            constraint_xfmr_winding_current_balance(pm, x, nw=nw)
        end

        # Fundamental mean-mode: fix load to constant power (as in dHHC_SOC).
        if h_idx == fundamental(pm) && nw_data["is_mean_mode"]
            for l in _PMs.ids(pm, :load, nw=nw)
                constraint_load_current(pm, l, nw=nw)
            end
        end

        # PCE SOC constraints for non-fundamental harmonics only.
        if h_idx ≠ fundamental(pm)
            constraint_pce_soc_voltage(pm, pce, nw, h_idx, n_harm)
            constraint_pce_soc_current(pm, pce, nw, h_idx, n_harm)
            constraint_pce_soc_branch_current(pm, pce, nw, h_idx, n_harm)
        end
    end

    # ── Fairness principle constraints (mean-mode networks only) ──────────
    for h_idx in h_harm
        nw_0 = nw_id(h_idx, 0, pce.P_size)
        constraint_fairness_principle(pm, nw=nw_0)
    end

    # ── Angle-to-rectangular: couple :cmd budget to all PCE modes ─────────
    for h_idx in h_harm
        nw_0 = nw_id(h_idx, 0, pce.P_size)
        for u_id in _PMs.ids(pm, :load, nw=nw_0)
            constraint_pce_angle_to_rectangular(pm, pce, h_idx, u_id, n_harm)
        end
    end

    # ── Chance constraints ────────────────────────────────────────────────
    isempty(h_harm) && return
    nw_0_ref = nw_id(h_harm[1], 0, pce.P_size)

    # Reference buses have harmonic voltages fixed to zero by constraint_voltage_ref_bus,
    # so the IHD/THD/RMS constraints are trivially satisfied and degenerate there.
    # Skipping them removes the LICQ failure that prevents Ipopt from converging.
    ref_bus_ids = Set(_PMs.ids(pm, :ref_buses, nw=nw_0_ref))

    for n_id in _PMs.ids(pm, :bus, nw=nw_0_ref)
        n_id in ref_bus_ids && continue
        for h_idx in h_harm
            constraint_chance_ihd(pm, pce, h_idx, n_id, n_harm)
        end
        constraint_chance_thd(pm, pce, h_harm, n_id, n_harm)
        constraint_chance_rms_voltage(pm, pce, h_harm, n_id, n_harm)
    end

    for b_id in _PMs.ids(pm, :branch, nw=nw_0_ref)
        constraint_chance_rms_current(pm, pce, h_harm, b_id, n_harm)
    end
end
