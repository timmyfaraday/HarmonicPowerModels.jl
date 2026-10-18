################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA - this needs to be reworked asap                       #
################################################################################

""
const ihd_limits = Dict(
    "Clean Bus" =>                 [1.0, zeros(49)...],
    "IEC61000-2-4:2002, Cl. 2" =>  [1.00000, 0.02000, 0.05000, 0.01000, 0.06000, 
                                    0.00500, 0.05000, 0.00500, 0.01500, 0.00500, 
                                    0.03500, 0.00458, 0.03000, 0.00429, 0.00400, 
                                    0.00406, 0.02000, 0.00389, 0.01761, 0.00375, 
                                    0.00200, 0.00364, 0.01408, 0.00354, 0.01275, 
                                    0.00346, 0.00200, 0.00339, 0.01061, 0.00333, 
                                    0.00975, 0.00328, 0.00200, 0.00324, 0.00833, 
                                    0.00319, 0.00773, 0.00316, 0.00200, 0.00313, 
                                    0.00671, 0.00310, 0.00627, 0.00307, 0.00200, 
                                    0.00304, 0.00551, 0.00302, 0.00518, 0.00300],
    "IEC61000-3-6:2008" =>         [1.00000, 0.01400, 0.02000, 0.00800, 0.02000,
                                    0.00400, 0.02000, 0.00400, 0.01000, 0.00350,
                                    0.01500, 0.00318, 0.01500, 0.00296, 0.00300,
                                    0.00279, 0.01200, 0.00266, 0.01074, 0.00255,
                                    0.00200, 0.00246, 0.00887, 0.00239, 0.00816,
                                    0.00233, 0.00200, 0.00228, 0.00703, 0.00223,
                                    0.00658, 0.00219, 0.00200, 0.00216, 0.00583,
                                    0.00213, 0.00551, 0.00210, 0.00200, 0.00208,
                                    0.00498, 0.00205, 0.00474, 0.00203, 0.00200,
                                    0.00201, 0.00434, 0.00200, 0.00416, 0.00198],
    "AS/NZS61000-3-6" =>           [1.00000, 0.02000, 0.05000, 0.01000, 0.05000,
                                    0.00500, 0.05000, 0.00500, 0.01500, 0.00500,
                                    0.03500, 0.00200, 0.03000, 0.00200, 0.00300,
                                    0.00200, 0.02000, 0.00200, 0.01500, 0.00200,
                                    0.00200, 0.00200, 0.01500, 0.00200, 0.01500,
                                    0.00200, 0.00200, 0.00200, 0.01483, 0.00200,
                                    0.01087, 0.00200, 0.00200, 0.00200, 0.00986,
                                    0.00200, 0.00943, 0.00200, 0.00200, 0.00200,
                                    0.00871, 0.00200, 0.00840, 0.00200, 0.00200,
                                    0.00200, 0.00785, 0.00200, 0.00761, 0.00200],
    "IEEE519-2022-1/69kV" =>       [1.0, 0.030 .* ones(49)...],
    "IEEE519-2022-69/161kV" =>     [1.0, 0.015 .* ones(49)...])

const thd_limits = Dict(
    "Clean Bus" =>                  0.00000,
    "IEC61000-2-4:2002, Cl. 2" =>   0.08000,
    "IEC61000-3-6:2008" =>          0.08000,
    "AS/NZS61000-3-6" =>            0.08000,
    "IEEE519-2022-1/69kV" =>        0.05000,
    "IEEE519-2022-69/161kV" =>      0.02500) 

"""
    HarmonicPowerModels.parse_file(path; kwargs...)

Parse a MATPOWER/PowerModels case while preserving the HarmonicPowerModels
"clean bus" convention (`bus_type == 4`).

A clean bus is an ideal source: it acts as a voltage reference at the fundamental
and is forced to V=0 (ideal ground) at every harmonic.  MATPOWER reserves
`type 4` for isolated buses and PowerModels' validation rewrites/deactivates such
buses, so this wrapper parses with `validate=false`, records the marking in
`bus["is_clean"]`, remaps the clean buses to `bus_type = 3` (fundamental
reference), and only then runs the standard network corrections.

All other buses receive `bus["is_clean"] = false`.
"""
function parse_file(path::String; kwargs...)
    data = _PMs.parse_file(path; validate=false, kwargs...)
    for (_, bus) in data["bus"]
        if bus["bus_type"] == 4
            bus["is_clean"]  = true
            bus["bus_type"]  = 3        # ideal source ⇒ fundamental reference
        else
            bus["is_clean"]  = false
        end
    end
    _PMs.correct_network_data!(data)
    return data
end

""
function ihd_current_limit_ieee_519(voltage, i_ratio, harmonic)
    if voltage < 69.0
        if i_ratio == Inf
            return 0.0
        elseif i_ratio < 20
            if 2 <= harmonic < 11
                return 0.04
            elseif 11 <= harmonic < 17
                return 0.02
            elseif 17 <= harmonic < 23
                return 0.015
            elseif 23 <= harmonic < 35
                return 0.006
            elseif 35 <= harmonic <= 50
                return 0.003
            else
                return 0.0
            end
        elseif 20.0 <= i_ratio < 50.0
            if 2 <= harmonic < 11
                return 0.07
            elseif 11 <= harmonic < 17
                return 0.035
            elseif 17 <= harmonic < 23
                return 0.025
            elseif 23 <= harmonic < 35
                return 0.010
            elseif 35 <= harmonic <= 50
                return 0.005
            else
                return 0.0
            end
        elseif 50.0 <= i_ratio < 100.0
            if 2 <= harmonic < 11
                return 0.10
            elseif 11 <= harmonic < 17
                return 0.045
            elseif 17 <= harmonic < 23
                return 0.040
            elseif 23 <= harmonic < 35
                return 0.015
            elseif 35 <= harmonic <= 50
                return 0.007
            else
                return 0.0
            end
        elseif 100.0 <= i_ratio < 1000.0
            if 2 <= harmonic < 11
                return 0.12
            elseif 11 <= harmonic < 17
                return 0.055
            elseif 17 <= harmonic < 23
                return 0.050
            elseif 23 <= harmonic < 35
                return 0.020
            elseif 35 <= harmonic <= 50
                return 0.010
            else
                return 0.0
            end
        else
            if 2 <= harmonic < 11
                return 0.15
            elseif 11 <= harmonic < 17
                return 0.07
            elseif 17 <= harmonic < 23
                return 0.06
            elseif 23 <= harmonic < 35
                return 0.025
            elseif 35 <= harmonic <= 50
                return 0.014
            else
                return 0.0
            end
        end
    else
        if i_ratio == Inf
            return 0.0
        elseif i_ratio < 20       
            if 2 <= harmonic < 11
                return 0.02
            elseif 11 <= harmonic < 17
                return 0.01
            elseif 17 <= harmonic < 23
                return 0.0075
            elseif 23 <= harmonic < 35
                return 0.003
            elseif 35 <= harmonic <= 50
                return 0.0015
            else
                return 0.0
            end
        elseif 20.0 <= i_ratio < 50.0
            if 2 <= harmonic < 11
                return 0.035
            elseif 11 <= harmonic < 17
                return 0.0175
            elseif 17 <= harmonic < 23
                return 0.0125
            elseif 23 <= harmonic < 35
                return 0.005
            elseif 35 <= harmonic <= 50
                return 0.0025
            else
                return 0.0
            end
        elseif 50.0 <= i_ratio < 100.0
            if 2 <= harmonic < 11
                return 0.05
            elseif 11 <= harmonic < 17
                return 0.0225
            elseif 17 <= harmonic < 23
                return 0.02
            elseif 23 <= harmonic < 35
                return 0.0075
            elseif 35 <= harmonic <= 50
                return 0.0035
            else
                return 0.0
            end
        elseif 100.0 <= i_ratio < 1000.0
            if 2 <= harmonic < 11
                return 0.06
            elseif 11 <= harmonic < 17
                return 0.0275
            elseif 17 <= harmonic < 23
                return 0.025
            elseif 23 <= harmonic < 35
                return 0.01
            elseif 35 <= harmonic <= 50
                return 0.005
            else
                return 0.0
            end
        else
            if 2 <= harmonic < 11
                return 0.075
            elseif 11 <= harmonic < 17
                return 0.035
            elseif 17 <= harmonic < 23
                return 0.03
            elseif 23 <= harmonic < 35
                return 0.0125
            elseif 35 <= harmonic <= 50
                return 0.007
            else
                return 0.0
end end end end

# Internal: resolve the (v_min, v_max) angle range for harmonic unit `uid` at
# harmonic `h` from a `pce_params["harmonics"][h]` entry that is either a bare
# `(v_min, v_max)` tuple (shared by all units) or a `Dict{unit_id => (v_min, v_max)}`.
function _resolve_angle_range(entry, uid::Int)
    if entry isa AbstractDict
        haskey(entry, uid) || error("update_hdata_with_pce_data!: no angle range " *
                                    "for harmonic unit $uid in per-unit harmonics entry")
        return entry[uid]
    else
        return entry            # bare (v_min, v_max), shared across units
    end
end

"""
    update_hdata_with_pce_data!(hdata, pce_params)

Precompute all polynomial chaos expansion (PCE) quantities required by
`build_shhc` and store them in `hdata["pce"]`, making them accessible inside
the PowerModels solve loop as `pm.data["pce"]`.

# Model — independent germ per harmonic unit
Each *harmonic unit* (a load designated as a harmonic current source) carries
its **own** independent Beta(α̂, β̂) germ; the germ space has dimension
`N` = number of harmonic units and the multivariate basis has
`P = C(N+δ, δ)` terms. A unit's germ drives its angle at every harmonic, while
the angle *range* may differ per `(unit, harmonic)`. Loads not designated as
harmonic units are "normal" loads: they receive zero trig PCE coefficients
(no harmonic injection) and occupy no germ dimension.

# Arguments
- `hdata`      — multinetwork data dict, as returned by `replicate`.
- `pce_params` — Dict with the following keys:

| Key | Type | Description |
|---|---|---|
| `"alpha_hat"` | Float64 | Beta shape parameter α̂ |
| `"beta_hat"`  | Float64 | Beta shape parameter β̂ |
| `"delta"`     | Int     | PCE degree (validated only for δ=2) |
| `"n_terms"`   | Int     | Accepted and ignored (superseded Taylor trig PCE; see `affine_trig_pce`) |
| `"epsilon"`   | Float64 | Risk level ε for the chance-constraint bound |
| `"bound"`     | String (optional) | Chance-constraint bound: `"cantelli"` (default, distribution-free √((1-ε)/ε)) or `"gaussian"` (Φ⁻¹(1-ε)) |
| `"harmonics"` | Dict{Int, T} | Maps harmonic order h to its angle range in **radians**, where `T` is either a bare `(v_min, v_max)` tuple (shared by all units) or a `Dict{unit_id => (v_min, v_max)}` (per-unit) |
| `"angle_centre"` | String (optional) | Where each unit's angle support is centred: `"range"` (default) keeps the supplied range as given, so every unit shares one mean angle; `"ref_angle"` re-centres it on each unit's own `load["ref_angle"]` (see below) |
| `"harmonic_units"` | Vector{Int} (optional) | Load ids designated as harmonic sources; defaults to **all** loads |
| `"type_germs"` | Dict{String, Dict} (optional) | Germs shared by groups of units (see below): name → `Dict("units" => Vector{Int}, "shift" => Dict{Int => (s_min, s_max)})`, shift in **radians** per harmonic |

# Shared type germs
A group of units (a converter type) can share one extra germ ω_t, which adds a
common-mode angle shift to every unit of the group:

    φ_{u,h} = v_min + (v_max − v_min)·ω_u + s_min + (s_max − s_min)·ω_t

e.g., a firing delay common to all phase-controlled units, shifting harmonic h
by −h·α_f. Each type adds one germ (dimension N + t), so the basis grows from
C(N+δ, δ) to C(N+T+δ, δ); a unit belongs to at most one type. Units outside any
type keep their single own germ.

# Per-unit mean angles
With `"angle_centre" => "ref_angle"` the supplied range supplies only the
*dispersion*, and each unit's mean angle becomes its own `load["ref_angle"]` —
the rotation between that unit and the network reference, which `replicate`
resolves per (load, harmonic) from `bus_harmonics`, with the sign set by the
harmonic's sequence (positive keeps it, negative flips it, zero-sequence is 0).
Mean angles that differ across units are the mechanism by which contributions at
a shared bus stop being collinear, so this is what lets the aggregate distortion
fall below the arithmetic sum of its parts.

A unit with no `ref_angle` in the data is centred on 0, and a network whose
`bus_harmonics` gives every bus the same reference angle gains no diversity from
this setting — check the data before expecting an effect. For centres that do not
come from `bus_harmonics` (a converter topology or pulse number, say), build the
per-unit ranges directly with [`recentre_angle_ranges`](@ref) and pass them in the
per-unit `Dict` form instead.

# Stored result — `hdata["pce"]`

| Key | Content |
|---|---|
| `"basis"`          | multivariate `PCEBasis` (N germs, monic Beta01, probability measure) |
| `"M"`              | P×P×P multiplication tensor (dense) |
| `"lambda"`         | Cantelli safety factor λ(ε) = √((1-ε)/ε) |
| `"psi_norms"`      | Vector of ⟨Ψ_k,Ψ_k⟩, length P; `psi_norms[1] == 1` |
| `"K"`              | full PCE index range 0:(P-1) (lifted variables, chance constraints) |
| `"K_phys"`         | reduced single-germ index values (physical voltage/current variables) |
| `"cos_phi"`        | `Dict{Int,Dict{Int,Vector}}` indexed by [nw_id][load_id]; length-P per load |
| `"sin_phi"`        | `Dict{Int,Dict{Int,Vector}}` indexed by [nw_id][load_id]; length-P per load |
| `"harmonic_units"` | sorted Vector of harmonic-unit load ids (germ dimension = position) |
| `"germ_dim"`       | `Dict{Int,Int}` mapping harmonic-unit id → germ dimension 1:N |
| `"type_dim"`       | `Dict{String,Int}` mapping type name → germ dimension N+1:N+T |
| `"unit_type"`      | `Dict{Int,String}` mapping a typed unit → its type |
| `"type_shifts"`    | `Dict{Int,Dict{String,Tuple}}` per-harmonic type shift (s_min, s_max) [rad] |

The basis is orthogonal w.r.t. the Beta **probability** measure, so
`psi_norms[1] == 1` and the consumer formulas E[X] = x₀ and
Var[X] = Σ_{k≥1} psi_norms[k+1]·x_k² hold without total-measure corrections.

Each harmonic unit `d` occupies the **same** germ dimension `germ_dim[d]` in
every harmonic network, so a unit's harmonics are perfectly correlated (one
germ) while distinct units are independent. This consistency is what makes the
cross-harmonic THD/RMS chance constraints (which sum the k-th PCE coefficient
over harmonics) well-defined.

All Dict keys use the same integer types that PowerModels uses internally for
network and load ids (matching `nw::Int` and `d::Int` in the constraint functions).
"""
function update_hdata_with_pce_data!(hdata::Dict, pce_params::Dict)
    α̂       = pce_params["alpha_hat"]
    β̂       = pce_params["beta_hat"]
    δ       = pce_params["delta"]
    ε       = pce_params["epsilon"]
    h_bounds = pce_params["harmonics"]      # Dict{Int => range} in radians

    δ != 2 && @warn "update_hdata_with_pce_data!: δ=$δ detected; " *
                    "PCE constraints are validated only for δ=2."

    # ── 1. Harmonic-unit designation and germ-dimension map (D3/D6) ───────────
    all_loads = sort(unique(parse(Int, d)
                            for (nw_str, ntw) in hdata["nw"] if nw_str != "1"
                            for d in keys(ntw["load"])))
    harmonic_units = sort(collect(get(pce_params, "harmonic_units", all_loads)))
    germ_dim       = Dict(uid => i for (i, uid) in enumerate(harmonic_units))
    N              = length(harmonic_units)
    N ≥ 1 || error("update_hdata_with_pce_data!: no harmonic units designated")

    # shared type germs: dimensions N+1 … N+T, one per type (sorted by name)
    type_germs = get(pce_params, "type_germs", Dict{String,Any}())
    type_names = sort(collect(keys(type_germs)))
    type_dim   = Dict(t => N + i for (i, t) in enumerate(type_names))
    unit_type  = Dict{Int,String}()
    for t in type_names, d in type_germs[t]["units"]
        haskey(germ_dim, d) || error("update_hdata_with_pce_data!: type \"$t\" lists " *
                                     "unit $d, which is not a harmonic unit")
        haskey(unit_type, d) && error("update_hdata_with_pce_data!: unit $d belongs " *
                                      "to types \"$(unit_type[d])\" and \"$t\"")
        unit_type[d] = t
    end
    Ngerm = N + length(type_names)

    # ── 2. Multivariate PCE workspace (N+T germs, monic Beta01, prob. measure) ──
    pb = PCEBasis(Ngerm, α̂, β̂, δ)
    P  = pb.dim

    # ── 3. Multiplication tensor (dense; consumers index M[j,l,k] directly) ───
    M = dense_tensor(pb)

    # ── 4. Chance-constraint factor and basis norms (psi_norms[1] == 1) ───────
    # bound ∈ {"cantelli" (default, distribution-free), "gaussian" (N(0,1)),
    #          "empirical" (per-constraint λ* calibrated from the sampled output)}.
    # For :empirical the scalar λ is a fallback (=1.0); per-constraint factors are
    # supplied via pce_params["empirical_lambdas"] and stored as pce["lambdas"].
    bound  = Symbol(get(pce_params, "bound", "cantelli"))
    lambda = risk_factor(ε, bound)

    # Where each unit's angle support is centred (see the loop below).
    centre_mode = get(pce_params, "angle_centre", "range")
    centre_mode in ("range", "ref_angle") || error(
        "update_hdata_with_pce_data!: angle_centre = \"$centre_mode\" is not " *
        "recognised; use \"range\" (default) or \"ref_angle\".")
    Eω = α̂ / (α̂ + β̂)          # Beta mean fraction, for re-centring

    emp_lambdas = get(pce_params, "empirical_lambdas", nothing)
    psi_norms = pb.sp2

    # ── 5. Index sets: full (lifted vars) and reduced single-germ (physical) ──
    # Physical voltage/current PCE coefficients are nonzero only on multi-indices
    # supported on a single germ (each is linear in one unit's cos/sin); cross-
    # germ indices are forced to zero by KCL. The lifted |·|² variables, however,
    # populate the full basis (products of two single-germ terms). Exposing the
    # reduced set K_phys lets the QCQP create physical variables only where they
    # can be nonzero (0-based index values, a subset of K).
    # With type germs a unit's cos/sin also live on (own germ, type germ) pairs,
    # so those mixed indices are physical too.
    K      = 0:(P - 1)
    pairs  = Set((germ_dim[d], type_dim[t]) for (d, t) in unit_type)
    function _physical(k)
        nz = findall(!iszero, view(pb.ind, k, :))
        length(nz) <= 1 || (length(nz) == 2 && (nz[1], nz[2]) in pairs)
    end
    K_phys = [k - 1 for k in 1:P if _physical(k)]

    cos_phi = Dict{Int, Dict{Int, Vector{Float64}}}()
    sin_phi = Dict{Int, Dict{Int, Vector{Float64}}}()
    # per-unit, per-harmonic angle support [rad] — retained for empirical-λ sampling
    angle_ranges = Dict{Int, Dict{Int, Tuple{Float64,Float64}}}()
    # per-harmonic type shift [rad] — retained for empirical-λ sampling
    type_shifts  = Dict{Int, Dict{String, Tuple{Float64,Float64}}}()
    zerovec = zeros(P)   # normal loads: no harmonic injection

    for (nw_str, ntw) in hdata["nw"]
        nw_int = parse(Int, nw_str)

        # Fundamental network has deterministic voltage — no stochastic angle.
        nw_int == 1 && continue

        # Skip harmonics for which no angle distribution was provided.
        haskey(h_bounds, nw_int) || continue
        entry = h_bounds[nw_int]

        cos_phi[nw_int] = Dict{Int, Vector{Float64}}()
        sin_phi[nw_int] = Dict{Int, Vector{Float64}}()
        angle_ranges[nw_int] = Dict{Int, Tuple{Float64,Float64}}()
        type_shifts[nw_int] = Dict{String, Tuple{Float64,Float64}}(
            t => Tuple(Float64.(type_germs[t]["shift"][nw_int])) for t in type_names)

        for d_str in keys(ntw["load"])
            d_int = parse(Int, d_str)
            if haskey(germ_dim, d_int)
                # harmonic unit: own germ dimension, support-centred trig PCE
                v_min, v_max = _resolve_angle_range(entry, d_int)

                # Per-unit mean angle. With centre_mode == "ref_angle" the support
                # is re-centred on this unit's own network reference angle — the
                # vector-group rotation between the unit and the reference, which
                # `replicate` already resolved per (load, harmonic) with the right
                # sequence sign. The supplied range then contributes only the
                # dispersion, and E[φ] equals the angle the deterministic model
                # injects, so the dHHC ≡ sHHC point-mass equivalence holds unit by
                # unit rather than only on average.
                if centre_mode == "ref_angle"
                    centre = deg2rad(get(ntw["load"][d_str], "ref_angle", 0.0))
                    v_min, v_max = recentre_range(v_min, v_max, centre, Eω)
                end

                # Local angle reference: the support is given relative to the
                # unit's own fundamental bus voltage and is rotated into the
                # network frame by load["angle_offset"] = h·θ_i (see
                # `update_hdata_with_angle_offsets!`); zero for the deprecated
                # "global" reference.
                offset = deg2rad(get(ntw["load"][d_str], "angle_offset", 0.0))
                v_min, v_max = v_min + offset, v_max + offset

                cc, sc = if haskey(unit_type, d_int)
                    t = unit_type[d_int]; s_min, s_max = type_shifts[nw_int][t]
                    affine2_trig_pce(pb, germ_dim[d_int], v_min, v_max, type_dim[t], s_min, s_max)
                else
                    affine_trig_pce(pb, germ_dim[d_int], v_min, v_max)
                end
                cos_phi[nw_int][d_int] = cc
                sin_phi[nw_int][d_int] = sc
                angle_ranges[nw_int][d_int] = (Float64(v_min), Float64(v_max))
            else
                # normal load: zero harmonic injection
                cos_phi[nw_int][d_int] = zerovec
                sin_phi[nw_int][d_int] = zerovec
            end
        end
    end

    # ── 6. Store in hdata ──────────────────────────────────────────────────────
    hdata["pce"] = Dict(
        "basis"          => pb,
        "M"              => M,
        "lambda"         => lambda,
        "psi_norms"      => psi_norms,
        "K"              => K,
        "K_phys"         => K_phys,
        "cos_phi"        => cos_phi,
        "sin_phi"        => sin_phi,
        "harmonic_units" => harmonic_units,
        "germ_dim"       => germ_dim,
        "type_dim"       => type_dim,
        "unit_type"      => unit_type,
        "type_shifts"    => type_shifts,
        "bound"          => bound,
        "epsilon"        => ε,
        "alpha_beta"     => (Float64(α̂), Float64(β̂)),
        "angle_ranges"   => angle_ranges,
        "lambdas"        => emp_lambdas,   # per-constraint λ* (or nothing)
    )

    return hdata
end

"""
    _warn_floating_harmonic_components(ntw, h)

Warn if harmonic network `ntw` (order `h`) contains a connected component with no
path to ground.  Now that no harmonic bus voltage is forced, the only harmonic
ground paths are the generator subtransient admittance (gsc/bsc), the line/branch
charging (g_fr/b_fr, g_to/b_to), bus shunts (gs/bs), and ideal ("inf") sources.
A component lacking all of these makes the nodal admittance Y_h singular.
"""
function _warn_floating_harmonic_components(ntw::Dict, h::Int)
    adj = Dict{Int, Vector{Int}}()
    for (_, b) in ntw["bus"]
        get!(adj, b["bus_i"], Int[])
    end
    addedge!(a, b) = (push!(get!(adj, a, Int[]), b); push!(get!(adj, b, Int[]), a))
    for (_, br) in ntw["branch"]
        get(br, "br_status", 1) == 0 && continue
        addedge!(br["f_bus"], br["t_bus"])
    end
    if haskey(ntw, "xfmr")
        for (_, x) in ntw["xfmr"]
            addedge!(x["f_bus"], x["t_bus"])
        end
    end

    grounded = Set{Int}()
    for (_, b) in ntw["bus"]
        get(b, "is_clean", false) && push!(grounded, b["bus_i"])   # ideal harmonic ground
    end
    for (_, br) in ntw["branch"]
        (get(br, "b_fr", 0.0) != 0.0 || get(br, "g_fr", 0.0) != 0.0) && push!(grounded, br["f_bus"])
        (get(br, "b_to", 0.0) != 0.0 || get(br, "g_to", 0.0) != 0.0) && push!(grounded, br["t_bus"])
    end
    for (_, s) in get(ntw, "shunt", Dict{String,Any}())
        (get(s, "gs", 0.0) != 0.0 || get(s, "bs", 0.0) != 0.0) && push!(grounded, s["shunt_bus"])
    end
    for (_, g) in ntw["gen"]
        if !iszero(get(g, "inf", 0))
            push!(grounded, g["gen_bus"])                  # ideal source grounds the bus
        elseif get(g, "gsc", 0.0) != 0.0 || get(g, "bsc", 0.0) != 0.0
            push!(grounded, g["gen_bus"])                  # subtransient admittance
        end
    end

    seen = Set{Int}()
    for s in keys(adj)
        s in seen && continue
        comp = Int[]; stack = [s]; push!(seen, s)
        while !isempty(stack)
            u = pop!(stack); push!(comp, u)
            for v in adj[u]
                v in seen || (push!(seen, v); push!(stack, v))
            end
        end
        if isempty(intersect(comp, grounded))
            @warn "harmonic $h: connected component $(sort(comp)) has no path to " *
                  "ground (no generator subtransient impedance, line charging, or " *
                  "bus shunt); the harmonic admittance Y_$h is singular."
        end
    end
end

"""
    HarmonicPowerModels.replicate
"""
function _HPM.replicate(data::Dict{String, Any};
                            bus_id::Int=1,
                            H::Array{Int}=Int[], 
                            xfmr_magn::Dict{String,Any}=Dict{String,Any}())
    ### Add entries to the fundamental data ####################################
    # set the branch current rating
    if haskey(data,"branch") 
        for branch in values(data["branch"])
            f_bus = data["bus"][string(branch["f_bus"])]
            t_bus = data["bus"][string(branch["t_bus"])]
            
            vmmin = min(f_bus["vmin"], t_bus["vmin"])
            
            # PowerModels per-unit uses S = V·I (three-phase folded into the
            # base), so the current limit is rate_a / V — no √3 factor. vmmin
            # (worst-case lowest voltage → highest current) keeps it conservative.
            branch["c_rating"] = branch["rate_a"] / vmmin
    end end

    # add xfmr current rating
    if haskey(data,"xfmr") 
        for xfmr in values(data["xfmr"])
            f_bus = data["bus"][string(xfmr["f_bus"])]
            t_bus = data["bus"][string(xfmr["t_bus"])]
            
            vmmin = min(f_bus["vmin"], t_bus["vmin"])
            
            xfmr["c_rating"] = xfmr["rateA"] / data["baseMVA"] / vmmin        # @Hakan: klopt dit
    end end

    # add the thd limits based on standard, if available
    for bus in values(data["bus"])
        if haskey(bus, "standard") && bus["standard"] in keys(thd_limits)
            bus["thdmax"] = thd_limits[bus["standard"]]
    end end

    ### create multi-network data structure, only keep ntws in H ###############
    hdata = _PMs.replicate(data, last(H))
    for nh in 1:last(H) if nh ∉ H delete!(hdata["nw"],"$nh") end end

    # add bus_id
    hdata["bus_id"] = bus_id

    # rename the case
    hdata["name"] = data["name"]
    haskey(data, "principle") ? hdata["principle"] = data["principle"] : ~ ;
    haskey(data, "angle_reference") ? hdata["angle_reference"] = data["angle_reference"] : ~ ;

    ### add entries to the harmonic data #######################################
    # add the xfmr magnetizing current
    if !isempty(xfmr_magn)
        sample_magnetizing_current(hdata, xfmr_magn)
    end

    # re-evaluate the data for each harmonic 
    for (nw,ntw) in hdata["nw"]
        # translate the ntw-id to an Int for the harmonic number 
        nh = parse(Int,nw)

        # re-evaluate load power
        for load in values(ntw["load"])
            bus = ntw["bus"]["$(load["source_id"][2])"]

            mult = haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 1.0 ;

            # Harmonic current rating of the unit, i.e. the cap on the harmonic
            # current it may inject (bound on crd/cid/cmd, and on the hosting
            # capacity I_h of the stochastic models).  Derived from the unit's
            # RATED FUNDAMENTAL current |S_d| / V_min at its bus, so that a unit
            # is never allowed to inject more harmonic current than the current
            # it is built for.  Computed BEFORE the harmonic power scaling below,
            # which is why it sits here.
            #
            # A "c_rating" already present in the input data wins: that is where
            # a converter nameplate, or a standards-based emission allowance such
            # as the IEEE-519 current-distortion limit (see
            # `ihd_current_limit_ieee_519`, which needs the site's I_sc/I_L ratio
            # and is therefore not derivable from the network data alone), is
            # supplied per load.
            if !haskey(load, "c_rating")
                vmin_d = get(bus, "vmin", 1.0)
                vmin_d > 0.0 || (vmin_d = 1.0)
                load["c_rating"] = sqrt(get(load, "pd", 0.0)^2 +
                                        get(load, "qd", 0.0)^2) / vmin_d
            end

            haskey(load, "pd") ? load["pd"] *= mult : ~ ;
            haskey(load, "qd") ? load["qd"] *= mult : ~ ;
            load["multiplier"] = mult

            if haskey(bus, "ref_angle")
                if is_pos_sequence(nh)
                    load["ref_angle"] = bus["ref_angle"]
                elseif is_neg_sequence(nh)
                    load["ref_angle"] = -bus["ref_angle"]
                elseif is_zero_sequence(nh)
                    load["ref_angle"] = 0.0
            end end

            if haskey(bus, "angle_range")
                load["harmonic_angle_range"] = bus["angle_range"]
            end
        end

        # re-evaluate gen 
        for gen in values(ntw["gen"])
            if haskey(gen, "isfilter") && gen["isfilter"] == 1
                # do nothing, is handled by constraints constraint_active_filter
                gen["pmin"] = -abs(gen["pmax"])
                gen["qmin"] = -abs(gen["qmax"])
            else #is true generator
                if nw != "1" #cost of harmonics set to 0 
                    if haskey(gen, "cost")
                        gen["cost"] *= 0 
                    end
                    #harmonics can be injected/absorbed to match load 
                    gen["pmin"] = -abs(gen["pmax"])
                    gen["qmin"] = -abs(gen["qmax"])
            end end
            gen_bus = gen["gen_bus"]
            gen["c_rating"] = sqrt(max(abs(gen["pmin"]), abs(gen["pmax"]))^2 + max(abs(gen["qmin"]), abs(gen["qmax"]))^2) / (sqrt(3) * ntw["bus"]["$gen_bus"]["vmin"])
        end

        # re-evaluate the bus data 
        for bus in values(ntw["bus"])
            # use fundamental as limit for rms
            bus["vminrms"] = bus["vmin"]
            bus["vmaxrms"] = bus["vmax"]
            
            # true harmonics don't have minimum voltage
            if nw != "1" 
                bus["vmin"] = 0.0 
            end

            # add the ihd limits based on standard, if available
            if haskey(bus, "standard") 
                std = bus["standard"]
                if std in keys(ihd_limits)
                    if nh <= length(ihd_limits[std])
                        bus["ihdmax"] = ihd_limits[std][nh]
                    else
                        @warn "harmonic $nh not included in $std"
        end end end end

        # re-evaluate the branch data
        for branch in values(ntw["branch"])
            haskey(branch, "br_r") ? branch["br_r"] *= sqrt(nh) : ~ ;
            haskey(branch, "br_x") ? branch["br_x"] *= nh : ~ ;
            haskey(branch, "b_fr") ? branch["b_fr"] *= nh : ~ ;
            haskey(branch, "b_to") ? branch["b_to"] *= nh : ~ ;
        end

        # re-evaluate the bus shunt data: a shunt is a fixed physical element, so
        # its ADMITTANCE is frequency-dependent even though its rating is not.
        # In the MATPOWER sign convention bs > 0 is capacitive and bs < 0 is
        # inductive, so at harmonic h
        #   capacitor: B(h) = ω(h)·C       = bs·h
        #   reactor:   B(h) = −1/(ω(h)·L)  = bs/h
        # A conductive shunt (gs) is frequency-independent. Leaving bs unscaled
        # understates a capacitor bank by a factor h and overstates a reactor by
        # the same factor, which is exactly where harmonic resonances live.
        for shunt in values(get(ntw, "shunt", Dict{String,Any}()))
            haskey(shunt, "bs") || continue
            shunt["bs"] = shunt["bs"] >= 0.0 ? shunt["bs"] * nh : shunt["bs"] / nh
        end

        # re-evaluate the gen data: the per-generator short-circuit impedance
        # fields (rsc, xsc) are the fundamental per-unit subtransient impedance
        # z_sc = rsc + j*xsc (see mpc.gen_sc).  At harmonic h the generator
        # presents the harmonic-scaled impedance
        #   z_sc(h) = rsc*sqrt(h) + j*xsc*h ,
        # so the shunt admittance to ground used by the harmonic constraints is
        #   Y_g(h) = gsc + j*bsc = 1/z_sc(h).
        # rsc = xsc = 0 (clean machine, no gen_sc row) ⇒ no ground path.
        for gen in values(ntw["gen"])
            rsc = get(gen, "rsc", 0.0)
            xsc = get(gen, "xsc", 0.0)
            if !iszero(rsc) || !iszero(xsc)
                Yh = 1 / (rsc * sqrt(nh) + im * xsc * nh)        # Y_g(h)
                gen["gsc"] = real(Yh)
                gen["bsc"] = imag(Yh)
            end
        end

        # re-evaluate the transformer data
        if haskey(ntw, "xfmr") 
            for xfmr in values(ntw["xfmr"])
                haskey(xfmr, "xsc") ? xfmr["xsc"] *= nh : ~ ;
                haskey(xfmr, "r1")  ? xfmr["r1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "r2")  ? xfmr["r2"] *= sqrt(nh) : ~ ;

                haskey(xfmr, "xe1") ? xfmr["xe1"] *= nh : ~ ;
                haskey(xfmr, "xe2") ? xfmr["xe2"] *= nh : ~ ;
                haskey(xfmr, "re1") ? xfmr["re1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "re2") ? xfmr["re2"] *= sqrt(nh) : ~ ;

                xfmr["cnf1"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][1]) : 'Y' ;
                xfmr["cnf2"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][2]) : 'Y' ;

                shift = haskey(xfmr, "vg") ? parse(Int, xfmr["vg"][3:end]) : 0 ;
                if is_pos_sequence(nh)
                    xfmr["tr"] = cosd(-30.0 * shift)
                    xfmr["ti"] = sind(-30.0 * shift)
                elseif is_neg_sequence(nh)
                    xfmr["tr"] = cosd(30.0 * shift)
                    xfmr["ti"] = sind(30.0 * shift)
                elseif is_zero_sequence(nh)
                    xfmr["tr"] = 1.0
                    xfmr["ti"] = 0.0
        end end end
    end

    # warn about floating (ungrounded) harmonic components ⇒ singular Y_h
    for (nw, ntw) in hdata["nw"]
        nw == "1" && continue
        _warn_floating_harmonic_components(ntw, parse(Int, nw))
    end

    return hdata
end
