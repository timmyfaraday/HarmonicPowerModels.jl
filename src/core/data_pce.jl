################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.1 - added PCE data layer                                                 #
################################################################################

# ── nw indexing helpers ──────────────────────────────────────────────────────
# Flat (h_idx, k) → nw_id encoding.
# h_idx is the harmonic number (e.g. 1, 3, 5); k is the 0-based PCE mode.

nw_id(h_idx::Int, k::Int, P_size::Int)     = (h_idx - 1) * P_size + k + 1
get_pce_mode(nw::Int, P_size::Int)         = (nw - 1) % P_size
get_harmonic_idx(nw::Int, P_size::Int)     = div(nw - 1, P_size) + 1

# ── Input validation ─────────────────────────────────────────────────────────

function check_shhc_data(data::Dict)
    haskey(data, "sdata")               || error("Missing data[\"sdata\"]")
    haskey(data["sdata"], "alpha")      || error("Missing sdata[\"alpha\"]")
    haskey(data["sdata"], "beta")       || error("Missing sdata[\"beta\"]")
    haskey(data["sdata"], "epsilon")    || error("Missing sdata[\"epsilon\"]")
    ε = data["sdata"]["epsilon"]
    (0 < ε < 1)                         || error("epsilon must be in (0,1), got $ε")

    # Validate per-load stochastic angle data.
    for (u_id, u_data) in data["load"]
        haskey(u_data, "sdata") || error("Load $u_id missing \"sdata\"")
        for h in keys(u_data["sdata"])
            sd = u_data["sdata"][h]
            haskey(sd, "angle_min_deg") ||
                error("Load $u_id harmonic $h missing angle_min_deg")
            haskey(sd, "angle_max_deg") ||
                error("Load $u_id harmonic $h missing angle_max_deg")
        end
    end
end

# ── Lambda computation ───────────────────────────────────────────────────────

function compute_lambda!(data::Dict, ::PCEData)
    ε = data["sdata"]["epsilon"]
    λ = sqrt((1.0 - ε) / ε)          # Cantelli one-sided bound
    data["sdata"]["lambda"] = λ

    # Store λ at bus/harmonic level (if bus has a "harmonics" key).
    for (_, n_data) in data["bus"]
        if haskey(n_data, "harmonics")
            for h in n_data["harmonics"]
                h_str = string(h)
                if !haskey(n_data, "harmonic_sdata")
                    n_data["harmonic_sdata"] = Dict()
                end
                n_data["harmonic_sdata"][h_str] = Dict(
                    "lambda_ihd" => λ,
                    "lambda_thd" => λ
                )
            end
        end
    end

    # Store λ at branch/harmonic level (if branch has a "harmonics" key).
    for (_, b_data) in data["branch"]
        if haskey(b_data, "harmonics")
            for h in b_data["harmonics"]
                if !haskey(b_data, "harmonic_sdata")
                    b_data["harmonic_sdata"] = Dict()
                end
                b_data["harmonic_sdata"][string(h)] = Dict(
                    "lambda_current" => λ
                )
            end
        end
    end
end

# ── Multi-network PCE expansion ──────────────────────────────────────────────

"""
    build_mn_pce_data(data, pce; bus_id, H, xfmr_magn) -> Dict

Build the (harmonic × PCE-mode) expanded multi-network dict for the sHHC problem.

The returned dict has `n_harm * pce.P_size` sub-networks, keyed by the integer
`nw_id(h_idx, k, pce.P_size)` where `h_idx` is the harmonic number and `k` is
the 0-based PCE mode index.  Mode `k == 0` is the mean (deterministic) mode.

`data` must be a flat (single-network) PowerModels data dict augmented with:
  - `data["sdata"]` containing `"alpha"`, `"beta"`, `"epsilon"`
  - `data["load"][id]["sdata"][h]` containing `"angle_min_deg"`, `"angle_max_deg"`

`H`, `bus_id`, and `xfmr_magn` are forwarded verbatim to `replicate`.
"""
function build_mn_pce_data(data::Dict, pce::PCEData;
                            bus_id::Int = 1,
                            H::Array{Int} = Int[],
                            xfmr_magn::Dict{String,Any} = Dict{String,Any}())::Dict
    check_shhc_data(data)
    compute_lambda!(data, pce)

    # Step 1: Build per-harmonic multi-network using the existing replicate function.
    mn_harm = replicate(data; bus_id = bus_id, H = H, xfmr_magn = xfmr_magn)

    # Step 2: Determine the number of harmonics from the existing mn dict.
    n_harm  = length(mn_harm["nw"])

    # Step 3: Build the (h_idx, k) expanded multi-network.
    mn = Dict{String, Any}(
        "multinetwork" => true,
        "pce"          => pce,
        "P_size"       => pce.P_size,
        "n_harmonics"  => n_harm,
        "nw"           => Dict{String, Any}()
    )

    # Copy top-level non-nw keys from mn_harm.
    for (k, v) in mn_harm
        k == "nw"          && continue
        k == "multinetwork" && continue
        mn[k] = v
    end

    # Expand each harmonic sub-network into P_size PCE-mode copies.
    for (h_str, h_nw) in mn_harm["nw"]
        h_idx = parse(Int, h_str)
        for k in 0:pce.deg
            new_id             = string(nw_id(h_idx, k, pce.P_size))
            sub_nw             = deepcopy(h_nw)
            sub_nw["harmonic_idx"]  = h_idx
            sub_nw["pce_mode"]      = k
            sub_nw["is_mean_mode"]  = (k == 0)
            mn["nw"][new_id]        = sub_nw
        end
    end

    return mn
end
