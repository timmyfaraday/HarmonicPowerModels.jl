# HarmonicPowerModels.jl                                                       #
# Phase-0 sanity check: validate the LP-reduction PCE pipeline against a       #
# brute-force Monte-Carlo simulation that bypasses PCE entirely.               #
#                                                                               #
# For each network and harmonic h, with unit injected current per load:        #
#   φ_u = v_min + (v_max − v_min)·ω_u,   ω_u ~ Beta(α̂, β̂) iid (independent) #
#                                          or one shared ω (shared germ)        #
#   V_n = Σ_u Z_h[n,u]·e^{jφ_u},   S_n = |V_n|²                              #
#   C_mc[n] = E[S_n] + λ·√Var[S_n]      (Cantelli coefficient)                #
#   I*_mc   = min_n √(V_lim²[n] / C_mc[n])                                    #
#                                                                               #
# These are compared against compute_C_scalar / lp_solution outputs            #
# (i.e. the values snapshotted in the Phase-0 goldens). Agreement within a     #
# few percent validates the δ=2 PCE + degree-4 Taylor pipeline; the residual   #
# is genuine truncation error, not a bug.                                      #
#                                                                               #
# Run: julia --project test/golden/validate_pce_mc.jl                          #
# Author: Tom Van Acker                                                         #

using HarmonicPowerModels, PowerModels
using Distributions, Random, Printf
import Statistics

include(joinpath(@__DIR__, "..", "..", "examples", "shhc", "pce_params.jl"))

const PMs = PowerModels
const HPM = HarmonicPowerModels

const HARM    = [5, 7, 11, 13]
const N_MC    = 200_000
const RNG     = MersenneTwister(20260612)

function validate(tag::String, casefile::String)
    data = PMs.parse_file(joinpath(@__DIR__, "..", "..", "examples", "shhc",
                                   "data", casefile))

    α̂ = PCE_PARAMS["alpha_hat"]
    β̂ = PCE_PARAMS["beta_hat"]
    δ = PCE_PARAMS["delta"]
    ε = PCE_PARAMS["epsilon"]
    λ = sqrt((1 - ε) / ε)

    beta_params = (alpha_shape = α̂, beta_shape = β̂, delta = δ)
    dist        = Beta(α̂, β̂)

    load_entries = sort([(id = parse(Int, k), bus = v["load_bus"])
                         for (k, v) in data["load"]], by = x -> x.id)
    unit_ids     = [l.id for l in load_entries]
    load_bus_ids = [l.bus for l in load_entries]
    N            = length(unit_ids)

    println("\n── $tag  (N=$N loads, $N_MC samples) ──────────────────────────")
    @printf("%-3s  %-12s %-12s %-7s │ %-12s %-12s %-7s\n",
            "h", "I*_ind PCE", "I*_ind MC", "Δ [%]", "I*_sh PCE", "I*_sh MC", "Δ [%]")

    for h in HARM
        hb = PCE_PARAMS["harmonics"][h]
        v_min, v_max = hb.v_min, hb.v_max
        angle_params_h = Dict(uid => (rad2deg(v_min), rad2deg(v_max))
                              for uid in unit_ids)

        # PCE pipeline (same call the golden snapshot uses)
        r = HPM.compare_shared_vs_independent(data, h, ε, beta_params,
                                              angle_params_h)

        # Monte-Carlo
        Z_h     = HPM.harmonic_transfer_impedance(data, h, load_bus_ids)
        active  = HPM._active_bus_ids_lp(data, h)
        limits  = HPM.bus_voltage_limits(data, h)
        n_buses = size(Z_h, 1)

        S_ind = zeros(N_MC, n_buses)
        S_sh  = zeros(N_MC, n_buses)
        for s in 1:N_MC
            ω_ind = rand(RNG, dist, N)
            ω_sh  = rand(RNG, dist)
            e_ind = cis.(v_min .+ (v_max - v_min) .* ω_ind)
            e_sh  = cis(v_min + (v_max - v_min) * ω_sh)
            for n in 1:n_buses
                v_i = sum(Z_h[n, u] * e_ind[u] for u in 1:N)
                v_s = e_sh * sum(Z_h[n, u] for u in 1:N)
                S_ind[s, n] = abs2(v_i)
                S_sh[s, n]  = abs2(v_s)
            end
        end

        function istar_mc(S)
            I, bus = Inf, -1
            for (n, bus_id) in enumerate(active)
                lim = get(limits, bus_id, 0.0)
                lim > 0.0 || continue
                μ, σ² = Statistics.mean(S[:, n]), Statistics.var(S[:, n])
                C = μ + λ * sqrt(σ²)
                cand = sqrt(lim / C)
                cand < I && ((I, bus) = (cand, bus_id))
            end
            return I, bus
        end

        I_ind_mc, b_ind_mc = istar_mc(S_ind)
        I_sh_mc,  b_sh_mc  = istar_mc(S_sh)

        Δi = (r.I_star_independent / I_ind_mc - 1) * 100
        Δs = (r.I_star_shared      / I_sh_mc  - 1) * 100
        flag = max(abs(Δi), abs(Δs)) > 5 ? "  ⚠" : ""
        @printf("%-3d  %-12.6e %-12.6e %+-7.2f │ %-12.6e %-12.6e %+-7.2f%s\n",
                h, r.I_star_independent, I_ind_mc, Δi,
                r.I_star_shared, I_sh_mc, Δs, flag)
        if r.binding_bus_independent != b_ind_mc || r.binding_bus_shared != b_sh_mc
            println("     binding bus mismatch: PCE ($(r.binding_bus_independent), " *
                    "$(r.binding_bus_shared)) vs MC ($b_ind_mc, $b_sh_mc)")
        end
    end
end

validate("industrial", "industrial_network_shhc.m")
validate("ieee30",     "30bus_network_shhc.m")
