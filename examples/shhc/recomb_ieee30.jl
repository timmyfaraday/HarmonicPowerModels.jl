################################################################################
# HarmonicPowerModels.jl — feasible recombination vs the summation law, IEEE30  #
#                                                                               #
# The individual-harmonic voltage limits bind in every scenario (see            #
# ieee30_binding.txt), so per harmonic the absolute-equality hosting current is #
# the tightest bus, I_h = min_i L_ih / A_ih, with L_ih = ihdmax·|V_i,1| and A_ih #
# the aggregate voltage per unit of current from the units' contributions        #
# c_ij = Z_h[i, bus_j] e^{jφ_j}. Six aggregation rules:                          #
#   arith  Σ_j |c_ij|                         no recombination (α = 1)           #
#   rob    max over the angle support of |Σ_j c_ij|                              #
#   det    |Σ_j c_ij| at E[φ]                 the deterministic angle            #
#   q95    95th percentile of |Σ_j c_ij|      exact per-constraint chance bound  #
#   iec    (Σ_j |c_ij|^α)^(1/α), α = 1.4 for h < 10 and 2.0 for h > 10           #
#          (IEC 61000-3-6 indicative summation exponents)                        #
#   rss    (Σ_j |c_ij|^2)^(1/2)               random, uniformly distributed      #
#          angles (α = 2)                                                        #
# For det and q95 the equivalent summation exponent α_eq is also reported: the  #
# α for which the summation law gives the same network hosting current.          #
#                                                                               #
# Output: output/ieee30_recombination.txt                                       #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "diag_binding_ieee30.jl"))

iec_alpha(h) = h < 10 ? 1.4 : 2.0
lawsum(m, α) = sum(m .^ α)^(1 / α)

"""
    recombination(scen, N, aref; ranges = nothing, types = nothing)

`ranges[h][d] = (lo, hi)` [rad] sets a per-unit angle support (e.g., from
`mix_ranges` in converter_types.jl); by default every unit takes the §II range.
`eval_I[h]` evaluates the violation probability of the binding individual-harmonic
voltage limit at the hosting current `eval_I[h]` (returned as `:viol`).
`types` adds shared type germs in the `pce_params["type_germs"]` format
(name => Dict("units" => …, "shift" => Dict(h => (s_lo, s_hi)))): every unit of a
type gets the extra shift s_lo + (s_hi − s_lo)·ω_t with one ω_t per type. The
robust bound then maximises over ω_t on a grid as well.
"""
function recombination(scen, N, aref; ranges = nothing, types = nothing, data = nothing, eval_I = nothing, prob = 0.95)
    data = data === nothing ? ieee30_scenario(scen, ipopt) : data
    units = sort(parse.(Int, collect(keys(data["load"]))))[1:N]
    hd = HPM.replicate(deepcopy(data), H = H)
    hd["angle_reference"] = aref
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    HPM.update_hdata_with_fundamental_hpf_results!(hd, dHHC_NLP, ipopt)
    Z, bus_pos = HPM.reduced_harmonic_impedance(hd)
    f = hd["nw"]["1"]
    ubus = Dict(d => f["load"]["$d"]["load_bus"] for d in units)
    beta = Beta(α̂, β̂); Random.seed!(1)
    res = Dict{Int,Dict{Symbol,Float64}}()
    for h in harm
        off = [deg2rad(get(hd["nw"]["$h"]["load"]["$d"], "angle_offset", 0.0)) for d in units]
        lo_u = [ranges === nothing ? rng[h][1] : ranges[h][d][1] for d in units]
        hi_u = [ranges === nothing ? rng[h][2] : ranges[h][d][2] for d in units]
        lo_j = lo_u .+ off; hi_j = hi_u .+ off; wd = hi_u .- lo_u
        # shared type germs: per type, the units' indicator and the shift range
        tg = types === nothing ? [] :
             [([d in tt["units"] for d in units], tt["shift"][h]) for tt in values(types)]
        tshift(ζ) = isempty(tg) ? zeros(N) :
                    sum(mask .* (sh[1] + (sh[2] - sh[1]) * ζ[i]) for (i, (mask, sh)) in enumerate(tg))
        φE = lo_j .+ wd .* Eω .+ tshift(fill(Eω, length(tg)))
        W  = [lo_j .+ wd .* rand(beta, N) .+ tshift(rand(beta, length(tg))) for _ in 1:20_000]
        # robust bound: exact over a grid of the type germ for at most one type; with
        # several types the type shifts are bounded per unit (a conservative bound)
        function robust_types(z)
            length(tg) <= 1 && return maximum(robust_abs(z, lo_j .+ tshift(ζ), hi_j .+ tshift(ζ))
                                              for ζ in (isempty(tg) ? [Float64[]] : [[x] for x in range(0, 1, length = 41)]))
            slo = sum(mask .* min(sh...) for (mask, sh) in tg); shi = sum(mask .* max(sh...) for (mask, sh) in tg)
            return robust_abs(z, lo_j .+ slo, hi_j .+ shi)
        end
        cap = minimum(hd["nw"]["$h"]["load"]["$d"]["c_rating"] for d in units)
        buses = [(L = hd["nw"]["$h"]["bus"][b]["ihdmax"] * bus["vm"],
                  z = [Z[h][bus_pos[parse(Int, b)], bus_pos[ubus[d]]] for d in units])
                 for (b, bus) in f["bus"] if haskey(bus_pos, parse(Int, b))]
        host(A) = min(cap, minimum(bb.L / A(bb.z) for bb in buses))
        r = Dict(
            :arith => host(z -> sum(abs.(z))),
            :rob   => host(z -> robust_types(z)),
            :det   => host(z -> abs(sum(z .* cis.(φE)))),
            :q95   => host(z -> quantile([abs(sum(z .* cis.(w))) for w in W], prob)),   # (1-ε) quantile
            :iec   => host(z -> lawsum(abs.(z), iec_alpha(h))),
            :rss   => host(z -> lawsum(abs.(z), 2.0)),
            :capped => cap)
        # violation probability of the individual harmonic voltage limits at a given
        # hosting current (e.g., a second-order cone solution): worst bus, per sample set
        if eval_I !== nothing
            Ih = eval_I[h]
            r[:viol] = maximum(count(w -> Ih * abs(sum(bb.z .* cis.(w))) > bb.L, W) / length(W) for bb in buses)
        end
        # equivalent summation exponent: α with host(law α) = target (bisection)
        for key in (:det, :q95)
            target = r[key]
            if target >= cap * (1 - 1e-9)
                r[Symbol(key, :_alpha)] = NaN            # rating-capped, exponent undefined
            elseif target <= r[:arith]
                r[Symbol(key, :_alpha)] = 1.0
            else
                a, b = 1.0, 8.0
                host(z -> lawsum(abs.(z), b)) < target && (r[Symbol(key, :_alpha)] = Inf; continue)
                for _ in 1:50
                    m = (a + b) / 2
                    host(z -> lawsum(abs.(z), m)) < target ? (a = m) : (b = m)
                end
                r[Symbol(key, :_alpha)] = (a + b) / 2
            end
        end
        res[h] = r
    end
    return res
end

fmt(x) = isnan(x) ? "  cap" : isinf(x) ? "  >8 " : @sprintf("%5.2f", x)
pct(a, b) = @sprintf("%+6.1f%%", 100(a / b - 1))

if abspath(PROGRAM_FILE) == @__FILE__
    open(joinpath(@__DIR__, "output", "ieee30_recombination.txt"), "w") do io
        println(io, "IEEE30 — feasible recombination vs the IEC 61000-3-6 summation law (binding: IHD).")
        println(io, "α_eq: summation exponent reproducing the network hosting current (IEC: 1.4 for h=5,7; 2.0 for h=11,13).")
        println(io, "Hosting Σ_h I_h relative to the deterministic angle E[φ]; 'cap' = unit rating binds.\n")
        @printf(io, "%-10s %2s %-6s | %-23s | %-23s | %8s %8s %8s %8s %8s\n", "scenario", "N", "ref",
                "α_eq det  h=5/7/11/13", "α_eq q95  h=5/7/11/13", "arith", "robust", "q95", "IEC", "rss")
        for scen in ["base", "banks", "cable33", "cable132c", "cableallc"], N in (8, 21), aref in ("global", "local")
            r = recombination(scen, N, aref)
            tot(k) = sum(r[h][k] for h in harm)
            @printf(io, "%-10s %2d %-6s | %s | %s | %8s %8s %8s %8s %8s\n", scen, N, aref,
                    join([fmt(r[h][:det_alpha]) for h in harm], " "),
                    join([fmt(r[h][:q95_alpha]) for h in harm], " "),
                    pct(tot(:arith), tot(:det)), pct(tot(:rob), tot(:det)), pct(tot(:q95), tot(:det)),
                    pct(tot(:iec), tot(:det)), pct(tot(:rss), tot(:det)))
            flush(io)
        end
    end
    println("✓ ieee30_recombination.txt")
end
