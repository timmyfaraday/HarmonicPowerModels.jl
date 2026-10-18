################################################################################
# HarmonicPowerModels.jl — binding-constraint diagnostic, IEEE30 scenarios      #
#                                                                               #
# For each (scenario, N, angle reference):                                      #
#   1. solve dHHC_SOC at the expected angles E[φ] and list, per harmonic, the   #
#      constraints with the highest utilisation (IHD, THD, RMS voltage, branch  #
#      RMS current, generator RMS current, unit rating);                        #
#   2. at the binding voltage constraint (bus i*, harmonic h), take the units'  #
#      contributions c_j = Z_h[i*, bus_j] · e^{jφ_j} per unit of current and     #
#      report:                                                                   #
#        spread   smallest arc containing the contribution angles at E[φ]       #
#        α(E)     summation exponent at E[φ]: |Σc_j| = (Σ|c_j|^α)^(1/α)          #
#        α(95)    the same with |Σc_j| replaced by its 95th percentile over the  #
#                 angle distribution (Monte Carlo) — the IEC 61000-3-6 reading   #
#        rob/q95  worst case of |Σc_j| over the angle support (support function) #
#                 divided by its 95th percentile: the current-hosting gain a 95% #
#                 chance bound buys over a robust bound for this constraint      #
#        E/q95    |Σc_j| at E[φ] divided by its 95th percentile                  #
#                                                                               #
# Output: output/ieee30_binding.txt                                             #
# Author: Tom Van Acker                                                          #
################################################################################

using HarmonicPowerModels, PowerModels, Ipopt, MosekTools, JuMP
using Printf, Random, Distributions, LinearAlgebra
const HPM = HarmonicPowerModels
PowerModels.silence()
include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "ieee30_scenarios.jl"))

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
harm  = [5, 7, 11, 13]
H     = [1; harm]
α̂, β̂  = PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]
Eω    = α̂ / (α̂ + β̂)
rng   = Range = Dict(h => (PCE_PARAMS["harmonics"][h].v_min, PCE_PARAMS["harmonics"][h].v_max) for h in harm)
NMC   = 50_000

"Smallest arc [deg] containing all angles [rad]."
function arc(angs)
    a = sort(mod.(rad2deg.(angs), 360.0))
    360.0 - maximum(vcat(diff(a), a[1] + 360.0 - a[end]))
end

"Summation exponent α ≥ 1 with (Σ m_j^α)^(1/α) = target (bisection on [1, 20])."
function sum_exponent(m, target)
    f(α) = sum(m .^ α)^(1 / α) - target
    target >= sum(m) && return 1.0
    target <= maximum(m) && return Inf
    lo, hi = 1.0, 20.0
    for _ in 1:60
        mid = (lo + hi) / 2
        f(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"Worst case of |Σ z_j e^{jφ_j}| over φ_j ∈ [lo_j, hi_j] (support function, θ scan)."
function robust_abs(z, lo, hi)
    best = 0.0
    for θ in range(0, 2π, length = 1441)
        s = 0.0
        for j in eachindex(z)
            # max over φ ∈ [lo, hi] of |z| cos(arg z + φ − θ)
            a0 = angle(z[j]) - θ
            φstar = mod(-a0 - lo[j], 2π) + lo[j]          # φ with cos(a0+φ)=1, mapped ≥ lo
            v = φstar <= hi[j] ? abs(z[j]) :
                abs(z[j]) * max(cos(a0 + lo[j]), cos(a0 + hi[j]))
            s += v
        end
        best = max(best, s)
    end
    return best
end

function diagnose(io, scen, N, aref)
    data = ieee30_scenario(scen, ipopt)
    units = sort(parse.(Int, collect(keys(data["load"]))))[1:N]
    hd = HPM.replicate(deepcopy(data), H = H)
    hd["angle_reference"] = aref
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    # deterministic angles at E[φ] (relative to the reference); offsets added by the model
    for h in harm, l in values(hd["nw"]["$h"]["load"])
        lo, hi = rng[h]; l["ref_angle"] = rad2deg(lo + (hi - lo) * Eω)
    end
    sol = HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt)
    s, f = sol["solution"]["nw"], hd["nw"]["1"]
    I = Dict(h => s["$h"]["load"]["$(units[1])"]["cmd"] for h in harm)

    # ── utilisation of every constraint at the dHHC solution ──────────────────
    util = Tuple{Float64,String}[]
    for (b, bus) in f["bus"]
        vm = bus["vm"]; V2 = Dict(h => s["$h"]["bus"][b]["vr"]^2 + s["$h"]["bus"][b]["vi"]^2 for h in harm)
        for h in harm
            push!(util, (sqrt(V2[h]) / (hd["nw"]["$h"]["bus"][b]["ihdmax"] * vm), "IHD  bus $b h=$h"))
        end
        push!(util, (sqrt(sum(values(V2))) / (bus["thdmax"] * vm), "THD  bus $b"))
        push!(util, (sqrt(sum(values(V2)) / (bus["vmaxrms"]^2 - vm^2)), "RMSV bus $b"))
    end
    for (k, br) in f["branch"]
        for (side, cm) in (("fr", br["cm_fr"]), ("to", br["cm_to"]))
            i2 = sum(s["$h"]["branch"][k]["cr_$side"]^2 + s["$h"]["branch"][k]["ci_$side"]^2 for h in harm)
            push!(util, (sqrt(i2 / (br["c_rating"]^2 - cm^2)), "RMSI branch $(br["f_bus"])-$(br["t_bus"]) $side"))
        end
    end
    # generator RMS current is not a constraint of dHHC_SOC; reported for information
    gen_notes = String[]
    for (g, gen) in f["gen"]
        i2 = sum(s["$h"]["gen"][g]["crg"]^2 + s["$h"]["gen"][g]["cig"]^2 for h in harm)
        head = gen["c_rating"]^2 - gen["cm"]^2
        head > 0 ? push!(gen_notes, @sprintf("gen@%d %.3f", gen["gen_bus"], sqrt(i2 / head))) :
                   push!(gen_notes, @sprintf("gen@%d over rating at the fundamental (cm=%.3f > %.3f)",
                                             gen["gen_bus"], gen["cm"], gen["c_rating"]))
    end
    for h in harm
        cap = minimum(hd["nw"]["$h"]["load"]["$d"]["c_rating"] for d in units)
        push!(util, (I[h] / cap, "UNIT rating h=$h"))
    end
    sort!(util, rev = true)

    println(io, repeat('─', 92))
    @printf(io, "%s  N=%d  %s   dHHC@E[φ] Σ_h I_h = %.4e   I_h = %s\n", scen, N, aref, sum(values(I)),
            join([@sprintf("%d:%.3e", h, I[h]) for h in harm], " "))
    println(io, "  top constraints (utilisation):")
    for (u, nm) in util[1:8]
        @printf(io, "    %6.4f  %s\n", u, nm)
    end
    println(io, "  generator RMS current (not enforced in dHHC_SOC): ", join(gen_notes, ", "))

    # ── recombination geometry at the binding voltage constraints ─────────────
    Z, bus_pos = HPM.reduced_harmonic_impedance(hd)
    ubus = Dict(d => hd["nw"]["1"]["load"]["$d"]["load_bus"] for d in units)
    targets = unique([(parse(Int, m[1]), h) for (u, nm) in util[1:8] if u > 0.99
                      for m in (match(r"bus (\d+)", nm),) if m !== nothing
                      for h in (occursin("h=", nm) ? [parse(Int, match(r"h=(\d+)", nm)[1])] : harm)])
    isempty(targets) && println(io, "  no binding voltage constraint")
    beta = Beta(α̂, β̂); Random.seed!(1)
    for (i, h) in targets
        z  = [Z[h][bus_pos[i], bus_pos[ubus[d]]] for d in units]
        off = [deg2rad(get(hd["nw"]["$h"]["load"]["$d"], "angle_offset", 0.0)) for d in units]
        lo, hi = rng[h]
        lo_j = lo .+ off; hi_j = hi .+ off
        φE = lo_j .+ (hi - lo) * Eω
        cE = z .* cis.(φE)
        mags = abs.(z)
        samples = [abs(sum(z .* cis.(lo_j .+ (hi - lo) .* rand(beta, N)))) for _ in 1:NMC]
        q95 = quantile(samples, 0.95)
        rob = robust_abs(z, lo_j, hi_j)
        @printf(io, "  bus %-3d h=%-3d spread %6.1f°  α(E)=%5.3f  α(95)=%5.3f  rob/q95=%6.4f  E/q95=%6.4f  (Σ|c|/|ΣcE|=%.3f)\n",
                i, h, arc(angle.(cE)), sum_exponent(mags, abs(sum(cE))), sum_exponent(mags, q95),
                rob / q95, abs(sum(cE)) / q95, sum(mags) / abs(sum(cE)))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    open(joinpath(@__DIR__, "output", "ieee30_binding.txt"), "w") do io
        println(io, "IEEE30 binding-constraint diagnostic (dHHC_SOC at E[φ], Mosek; geometry from Z_h, $NMC MC samples)")
        for scen in ["base", "banks", "cable132c"], N in (8, 21), aref in ("global", "local")
            diagnose(io, scen, N, aref); flush(io)
        end
    end
    println("✓ ieee30_binding.txt")
end
