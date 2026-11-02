################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): the pessimizing oracle.             #
#                                                                               #
# For hosting variables g, find the parameter points of the box that violate   #
# the chance constraints most: scan the dictionary (compact form, cheap), then #
# refine the best points of the most violated constraints by projected-        #
# gradient ascent on the constraint itself. This is an approximate             #
# pessimizing oracle in the sense of Mutapcic & Boyd (2009, §4.2): it returns  #
# feasible parameter points but does not guarantee the maximum. Its effort is  #
# set by the number of constraints refined, the starts per constraint and the  #
# ascent iterations.                                                           #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

"""
    OracleResult

Result of [`pessimize!`](@ref): `violation[j]` is the worst normalized violation found for
constraint j (in the order of `dict.keys`), `point[j]` the index in the dictionary of the
point attaining it, and `worst` the index of the most violated constraint.
"""
struct OracleResult
    violation::Vector{Float64}
    point::Vector{Int}
    worst::Int
end

"Largest normalized violation found by the oracle."
max_violation(r::OracleResult) = r.violation[r.worst]

# every normalized violation at one latent point, from the full constraint records
function _violations_at(hdata, model, box, ξ, g)
    m = moments_at(hdata, model, latent(box, ξ))
    return [normalized_violation(c, g, m) for c in shhc_constraints(hdata, m)]
end

# Projected-gradient ascent of f over [-1, 1]^d from ξ0: central finite differences (step
# `fd`), a normalised step halved until f improves; stops after `iters` steps or when no
# step improves.
function _ascend(f, ξ0::AbstractVector; iters::Int, fd::Real)
    ξ, fξ = collect(Float64, ξ0), f(ξ0)
    d = length(ξ)
    for _ in 1:iters
        g = zeros(d)
        for i in 1:d
            e = zeros(d); e[i] = fd
            g[i] = (f(clamp.(ξ .+ e, -1, 1)) - f(clamp.(ξ .- e, -1, 1))) / (2fd)
        end
        norm(g) > 0 || break
        α, improved = 0.5, false
        while α > 1e-3
            x  = clamp.(ξ .+ α .* g ./ norm(g), -1, 1)
            fx = f(x)
            if fx > fξ
                ξ, fξ, improved = x, fx, true
                break
            end
            α /= 2
        end
        improved || break
    end
    return ξ, fξ
end

"""
    pessimize!(dict, hdata, model, g; ncons = 5, starts = 3, iters = 10, fd = 1e-3,
               margin = Inf) -> OracleResult

Approximate pessimizing oracle for the hosting variables `g` (h ⇒ g_h = I_h²):

1. evaluate every constraint at every point of `dict` (compact form);
2. for the `ncons` constraints with the largest violation, among those whose dictionary
   violation exceeds `-margin` (i.e., within `margin` of their limit somewhere in the
   dictionary), start a projected-gradient ascent on that constraint's normalized violation
   from its `starts` best points (`iters` steps, finite-difference step `fd` in latent
   coordinates);
3. append every point that improves on the dictionary to `dict` (label `:refined`), so the
   dictionary keeps what the oracle learned, and update every constraint's worst value with
   it.

The effort parameters are `ncons`, `margin`, `starts` and `iters`; the sign-off of the
cutting-set method should run at maximum effort (Mutapcic & Boyd 2009, §5.1).
"""
function pessimize!(dict::ParameterDictionary, hdata::Dict, model::AbstractNetworkModel, g;
                    ncons::Int = 5, starts::Int = 3, iters::Int = 10, fd::Real = 1e-3,
                    margin::Real = Inf)
    gv  = [g[h] for h in dict.harmonics]
    V   = dictionary_violations(dict, g)
    nc  = length(dict.keys)
    best = [maximum(view(V, j, :)) for j in 1:nc]
    arg  = [argmax(view(V, j, :)) for j in 1:nc]

    cache = Dict{Vector{Float64},Vector{Float64}}()
    viol(ξ) = get!(() -> _violations_at(hdata, model, dict.box, ξ, g), cache, collect(Float64, ξ))

    near = [j for j in sortperm(best; rev = true) if best[j] >= -margin]
    for j in near[1:min(ncons, length(near))]
        for s in sortperm(V[j, :]; rev = true)[1:min(starts, size(V, 2))]
            ξ, f = _ascend(x -> viol(x)[j], dict.points[s]; iters, fd)
            f > best[j] || continue
            idx = add_point!(dict, hdata, model, ξ, :refined)
            for i in 1:nc                       # the new point may raise other constraints
                v = _compact_violation(dict, i, dict.μ[idx], dict.M[idx], gv)
                v > best[i] && (best[i] = v; arg[i] = idx)
            end
        end
    end
    return OracleResult(best, arg, argmax(best))
end
