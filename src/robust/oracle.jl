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

# Projected-gradient ascent of f over [-1, 1]^d from ξ0: the gradient from `grad` (ξ ↦ ∇f)
# or, without it, from central finite differences (step `fd`); a normalised step halved
# until f improves; stops after `iters` steps or when no step improves.
function _ascend(f, ξ0::AbstractVector; iters::Int, fd::Real, grad = nothing)
    ξ, fξ = collect(Float64, ξ0), f(ξ0)
    d = length(ξ)
    for _ in 1:iters
        if grad === nothing
            g = zeros(d)
            _tforeach(1:d) do i                 # the 2d evaluations are independent
                e = zeros(d); e[i] = fd
                g[i] = (f(clamp.(ξ .+ e, -1, 1)) - f(clamp.(ξ .- e, -1, 1))) / (2fd)
            end
        else
            g = grad(ξ)
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
               margin = Inf, ce = nothing, gradient = :adjoint, branch_derivative = :central) -> OracleResult

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

The ascent's gradient is analytic (`gradient = :adjoint`, see
[`fast_violation_gradient`](@ref)) or by central finite differences of step `fd`
(`gradient = :fd`). `branch_derivative = :central` selects central differences of the
branch data; `:dual` uses ForwardDiff, and needs branch data that carry dual numbers through
(see [`fast_violation_gradient`](@ref)). It has an effect only with `gradient = :adjoint`.
Each ascent evaluates only its own constraint, through a [`ConstraintEvaluator`](@ref)
(`ce`, built here when not given): its harmonics and the transfer-impedance rows it needs,
with the admittance of the fixed network elements cached. The ascents, and the finite
differences within each ascent step, run in parallel on the available threads; the result
does not depend on the number of threads.

The effort parameters are `ncons`, `margin`, `starts` and `iters`; the sign-off of the
cutting-set method should run at maximum effort (Mutapcic & Boyd 2009, §5.1).
"""
function pessimize!(dict::ParameterDictionary, hdata::Dict, model::AbstractNetworkModel, g;
                    ncons::Int = 5, starts::Int = 3, iters::Int = 10, fd::Real = 1e-3,
                    margin::Real = Inf, ce = nothing, gradient::Symbol = :adjoint,
                    branch_derivative = :central)
    _check_branch_derivative(branch_derivative)
    gv  = [g[h] for h in dict.harmonics]
    V   = dictionary_violations(dict, g)
    nc  = length(dict.keys)
    best = [maximum(view(V, j, :)) for j in 1:nc]
    arg  = [argmax(view(V, j, :)) for j in 1:nc]

    # each ascent evaluates its own constraint only (harmonics and impedance rows it needs)
    ce === nothing && (ce = ConstraintEvaluator(hdata, model, dict.box))

    # the ascents are independent: run them in parallel, then add the improving points in
    # the order of `pairs`, so the result does not depend on the scheduling
    near  = [j for j in sortperm(best; rev = true) if best[j] >= -margin]
    pairs = [(j, s) for j in near[1:min(ncons, length(near))]
                    for s in sortperm(V[j, :]; rev = true)[1:min(starts, size(V, 2))]]
    asc   = _tmap(pairs) do (j, s)
        _ascend(x -> fast_violation(ce, j, x, g), dict.points[s]; iters, fd,
                grad = gradient === :adjoint ? (x -> last(fast_violation_gradient(ce, j, x, g; branch_derivative))) : nothing)
    end
    new   = [ξ for ((j, _), (ξ, f)) in zip(pairs, asc) if f > best[j]]
    for idx in add_points!(dict, hdata, model, new, :refined)
        for i in 1:nc                           # a new point may raise every constraint
            v = _compact_violation(dict, i, dict.μ[idx], dict.M[idx], gv)
            v > best[i] && (best[i] = v; arg[i] = idx)
        end
    end
    return OracleResult(best, arg, argmax(best))
end
