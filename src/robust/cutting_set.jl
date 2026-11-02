################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): the cutting-set method.             #
#                                                                               #
# The robust problem enforces the chance constraints of sHHC_SOC at every     #
# point of the parameter box. The cutting-set method (Mutapcic & Boyd 2009,    #
# §5.1) alternates                                                             #
#   optimization:  the master, i.e. sHHC_SOC with the constraints of a finite  #
#                  scenario set per constraint (starting at the nominal point),#
#   pessimization: the approximate oracle `pessimize!`,                        #
# and adds, for every violated constraint, its worst point to its scenarios.   #
#                                                                               #
# Bounds. The master relaxes the robust problem, so its value is an upper      #
# bound (§3.1). g = 0 is a Slater point whenever every limit is positive, so   #
# scaling the master solution down makes it feasible (§3.2); as every chance   #
# constraint is positively homogeneous of degree 1 in g, the scaling           #
# g/(1+v), i.e. I/√(1+v), with v the largest normalized violation, suffices.   #
# With an approximate oracle, v and the lower bound hold relative to the       #
# points evaluated (§5.1); the a-posteriori check is the evidence beyond them. #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

"""
    solve_rhhc(hdata, model, box, optimizer; pce_params, hpf_optimizer = optimizer,
               nbackground = 1024, use_seeds = true, effort = (ncons = 5, starts = 3, iters = 10),
               signoff = (ncons = typemax(Int), starts = 5, iters = 30), tol = 1e-4,
               maxiter = 50, verbose = false) -> Dict

Robust (stochastic-angle, worst-case-network) harmonic hosting capacity by the cutting-set
method, for the uncertain branches of `model` over the parameter box `box`.

- `pce_params` as for [`solve_shhc`](@ref) with `sHHC_SOC`; the bound must be `"cantelli"` or
  `"gaussian"` (the empirical λ* is calibrated on the nominal network only).
- The dictionary of the oracle holds the nominal point, the vertices (d ≤ 10), the modal seeds
  (`use_seeds`) and `nbackground` Sobol points.
- Each iteration runs the oracle at `effort`; when it finds no violation above `tol`, the
  sign-off reruns it at `signoff` effort before stopping.

Returns a Dict with `"termination_status"` of the last master, `"upper_bound"` (the master
value), `"lower_bound"` and `"solution"` (I_h of the scaled, robust-feasible decision),
`"master_solution"` (the master's I_h), `"max_violation"` (the oracle's last v at the master
solution), `"gap"`, `"iterations"`, `"history"` (per iteration: v, bounds, cuts added, time)
and `"dictionary"`.
"""
function solve_rhhc(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox, optimizer;
                    pce_params::Dict, hpf_optimizer = optimizer, nbackground::Int = 1024,
                    use_seeds::Bool = true,
                    effort = (ncons = 5, starts = 3, iters = 10),
                    signoff = (ncons = typemax(Int), starts = 5, iters = 30),
                    tol::Real = 1e-4, maxiter::Int = 50, verbose::Bool = false)
    principle = get(hdata, "principle", "absolute equality")
    principle in ("absolute equality", "weighted equality") || error(
        "solve_rhhc: fairness principle \"$principle\" is not representable in the second-" *
        "order cone program; use \"absolute equality\" or \"weighted equality\".")
    t0 = time()

    # fundamental power flow and PCE data, as sHHC_SOC (fundamental held at nominal)
    update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, hpf_optimizer)
    update_hdata_with_pce_data!(hdata, pce_params)
    _check_robust_bound(hdata)

    # master: sHHC_SOC at the nominal point
    m0     = moments_at(hdata, model, nominal_point(box))
    master = shhc_master(hdata, m0, optimizer)
    c0     = shhc_constraints(hdata, m0)
    for c in c0; add_shhc_constraint!(master, c, m0); end
    (; I_h, S) = master
    jm = master.model
    H = sort(m0.harm_nw_ints)

    # oracle dictionary
    seeds = use_seeds ? resonance_seeds(hdata, model, box) : Seed[]
    dict  = build_dictionary(hdata, model, box; seeds, nbackground)
    pos   = dict.limit .> 0.0

    # scenario sets: point indices per constraint; the nominal point (index 1) holds already
    scen   = [Set{Int}(1) for _ in dict.keys]
    mcache = Dict{Int,Any}(1 => m0)
    history = NamedTuple[]
    Iv = Dict{Int,Float64}(); v = Inf; UB = NaN; LB = NaN; status = nothing
    signed_off = false

    for it in 1:maxiter
        JuMP.optimize!(jm)
        status = JuMP.termination_status(jm)
        JuMP.has_values(jm) || break
        UB = S * JuMP.objective_value(jm)
        Iv = Dict(h => max(JuMP.value(I_h[h]), 0.0) for h in H)
        g  = Dict(h => Iv[h]^2 for h in H)

        r  = pessimize!(dict, hdata, model, g; effort...)
        v  = maximum(r.violation[pos]; init = 0.0)
        if v <= tol
            r = pessimize!(dict, hdata, model, g; signoff...)     # sign-off at maximum effort
            v = maximum(r.violation[pos]; init = 0.0)
            signed_off = v <= tol
        end
        LB = UB / sqrt(1.0 + max(v, 0.0))

        ncut = 0
        if !signed_off
            for j in eachindex(dict.keys)
                r.violation[j] > tol || continue
                s = r.point[j]
                s in scen[j] && continue
                m = get!(() -> moments_at(hdata, model, latent(box, dict.points[s])), mcache, s)
                c = shhc_constraints(hdata, m)[j]
                add_shhc_constraint!(master, c, m)
                push!(scen[j], s); ncut += 1
            end
        end
        push!(history, (iteration = it, max_violation = v, upper_bound = UB, lower_bound = LB,
                        cuts = ncut, time = time() - t0))
        verbose && @info "solve_rhhc" iteration=it max_violation=v upper_bound=UB lower_bound=LB cuts=ncut
        (signed_off || ncut == 0) && break
    end
    ncut_last = isempty(history) ? 0 : last(history).cuts
    !signed_off && ncut_last == 0 && v > tol &&
        @warn "solve_rhhc: violations above tol remain at points already in the scenario sets " *
              "(numerical tolerance of the master); the lower bound accounts for them."

    scale = 1.0 / sqrt(1.0 + max(v, 0.0))
    return Dict(
        "termination_status" => status,
        "upper_bound"        => UB,
        "lower_bound"        => LB,
        "gap"                => isfinite(UB) && UB > 0 ? (UB - LB) / UB : NaN,
        "solution"           => Dict(h => scale * I for (h, I) in Iv),
        "master_solution"    => Iv,
        "max_violation"      => v,
        "signed_off"         => signed_off,
        "iterations"         => length(history),
        "history"            => history,
        "scenarios"          => Dict(dict.keys[j] => sort(collect(scen[j])) for j in eachindex(scen)
                                     if length(scen[j]) > 1),
        "dictionary"         => dict,
        "solve_time"         => time() - t0,
    )
end
