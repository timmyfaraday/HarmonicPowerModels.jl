include(joinpath(@__DIR__, "run.jl"))

"""
    constraint_signs(hdata, model, box, solution)

Differentiate every positive-limit normalized constraint at the nominal point.
Use the shared half step on the signed latent scale and its sensitivity zero band.
Near a box endpoint the difference is one-sided. Logarithmic axes keep their scale.
"""
function constraint_signs(hdata, model, box, solution)
    nominal = nominal_point(box)
    moments = moments_at(hdata, model, nominal)
    constraints = filter(c -> c.limit > 0, shhc_constraints(hdata, moments))
    keys = [(c.kind, c.key) for c in constraints]
    g = Dict(h => current^2 for (h, current) in solution)
    derivatives = zeros(length(keys), length(box))
    for (j, name) in enumerate(box.names)
        lo, center, hi = box.lower[j], box.nominal[j], box.upper[j]
        lo == hi && continue
        if box.scale[j] == :log
            lo, center, hi = log(lo), log(center), log(hi)
        end
        step = METHOD.difference_step * (hi - lo) / 2
        lower, upper = max(lo, center - step), min(hi, center + step)
        a, b = box.scale[j] == :log ? (exp(lower), exp(upper)) : (lower, upper)
        left = moments_at(hdata, model, merge(nominal, NamedTuple{(name,)}((a,))))
        right = moments_at(hdata, model, merge(nominal, NamedTuple{(name,)}((b,))))
        lc = filter(c -> c.limit > 0, shhc_constraints(hdata, left))
        rc = filter(c -> c.limit > 0, shhc_constraints(hdata, right))
        [(c.kind, c.key) for c in lc] == keys == [(c.kind, c.key) for c in rc] ||
            error("constraint identities changed during the nominal difference")
        for k in eachindex(keys)
            derivatives[k, j] = (normalized_violation(rc[k], g, right) -
                                 normalized_violation(lc[k], g, left)) * (hi - lo) /
                                (2(upper - lower))
        end
    end
    signs = map(x -> abs(x) <= METHOD.zero_band ? 0 : Int(sign(x)), derivatives)
    return (; keys, derivatives, signs)
end

"""
    corner_adversary(hdata, model, box, solution)

Predict one adverse point per positive-limit constraint from its nominal signs.
Deduplicate the points and evaluate each through public HPM `posterior_check` on
a zero-width box. Coordinates with zero sensitivity keep their nominal value.
The result is sampled evidence. Nonlinear constraints can have worse points elsewhere.
"""
function corner_adversary(hdata, model, box, solution)
    sensitivities = constraint_signs(hdata, model, box, solution)
    corners = unique([Tuple(sign > 0 ? box.upper[j] :
                            sign < 0 ? box.lower[j] : box.nominal[j]
                      for (j, sign) in enumerate(row))
                      for row in eachrow(sensitivities.signs)])
    observations = map(corners) do corner
        fixed = ParameterBox(box.names, corner, corner, corner, box.scale)
        posterior_check(hdata, model, fixed, solution; n = 1)
    end
    worst = argmax([result.max_violation for result in observations])
    point = NamedTuple{Tuple(box.names)}(corners[worst])
    return (; max_violation = observations[worst].max_violation,
        constraint = observations[worst].constraint, point,
        count = length(corners), corners, observations, sensitivities)
end

function comparison(io, source_path, equivalent_path)
    println(io, "Source and equivalent uncertainty comparison")
    println(io, "HPM commit: ", readchomp(`git -C $STUDIES rev-parse HEAD`))
    println(io, "LCM commit: ", readchomp(`git -C $(pkgdir(LineCableModels)) rev-parse HEAD`))
    println(io, "Source: ", source_path, "; equivalent: ", equivalent_path)
    println(io, "Threads: ", Threads.nthreads(), "; solver: Clarabel; unchanged runner settings.")
    println(io, "Method parameters: ", METHOD)
    flush(io)
    source = study(source_path)
    equivalent = study(equivalent_path)
    println("Solving source case")
    first_run = runstudy(source, IOBuffer())
    println("Solving equivalent case")
    second_run = runstudy(equivalent, IOBuffer())
    println(io, "quantity | source | equivalent")
    for key in ("lower_bound", "upper_bound", "gap", "iterations", "signed_off", "termination_status")
        println(io, key, " | ", first_run.robust[key], " | ", second_run.robust[key])
    end
    println(io, "LCM solves (stored point-frequency pairs) | ",
        first_run.solves, " | ", second_run.solves)
    println(io, "solve seconds | ", first_run.elapsed, " | ", second_run.elapsed)
    println(io,
        "harmonic | source robust / nominal | equivalent robust / nominal | source current | equivalent current")
    for h in source.config.harmonics
        println(io, h, " | ",
            first_run.robust["solution"][h] / first_run.nominal["master_solution"][h],
            " | ", second_run.robust["solution"][h] /
                   second_run.nominal["master_solution"][h],
            " | ", first_run.robust["solution"][h], " | ", second_run.robust["solution"][h])
    end
    flush(io)
    println(io, "Posterior checks: 1024 random points, public HPM defaults (seed 2, tolerance 0.0001).")
    println(io, "solution | case | maximum normalized violation | constraint | violated fraction")
    combinations = (("source", "source", first_run, source, first_run),
        ("equivalent", "equivalent", second_run, equivalent, second_run),
        ("equivalent", "source", second_run, source, first_run),
        ("source", "equivalent", first_run, equivalent, second_run))
    for (label, target, solved, case, prepared) in combinations
        checked = posterior_check(prepared.hdata, case.model, case.box,
            solved.robust["solution"]; n = 1024)
        println(io, label, " | ", target, " | ",
            checked.max_violation, " | ", checked.constraint,
            " | ", checked.violated_fraction)
        flush(io)
    end
    println(io, "Corner adversary: nominal finite-difference signs, unique predicted corners.")
    println(io, "Positive values violate the constraint. The sign test is not a proof of monotonicity.")
    println(io, "A source-on-source violation refutes the source solution's reported feasibility.")
    println(io, "An equivalent-on-equivalent violation refutes the equivalent solution's reported feasibility.")
    println(io, "An equivalent-on-source violation refutes its transfer to the source case.")
    println(io, "A source-on-equivalent violation refutes its transfer to the independent equivalent box.")
    println(io, "Cross-case violations alone do not refute a solution's own-case certificate.")
    println(io,
        "solution | case | corners evaluated | maximum normalized violation | constraint | physical corner")
    for (label, target, solved, case, prepared) in combinations
        println("Corner adversary: ", label, " on ", target)
        adversary = corner_adversary(prepared.hdata, case.model, case.box, solved.robust["solution"])
        println(io, label, " | ", target, " | ",
            adversary.count, " | ", adversary.max_violation,
            " | ", adversary.constraint, " | ", adversary.point)
        flush(io)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 ||
        error("usage: compare.jl cases/<source>.toml cases/<equivalent>.toml")
    path = joinpath(STUDIES, "output", splitext(basename(first(ARGS)))[1] *
                                       "-comparison.txt")
    open(path, "w") do io
        comparison(io, ARGS...)
    end
end
