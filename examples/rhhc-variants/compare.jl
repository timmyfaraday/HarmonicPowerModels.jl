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

function comparison(io, paths...)
    println(io, "Uncertainty levels and transfer to the source case")
    println(io, "HPM commit: ", readchomp(`git -C $STUDIES rev-parse HEAD`))
    println(io, "LCM commit: ", readchomp(`git -C $(pkgdir(LineCableModels)) rev-parse HEAD`))
    println(io, "Cases: ", join(paths, ", "))
    println(io, "Threads: ", Threads.nthreads(), "; solver: Clarabel; unchanged runner settings.")
    println(io, "Method parameters: ", METHOD)
    println(io,
        "level | nominal capacity | reported lower bound | upper bound | gap | iterations | LCM stored pairs | solve seconds | signed off | termination")
    flush(io)
    cases = [study(path) for path in paths]
    runs = NamedTuple[]
    for case in cases
        label = splitext(basename(case.file))[1]
        println("Solving ", label)
        flush(stdout)
        own_report = joinpath(STUDIES, "output", label * ".txt")
        result = open(own_report, "w") do own_io
            runstudy(case, own_io)
        end
        push!(runs, result)
        res = result.robust
        println(io, label, " | ", result.nominal["upper_bound"],
            " | ", res["lower_bound"], " | ",
            res["upper_bound"], " | ", res["gap"], " | ", res["iterations"], " | ",
            result.solves, " | ", result.elapsed, " | ", res["signed_off"],
            " | ", res["termination_status"])
        flush(io)
    end
    println(io, "level | harmonic | robust / nominal | robust current")
    for (case, result) in zip(cases, runs), h in case.config.harmonics

        println(io, splitext(basename(case.file))[1], " | ", h, " | ",
            result.robust["solution"][h] / result.nominal["master_solution"][h], " | ",
            result.robust["solution"][h])
    end
    println(io, "Posterior checks: 1024 random points, public HPM defaults (seed 2, tolerance 0.0001).")
    println(io, "solution | target case | maximum normalized violation | constraint | violated fraction")
    for i in eachindex(cases), j in eachindex(cases)

        checked = posterior_check(runs[j].hdata, cases[j].model, cases[j].box,
            runs[i].robust["solution"]; n = 1024)
        println(io, splitext(basename(cases[i].file))[1],
            " | ", splitext(basename(cases[j].file))[1],
            " | ", checked.max_violation, " | ", checked.constraint,
            " | ", checked.violated_fraction)
        flush(io)
    end
    println(io, "Corner adversary: nominal finite-difference signs, unique predicted corners.")
    println(io, "Positive own-case violations refute reported feasibility and its lower-bound certificate.")
    println(io, "Positive transfer violations refute feasibility of the level's solution on the source box.")
    println(io, "A transfer failure alone does not refute an own-case certificate.")
    println(io, "The adversary is not a proof of feasibility: resonance constraints can be nonmonotone.")
    println(io,
        "solution | target case | corners evaluated | maximum normalized violation | constraint | case-coordinate point")
    for i in eachindex(cases), j in eachindex(cases)

        println("Corner adversary: ", basename(cases[i].file), " on ", basename(cases[j].file))
        flush(stdout)
        result = corner_adversary(runs[j].hdata, cases[j].model, cases[j].box, runs[i].robust["solution"])
        println(io, splitext(basename(cases[i].file))[1],
            " | ", splitext(basename(cases[j].file))[1],
            " | ", result.count, " | ", result.max_violation,
            " | ", result.constraint, " | ", result.point)
        flush(io)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: compare.jl cases/<source>.toml cases/<level>.toml ...")
    path = joinpath(STUDIES, "output", splitext(basename(first(ARGS)))[1] *
                                       "-comparison.txt")
    open(path, "w") do io
        comparison(io, ARGS...)
    end
end
