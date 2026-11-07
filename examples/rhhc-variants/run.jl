isdefined(@__MODULE__, :STUDIES) || include(joinpath(@__DIR__, "builders.jl"))

function report(io, args...)
    println(io, args...)
    println(args...)
    flush(io)
end

function runstudy(s, io)
    lines = IOBuffer()
    header(lines, s)
    declarations(lines, s)
    text = String(take!(lines))
    print(io, text)
    print(text)
    source = s.config.network
    label = source.kind == "scenario" ? "IEEE 30-bus, " * source.scenario : source.file
    report(io, label, ", ", length(s.details),
        " branches from LineCableModels, ", length(s.pce["harmonic_units"]), " units, ε = ",
        s.config.epsilon, ", Cantelli, ", Threads.nthreads(), " threads")
    box = s.box
    box0 = ParameterBox(box.names, box.nominal, box.nominal, box.nominal, box.scale)
    nom = solve_rhhc(deepcopy(s.hd), s.model, box0, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt,
        nbackground = 0, use_seeds = false)
    report(io, @sprintf("nominal stochastic: %.6f (%s)",
        nom["upper_bound"], nom["termination_status"]))
    isempty(nom["solution"]) && report(io, "nominal solution: ", nom["solution"])
    flush(io)
    report(io, "resonances at the nominal point: ",
        join(
            [@sprintf("h = %.2f", r.h)
             for r in modal_scan(deepcopy(s.hd), s.model, nominal_point(box))],
            ", "))
    flush(io)
    prepared = deepcopy(s.hd)
    t = @elapsed res = solve_rhhc(prepared, s.model, box, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt,
        posterior = s.config.posterior, verbose = true)
    report(io,
        @sprintf("robust: %.6f  [LB %.6f, UB %.6f, gap %.1e], %d iterations, %.0f s, %s",
            res["lower_bound"], res["lower_bound"], res["upper_bound"], res["gap"],
            res["iterations"], t, res["signed_off"] ? "signed off" : "NOT signed off"))
    report(io, "robust master termination: ", res["termination_status"])
    isempty(res["solution"]) && report(io, "solution: ", res["solution"])
    if isempty(res["solution"]) || isempty(nom["master_solution"])
        report(io, "Harmonic ratios not evaluated: an empty solution has no harmonic currents.")
    else
        report(io, "per harmonic, robust / nominal: ",
            join(
                [@sprintf("h%d %.3f", h, res["solution"][h] / nom["master_solution"][h])
                 for h in s.config.harmonics], ", "))
    end
    if isempty(res["solution"])
        report(io, "Posterior check not evaluated: the solution is empty.")
    else
        pr = res["posterior"]
        report(io,
            @sprintf("a posteriori (%d random points): worst %+.2e (%s), violated at %.1f %%",
                s.config.posterior, pr.max_violation, string(pr.constraint),
                100pr.violated_fraction))
    end
    solves = lcm_solves(s)
    report(io, "LineCableModels solves (sum of circuit table sizes): ", solves)
    flush(io)
    return (; robust = res, nominal = nom, hdata = prepared, elapsed = t, solves)
end

function main(path)
    mkpath(joinpath(STUDIES, "output"))
    dest = joinpath(STUDIES, "output", splitext(basename(path))[1] * ".txt")
    failed = open(dest, "w") do io
        try
            s = study(path)
            runstudy(s, io)
            return false
        catch err
            report(io, "Run failed:")
            showerror(io, err)
            println(io)
            showerror(stderr, err, catch_backtrace())
            println(stderr)
            return true
        end
    end
    return failed
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run.jl cases/<case>.toml")
    main(only(ARGS)) && exit(1)
end
