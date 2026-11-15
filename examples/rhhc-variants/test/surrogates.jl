using Test
isdefined(@__MODULE__, :reduction_inputs) || include(joinpath(@__DIR__, "..", "reduce.jl"))

mktempdir() do directory
    @testset "RLC surrogates" begin
        inputs = reduction_inputs("cases/ieee30-full-exact.toml")
        rederive = get(ENV, "RHHC_REDERIVE", "0") == "1"
        evidence = rederive ? joinpath(directory, "evidence") :
                   joinpath(STUDIES, "data", "reduction", "ieee30-full-exact")
        factors = rederive ? factor_study(inputs, evidence) :
                  stored_factors(inputs, evidence)
        fresh = generate_surrogates(factors; destination = joinpath(directory, "regenerated"))
        for path in keys(fresh.files)
            @test read(joinpath(STUDIES, path)) ==
                  read(joinpath(joinpath(directory, "regenerated"), path))
        end
        source = study("cases/ieee30-full-exact.toml")
        factor = study("cases/ieee30-full-exact-factor.toml")
        scaled = study("cases/ieee30-full-exact-scaled.toml")
        p, q, r = nominal_point(source.box), nominal_point(factor.box),
        nominal_point(scaled.box)
        @test length(factor.box) == length(scaled.box) == 6
        @test factor.box.lower == scaled.box.lower
        @test factor.box.upper == scaled.box.upper
        for b in circuits(source), f in (50.0, 650.0, 1250.0)

            a = source.model.branches[b].zy(p, f)
            c = factor.model.branches[b].zy(q, f)
            for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

                @test isequal(a[sequence][quantity], c[sequence][quantity])
            end
        end
        for b in uncertain_branches(source.model)
            @test branch_data(source.model, b, p, 1.0, :positive) ==
                  branch_data(scaled.model, b, r, 1.0, :positive)
        end
        for b in circuits(source)
            family = string(specification(source.details[b].spec).family_name)
            name = Symbol("kR_", family)
            index = findfirst(==(name), factor.box.names)
            value = factor.box.upper[index]
            point = merge(q, NamedTuple{(name,)}((value,)))
            reference = source.model.branches[b].zy(p, 650.0, :positive)
            changed = factor.model.branches[b].zy(point, 650.0, :positive)
            @test real(changed.z) == value * real(reference.z)
            @test imag(changed.z) == imag(reference.z)
            @test changed.y == reference.y
            base = branch_data(source.model, b, p, 1.0, :positive)
            data = branch_data(scaled.model, b, point, 13.0, :positive)
            @test data.br_r == value * base.br_r * sqrt(13.0)
            @test data.br_x == base.br_x * 13.0
            @test data.b_fr == base.b_fr * 13.0
        end
        @test !isempty(fresh.data["screen_loss_evidence"])
        for evidence in fresh.data["screen_loss_evidence"]
            @test evidence["peak_ratio"] >
                  max(evidence["lower_ratio"], evidence["upper_ratio"])
            @test evidence["screen_lower"] < evidence["operating_screen_resistance"] <
                  evidence["screen_upper"]
        end
        dataset = readrecord(joinpath(STUDIES, factor.config.factor_data))
        points = [latent(inputs.box, p)
                  for p in background_points(inputs.box, 256; seed = 3)]
        for result in factors.results
            @test isempty(result.checked.escaped)
            observed = rederive ?
                       hcat(result.data.observed, result.checked.observed) :
                       reduce(hcat,
                [last.(observations(result.row.problem(p, 50.0), result.row))
                 for p in [source_points(inputs.box); points]])
            for (k, interval) in enumerate(result.hulls.bounds)
                @test all(x -> interval.lower <= x <= interval.upper, observed[k, :])
            end
            for k in eachindex(result.hulls.cells)
                bad = filter(c -> !c.consistent, result.hulls.cells[k])
                isempty(bad) && continue
                @test all(
                    c -> c.depth == METHOD.subdivision_depth &&
                         c.method == "interior search", bad)
                checks = filter(c -> c.metric == k, result.interior)
                @test sum(c.samples for c in checks) >= METHOD.verification_samples
                @test all(c -> c.escapes == 0 && c.largest <= 0, checks)
            end
            for (j, f) in enumerate(factors.frequencies),
                (i, quantity) in enumerate(("R", "L", "C"))

                interval = getproperty(dataset.coordinates,
                    Symbol("k", quantity, "_", result.row.spec.family_name))
                values = observed[3(j - 1) + i, :]
                @test all(x -> interval.lower <= x <= interval.upper, values)
            end
        end
    end
end
