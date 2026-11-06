using Test
isdefined(@__MODULE__, :STUDIES) || include(joinpath(@__DIR__, "..", "reduce.jl"))

@testset "RLC surrogates" begin
    inputs = reduction_inputs("cases/ieee30-full.toml")
    factors = factor_study(inputs, "/tmp/hpm-surrogate-tests")
    fresh = generate_surrogates(factors; destination = "/tmp/hpm-surrogate-regenerated")
    for path in keys(fresh.files)
        @test read(joinpath(STUDIES, path)) ==
              read(joinpath("/tmp/hpm-surrogate-regenerated", path))
    end
    source = study("cases/ieee30-full.toml")
    factor = study("cases/ieee30-full-factor.toml")
    scaled = study("cases/ieee30-full-scaled.toml")
    p, q, r = nominal_point(source.box), nominal_point(factor.box),
    nominal_point(scaled.box)
    @test length(factor.box) == length(scaled.box) == 6
    @test factor.box.lower == scaled.box.lower
    @test factor.box.upper == scaled.box.upper
    for b in circuits(source), f in (50.0, 650.0)

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
        @test evidence["peak_ratio"] > max(evidence["lower_ratio"], evidence["upper_ratio"])
        @test evidence["screen_lower"] < evidence["operating_screen_resistance"] <
              evidence["screen_upper"]
    end
    dataset = readrecord(joinpath(STUDIES, factor.config.factor_data))
    for result in factors.results
        @test isempty(result.checked.escaped)
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

            interval = getproperty(dataset.coordinates, Symbol("k", quantity, "_", result.row.spec.family_name))
            values = [result.data.observed[3(j - 1) + i, :];
                      result.checked.observed[3(j - 1) + i, :]]
            @test all(x -> interval.lower <= x <= interval.upper, values)
        end
    end
end
