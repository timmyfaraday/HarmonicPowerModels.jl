@testset "Converter populations" begin
    reference = study("cases/ieee30-tom.toml")
    units = sort(reference.pce["harmonic_units"])
    expected = Dict{String, Any}(deepcopy(PCE_PARAMS)...)
    expected["recombination"] = true
    expected["harmonics"] = Dict{Int, Any}(h => Tuple(PCE_PARAMS["harmonics"][h])
    for h in reference.config.harmonics)
    expected["harmonic_units"] = reference.pce["harmonic_units"]
    expected["bound"] = reference.config.bound
    expected["epsilon"] = reference.config.epsilon
    @test reference.pce == expected
    @test !haskey(reference.pce, "unit_shapes")

    mktempdir() do directory
        path = joinpath(directory, "population.toml")
        for population in ("diode", "thyristor", "mix")
            rows = Dict("base" => "ieee30-tom",
                "population" => Dict("value" => population, "origin" => "test"))
            open(path, "w") do io
                TOML.print(io, rows)
            end
            s = study(path)
            if population == "diode"
                @test s.pce == reference.pce
                @test !haskey(s.pce, "unit_shapes")
            elseif population == "thyristor"
                @test s.pce["alpha_hat"] == PCE_PARAMS_THYRISTOR["alpha_hat"]
                @test s.pce["beta_hat"] == PCE_PARAMS_THYRISTOR["beta_hat"]
                for h in s.config.harmonics
                    @test s.pce["harmonics"][h] ==
                          Tuple(PCE_PARAMS_THYRISTOR["harmonics"][h])
                end
            else
                @test count(
                    d -> s.pce["unit_shapes"][d] ==
                         (PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]),
                    units) == 11
                @test count(
                    d -> s.pce["unit_shapes"][d] ==
                         (PCE_PARAMS_THYRISTOR["alpha_hat"],
                        PCE_PARAMS_THYRISTOR["beta_hat"]),
                    units) == 10
                for (i, d) in enumerate(units)
                    params = isodd(i) ? PCE_PARAMS : PCE_PARAMS_THYRISTOR
                    @test s.pce["unit_shapes"][d] ==
                          (params["alpha_hat"], params["beta_hat"])
                    for h in s.config.harmonics
                        @test s.pce["harmonics"][h][d] == Tuple(params["harmonics"][h])
                    end
                end
            end
        end
        # List order does not select the converter type: sorted load ids do.
        rows = Dict("base" => "ieee30-tom", "units" => Dict("value" => [9, 2, 7]),
            "population" => Dict("value" => "mix", "origin" => "test"))
        open(path, "w") do io
            TOML.print(io, rows)
        end
        s = study(path)
        for (i, d) in enumerate([2, 7, 9])
            params = isodd(i) ? PCE_PARAMS : PCE_PARAMS_THYRISTOR
            @test s.pce["unit_shapes"][d] == (params["alpha_hat"], params["beta_hat"])
        end
        # A variant replaces complete top-level records.
        rows["harmonics"] = Dict("value" => [5, 7], "origin" => "test")
        rows["epsilon"] = Dict("value" => 0.1, "origin" => "test")
        open(path, "w") do io
            TOML.print(io, rows)
        end
        s = study(path)
        @test s.config.harmonics == (5, 7)
        @test s.pce["epsilon"] == 0.1
        @test Set(keys(s.pce["harmonics"])) == Set([5, 7])

        for population in ("diode", "thyristor", "mix")
            for extra in (Dict("recombination" => false),
                Dict("harmonics" => Dict("value" => [3, 5])),
                Dict("base" => "ieee30-tom-mix"))
                rows = merge(
                    Dict("base" => "ieee30-tom",
                        "population" => Dict("value" => population, "origin" => "test")),
                    extra)
                open(path, "w") do io
                    TOML.print(io, rows)
                end
                @test_throws ArgumentError study(path)
            end
        end
        for population in ("unknown", "")
            open(path, "w") do io
                TOML.print(io, Dict("base" => "ieee30-tom", "population" =>
                    Dict("value" => population)))
            end
            @test_throws ArgumentError study(path)
        end
    end
    for value in ("thyristor", "unknown", "")
        withenv("SHHC_UNIT" => value) do
            @test_throws ArgumentError study("cases/ieee30-tom.toml")
        end
    end
    withenv("SHHC_UNIT" => "diode") do
        @test study("cases/ieee30-tom.toml").pce == reference.pce
    end
end

@testset "Population variants preserve the network" begin
    files = sort(filter(endswith(".toml"), readdir(joinpath(STUDIES, "cases"))))
    bases = filter(file -> !haskey(TOML.parsefile(joinpath(STUDIES, "cases", file)), "base"), files)
    variants = setdiff(files, bases)
    expected = [replace(file, ".toml" => "-mix.toml")
                for file in bases
                if get(readrecord(joinpath(STUDIES, "cases", file)), :recombination, true)]
    @test variants == sort(expected)
    for file in variants
        @testset "$file" begin
            rows = TOML.parsefile(joinpath(STUDIES, "cases", file))
            @test Set(keys(rows)) == Set(["base", "population"])
            @test rows["population"] ==
                  Dict("value" => "mix", "origin" => "examples/shhc/paper_mix_ab.jl")
            base = study("cases/" * rows["base"] * ".toml")
            variant = study("cases/" * file)
            for field in (:names, :nominal, :lower, :upper, :scale)
                @test getproperty(variant.box, field) == getproperty(base.box, field)
            end
            @test variant.hd == base.hd
            @test variant.details == base.details
            @test typeof(variant.model) == typeof(base.model)
            p = nominal_point(base.box)
            for b in uncertain_branches(base.model), h in (1, 5, 7, 11, 13),
                sequence in (:zero, :positive, :negative)
                if base.config.model == "scaled" && sequence === :zero
                    @test_throws ArgumentError branch_data(base.model, b, p, h, sequence)
                    @test_throws ArgumentError branch_data(variant.model, b, p, h, sequence)
                else
                    @test branch_data(variant.model, b, p, h, sequence) ===
                          branch_data(base.model, b, p, h, sequence)
                end
            end
        end
    end
end
