using Test
include(joinpath(@__DIR__, "..", "builders.jl"))

@testset "Exact variants preserve the network" begin
    files = sort(filter(endswith("-exact.toml"), readdir(joinpath(STUDIES, "cases"))))
    @test length(files) == 5
    mktempdir() do directory
        for file in files
            @testset "$file" begin
                rows = TOML.parsefile(joinpath(STUDIES, "cases", file))
                # The base uses the same loaded tables for this network comparison.
                base_rows = TOML.parsefile(joinpath(STUDIES, "cases", rows["base"] *
                                                                      ".toml"))
                base_rows["angle_tables"] = rows["angle_tables"]
                path = joinpath(directory, "base.toml")
                open(path, "w") do io
                    TOML.print(io, base_rows)
                end
                base = study(path)
                variant = study("cases/" * file)
                for field in (:names, :nominal, :lower, :upper, :scale)
                    @test getproperty(variant.box, field) == getproperty(base.box, field)
                end
                @test variant.hd["nw"]["1"] == base.hd["nw"]["1"]
                @test variant.details == base.details
                @test typeof(variant.model) == typeof(base.model)
                @test Set(uncertain_branches(variant.model)) ==
                      Set(uncertain_branches(base.model))
                p = nominal_point(base.box)
                for b in uncertain_branches(base.model), h in (1, 5, 25),
                    sequence in (:zero, :positive, :negative)
                    @test branch_data(variant.model, b, p, h, sequence) ===
                          branch_data(base.model, b, p, h, sequence)
                end
                harmonics = [5, 7, 11, 13, 17, 19, 23, 25]
                @test variant.config.harmonics == Tuple(harmonics)
                @test sort(collect(keys(variant.pce["harmonics"]))) == harmonics
                @test sort(parse.(Int, collect(keys(variant.hd["nw"])))) == [1; harmonics]
                @test variant.pce["moments"] == "exact"
            end
        end
    end
end

@testset "Exact uncertainty levels" begin
    source = study("cases/ieee30-full-exact.toml")
    for level in ("equivalent", "factor", "scaled")
        s = study("cases/ieee30-full-exact-$level.toml")
        @test s.config.base == "ieee30-full-exact"
        @test s.config.harmonics == source.config.harmonics
        @test s.config.angle_tables == source.config.angle_tables
        @test s.config.branches == source.config.branches
        @test s.pce == source.pce
        @test s.pce["moments"] == "exact"
        @test Set(keys(s.hd["nw"])) == Set(keys(source.hd["nw"]))
        @test Set(uncertain_branches(s.model)) == Set(uncertain_branches(source.model))
        @test length(s.box) == (level == "equivalent" ? 35 : 6)
        @test all(s.details[b].len == source.details[b].len for b in keys(source.details))
        mktempdir() do directory
            file = joinpath(directory, "mix.toml")
            open(file, "w") do io
                TOML.print(io, Dict("base" => "ieee30-full-exact-$level", "population" => "mix"))
            end
            mixed = study(file)
            @test mixed.pce["moments"] == "exact"
            @test mixed.pce["bound"] == source.pce["bound"]
            @test length(mixed.pce["unit_shapes"]) == length(source.pce["harmonic_units"])
            @test mixed.box.names == s.box.names
            @test mixed.hd == s.hd
            for b in uncertain_branches(s.model), h in (1, 5, 25)

                sequence = harmonic_sequence(h)
                @test branch_data(mixed.model, b, nominal_point(mixed.box), h, sequence) ===
                      branch_data(s.model, b, nominal_point(s.box), h, sequence)
            end
        end
        if level == "equivalent"
            Ext = Base.get_extension(HPM, :HarmonicPowerModelsLineCableModelsExt)
            points = (
                nominal_point(s.box),
                NamedTuple{Tuple(s.box.names)}(Tuple(s.box.upper))
            )
            for b in circuits(s), p in points

                c = s.model.branches[b].zy
                full = Ext.Circuit(c.problem, c.formulation, s.box)
                @test c(p, 1250.0) === compute(full, p, 1250.0)
                for name in s.box.names
                    name in s.details[b].coordinates && continue
                    value = getproperty(p, name) == 1.0 ? -1.0 : 1.0
                    point = merge(p, NamedTuple{(name,)}((value,)))
                    @test compute(full, point, 1250.0) === c(p, 1250.0)
                    @test c(point, 1250.0) === c(p, 1250.0)
                end
            end
        end
    end
end

include("equivalent.jl")
include("surrogates.jl")
