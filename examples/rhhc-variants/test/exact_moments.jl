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
