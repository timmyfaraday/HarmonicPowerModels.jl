@testset "Circuit coordinate declarations" begin
    Ext = Base.get_extension(HPM, :HarmonicPowerModelsLineCableModelsExt)
    for file in sort(readdir(joinpath(STUDIES, "cases")))
        endswith(file, ".toml") || continue
        haskey(TOML.parsefile(joinpath(STUDIES, "cases", file)), "base") && continue
        @testset "$file" begin
            s = study(joinpath("cases", file); all_coordinates = true)
            p = nominal_point(s.box)
            for b in circuits(s)
                c = s.model.branches[b].zy
                declared = s.details[b].coordinates
                full = Ext.Circuit(c.problem, c.formulation, s.box)
                for f in (50.0, 650.0)
                    reference = compute(full, p, f)
                    @test c(p, f) === reference
                    for (i, name) in enumerate(s.box.names)
                        name in declared && continue
                        for value in (-1.0, -0.5, 0.5, 1.0)
                            ξ = zeros(length(s.box))
                            ξ[i] = value
                            point = merge(p, NamedTuple{(name,)}((getproperty(latent(s.box, ξ), name),)))
                            # Bypass the table to test what the problem actually reads.
                            @test compute(full, point, f) === reference
                            @test c(point, f) === reference
                        end
                    end
                end
            end
        end
    end
end

@testset "Installed length" begin
    baseline = study("cases/ieee30-physical.toml")
    enabled = study("cases/ieee30-physical-length.toml")
    p = nominal_point(enabled.box)
    q = nominal_point(baseline.box)
    ids = sort([parse(Int, string(name))
                for (name, row) in pairs(enabled.config.branches) if haskey(row, :length_m)])
    names = [Symbol("length_", b) for b in ids]
    @test length(ids) == 22
    @test enabled.box.names == [baseline.box.names; names]
    for (i, b) in enumerate(ids)
        row = getproperty(enabled.config.branches, Symbol(b)).length_m
        name = names[i]
        index = length(baseline.box) + i
        @test enabled.box.nominal[index] === row.value
        @test enabled.box.lower[index] === row.lower
        @test enabled.box.upper[index] === row.upper
        @test enabled.box.scale[index] === :linear
        @test enabled.details[b].len === baseline.details[b].len
        @test name ∉ enabled.details[b].coordinates
        @test enabled.model.branches[b].len(p) === getproperty(p, name)
        for value in (row.lower, row.upper)
            point = merge(p, NamedTuple{(name,)}((value,)))
            @test enabled.model.branches[b].len(point) === value
            for sequence in (:zero, :positive, :negative)
                @test branch_data(enabled.model, b, point, 13, sequence) !=
                      branch_data(enabled.model, b, p, 13, sequence)
            end
        end
    end

    directory = mktempdir()
    try
        config = TOML.parsefile(joinpath(STUDIES, "cases", "ieee30-physical-length.toml"))
        for row in values(config["branches"])
            haskey(row, "length_m") || continue
            delete!(row["length_m"], "lower")
            delete!(row["length_m"], "upper")
        end
        path = joinpath(directory, "length-fixed.toml")
        open(path, "w") do io
            TOML.print(io, config)
        end
        fixed = study(path)
        for field in (:names, :nominal, :lower, :upper, :scale)
            @test getproperty(fixed.box, field) == getproperty(baseline.box, field)
        end
        for b in keys(baseline.details), h in (1.0, 13.0),
            sequence in (:zero, :positive, :negative)
            reference = branch_data(baseline.model, b, q, h, sequence)
            @test branch_data(fixed.model, b, q, h, sequence) === reference
            @test branch_data(enabled.model, b, p, h, sequence) === reference
        end

        # A fixed value overrides the source route without adding a coordinate.
        config["branches"]["17"]["length_m"]["value"] *= 1.01
        open(path, "w") do io
            TOML.print(io, config)
        end
        changed = study(path)
        @test changed.model.branches[17].len ===
              config["branches"]["17"]["length_m"]["value"]
        @test changed.box.names == baseline.box.names
        @test branch_data(changed.model, 17, q, 1, :positive) !=
              branch_data(baseline.model, 17, q, 1, :positive)

        config["branches"]["17"]["route"] = "retained"
        open(path, "w") do io
            TOML.print(io, config)
        end
        @test_throws ArgumentError("branch 17 declares both route and length_m") study(path)
    finally
        rm(directory; recursive = true)
    end
end
