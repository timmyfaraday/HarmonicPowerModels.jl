@testset "ownership, placement and names probes" begin
    planted, _ = probe("structure/planted", "StructureProbePlanted")
    clean, _ = probe("structure/clean", "StructureProbeClean")
    for (tree, mode) in ((planted, "planted"), (clean, "clean"))
        methods = Base.invokelatest(Guards.package_methods, tree)
        homes = Base.invokelatest(Guards.home_directories, tree)
        ownership = Base.invokelatest(Guards.ownership, methods, tree)
        placement = Base.invokelatest(Guards.placement, methods, tree, homes)
        names = Base.invokelatest(Guards.shared_names, tree)
        if mode == "planted"
            @test ownership == Dict(
                "StructureProbePlanted | branch_data | ext/HarmonicPowerModelsLineCableModelsExt/StructureProbePlanted.jl" =>
                    2,
                "StructureProbePlanted | compute | ext/HarmonicPowerModelsLineCableModelsExt/StructureProbePlanted.jl" =>
                    1,
                "StructureProbePlanted | length | ext/HarmonicPowerModelsLineCableModelsExt/StructureProbePlanted.jl" =>
                    1)
            @test placement ==
                  Dict("test/linecablemodels/misplaced.jl | StructureProbePlanted" => 1)
            @test names == Dict("StructureProbePlanted | length -> Base" => 1,
                "StructureProbePlanted | ParameterBox -> HarmonicPowerModels" => 1,
                "StructureProbePlanted | CoaxialFormulation -> LineCableModels" => 1)
        else
            @test ownership == Dict{String, Int}()
            @test placement == Dict{String, Int}()
            @test names == Dict{String, Int}()
        end
    end
end

@testset "naming probes" begin
    for mode in ("planted", "clean")
        name = mode == "planted" ? "NamingProbePlanted" : "NamingProbeClean"
        tree, _ = probe("naming/" * mode, name)
        methods = Base.invokelatest(Guards.package_methods, tree)
        names = Guards.TOML.parsefile(joinpath(tree.directory, "names.toml"))
        found = Base.invokelatest(Guards.naming_rule, methods, tree, names)
        if mode == "planted"
            @test found == Dict(
                "terms | function | NamingProbePlanted | unknown_name | unknown" => 1,
                "context | function | NamingProbePlanted | line_matrix | line" => 1)
        else
            @test found == Dict{String, Int}()
        end
    end
end

@testset "physical names and constants probes" begin
    for mode in ("planted", "clean")
        directory = joinpath(@__DIR__, "probes", "physics", mode)
        files = Guards.source_files(directory)
        reserved = Guards.reserved_vocabulary(files, directory, Guards.VOCABULARY)
        fingerprints = Guards.fingerprints(files, directory)
        if mode == "planted"
            @test reserved == Dict(
                "ext/HarmonicPowerModelsLineCableModelsExt/physics.jl | mu0" => 1,
                "ext/HarmonicPowerModelsLineCableModelsExt/physics.jl | eliminate!" => 1,
                "ext/HarmonicPowerModelsLineCableModelsExt/physics.jl | transposed!" => 1)
            @test fingerprints ==
                  Dict("ext/HarmonicPowerModelsLineCableModelsExt/physics.jl" => 3)
        else
            @test reserved == Dict{String, Int}()
            @test fingerprints == Dict{String, Int}()
        end
    end
end

@testset "file guards on the kit" begin
    for mode in ("planted", "clean")
        directory = joinpath(@__DIR__, "probes", "sources", mode)
        files = Guards.territory_files(directory)
        found = Guards.file_hygiene(files, directory)
        if mode == "planted"
            @test sum(values(found)) == 6
            @test found["test/linecablemodels/files.jl | homedir"] == 1
        else
            @test found == Dict{String, Int}()
        end
        directory = joinpath(@__DIR__, "probes", "narrow", mode)
        found = Guards.narrow_floats(directory)
        if mode == "planted"
            @test found == Dict(
                "ext/HarmonicPowerModelsLineCableModelsExt/values.jl | Float32" => 1,
                "test/linecablemodels/fixture.md | ComplexF32" => 1)
        else
            @test found == Dict{String, Int}()
        end
    end
end

@testset "public body clone probes" begin
    directory = normpath(joinpath(@__DIR__, "../../.."))
    tree = Guards.PackageTree(Ext, Module[], directory)
    methods = Guards.package_methods(tree)
    for mode in ("planted", "clean")
        directory = joinpath(@__DIR__, "probes", "clones", mode)
        files = Guards.source_files(directory)
        subject = Guards.PackageTree(Ext, Module[], directory)
        found = Guards.clones(files, methods, subject)
        if mode == "planted"
            @test found ==
                  Dict("ext/HarmonicPowerModelsLineCableModelsExt/clone.jl | copied ~ branch_data" =>
                1)
        else
            @test found == Dict{String, Int}()
        end
    end
end
