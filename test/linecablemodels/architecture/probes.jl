# LCM's negative controls keep their per-name counts after the path changes.
function probe(directory, name; extension = nothing)
    root_directory = joinpath(@__DIR__, "probes", directory)
    root = Base.include(Main, joinpath(root_directory, Guards.HOME, name * ".jl"))
    extensions = extension === nothing ? Module[] :
                 [Base.include(Main, joinpath(root_directory, "ext", extension, extension *
                                                                                ".jl"))]
    tree = Base.invokelatest(Guards.PackageTree, root, extensions, root_directory)
    return tree, extensions
end

@testset "negative controls" begin
    planted,
    planted_extensions = probe("imports/planted", "ImportsProbePlanted";
        extension = "ImportsProbePlantedExt")
    @test Base.invokelatest(Guards.explicit_imports, planted, planted_extensions) ==
          Dict{String, Int}(
        "no_implicit_imports | Consumer | Exporter" => 1,
        "no_implicit_imports | Consumer | implicit_value" => 1,
        "no_stale_explicit_imports | Consumer | public_value" => 1,
        "all_explicit_imports_via_owners | Consumer | relayed_value" => 1,
        "all_explicit_imports_are_public | Consumer | relayed_value" => 1,
        "all_explicit_imports_are_public | Consumer | private_import" => 1,
        "all_qualified_accesses_via_owners | Consumer | relayed_value" => 1,
        "all_qualified_accesses_are_public | Consumer | relayed_value" => 1,
        "all_qualified_accesses_are_public | Consumer | private_value" => 1,
        "all_qualified_accesses_are_public | ImportsProbePlantedExt | private_value" => 1,
        "no_self_qualified_accesses | Consumer | implicit" => 1)

    clean,
    clean_extensions = probe("imports/clean", "ImportsProbeClean";
        extension = "ImportsProbeCleanExt")
    @test Base.invokelatest(Guards.explicit_imports, clean, clean_extensions) ==
          Dict{String, Int}()

    directory = joinpath(@__DIR__, "probes", "sources", "planted")
    files = Guards.source_files(directory)
    @test Guards.validate_returns(files, directory) ==
          Dict{String, Int}("ext/HarmonicPowerModelsLineCableModelsExt/validate.jl" => 6)
    @test Guards.reserved_verbs(files, directory) == Dict{String, Int}(
        "ext/HarmonicPowerModelsLineCableModelsExt/validate.jl | __validate_shape" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/validate.jl | _check_input" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/validate.jl | ensure_ready" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/validate.jl | require_kind" => 1)
    @test Guards.storage_vocabulary(files, directory) == Dict{String, Int}(
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | AbstractStorage" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | AbstractBuffers" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | FormulaBuffers" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | KernelCache" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | cache" => 2,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | kernel_workspace" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | initialize_buffers" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | sample_storage" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | scratch" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | storage" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | trace_buffers" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/storage.jl | work" => 2)
    @test Guards.repeated_access_names(files, directory) ==
          Dict{String, Int}("ext/HarmonicPowerModelsLineCableModelsExt/repeats.jl" => 3)
    @test Guards.symbol_switches(files, directory) == Dict{String, Int}(
        "ext/HarmonicPowerModelsLineCableModelsExt/strings.jl | symbol" => 2,
        "ext/HarmonicPowerModelsLineCableModelsExt/switches.jl | applicable" => 2,
        "ext/HarmonicPowerModelsLineCableModelsExt/switches.jl | eval" => 2,
        "ext/HarmonicPowerModelsLineCableModelsExt/switches.jl | symbol" => 9,
        "ext/HarmonicPowerModelsLineCableModelsExt/switches.jl | type test" => 9)
    @test Guards.absences(files, directory) == Dict{String, Int}(
        "ext/HarmonicPowerModelsLineCableModelsExt/absence.jl | argument" => 9,
        "ext/HarmonicPowerModelsLineCableModelsExt/absence.jl | comparison" => 9,
        "ext/HarmonicPowerModelsLineCableModelsExt/absence.jl | ismissing" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/absence.jl | isnothing" => 4,
        "ext/HarmonicPowerModelsLineCableModelsExt/absence.jl | parameter" => 7,
        "ext/HarmonicPowerModelsLineCableModelsExt/absence.jl | something" => 1)
    @test Guards.file_hygiene(files, directory) == Dict{String, Int}(
        "ext/HarmonicPowerModelsLineCableModelsExt/extension_hygiene.jl | homedir" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/hygiene.jl | homedir" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/hygiene.jl | mktempdir" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/hygiene.jl | tempdir" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/hygiene.jl | tempname" => 1)

    directory = joinpath(@__DIR__, "probes", "sources", "clean")
    files = Guards.source_files(directory)
    @test Guards.validate_returns(files, directory) == Dict{String, Int}()
    @test Guards.reserved_verbs(files, directory) == Dict{String, Int}()
    @test Guards.storage_vocabulary(files, directory) == Dict{String, Int}()
    @test Guards.repeated_access_names(files, directory) == Dict{String, Int}()
    @test Guards.symbol_switches(files, directory) == Dict{String, Int}()
    @test Guards.absences(files, directory) == Dict{String, Int}()
    @test Guards.file_hygiene(files, directory) == Dict{String, Int}()

    planted, _ = probe("values/planted", "ValuesProbePlanted")
    planted_methods = Base.invokelatest(Guards.package_methods, planted)
    @test Base.invokelatest(Guards.boxed_captures, planted_methods, planted) ==
          Dict{String, Int}(
        "ValuesProbePlanted | boxed | ext/HarmonicPowerModelsLineCableModelsExt/ValuesProbePlanted.jl" =>
            1,
        "ValuesProbePlanted | boxed_keyword | ext/HarmonicPowerModelsLineCableModelsExt/ValuesProbePlanted.jl" =>
            1)
    @test Base.invokelatest(Guards.typed_values, planted) == Dict{String, Int}(
        "ValuesProbePlanted | Blend | x" => 1,
        "ValuesProbePlanted | Either | x" => 1,
        "ValuesProbePlanted | Loose | constructor | ext/HarmonicPowerModelsLineCableModelsExt/ValuesProbePlanted.jl" =>
            1,
        "ValuesProbePlanted | Loose | x" => 1,
        "ValuesProbePlanted | Nested | x" => 1,
        "ValuesProbePlanted | Nested | y" => 1,
        "ValuesProbePlanted | Open | constructor | ext/HarmonicPowerModelsLineCableModelsExt/ValuesProbePlanted.jl" =>
            1,
        "ValuesProbePlanted | Open | x" => 1,
        "ValuesProbePlanted | Plain | constructor | ext/HarmonicPowerModelsLineCableModelsExt/ValuesProbePlanted.jl" =>
            3,
        "ValuesProbePlanted | Wide | constructor | ext/HarmonicPowerModelsLineCableModelsExt/ValuesProbePlanted.jl" =>
            1,
        "ValuesProbePlanted | Wide | x" => 1)

    clean, _ = probe("values/clean", "ValuesProbeClean")
    clean_methods = Base.invokelatest(Guards.package_methods, clean)
    @test Base.invokelatest(Guards.boxed_captures, clean_methods, clean) ==
          Dict{String, Int}()
    @test Base.invokelatest(Guards.typed_values, clean) == Dict{String, Int}()

    planted, _ = probe("tags/planted", "TagsProbePlanted")
    planted_methods = Base.invokelatest(Guards.package_methods, planted)
    planted_names = Guards.TOML.parsefile(joinpath(planted.directory, "names.toml"))
    @test Base.invokelatest(Guards.val_tags, planted_methods, planted, planted_names) ==
          Dict{String, Int}(
        "Laws | choose | Val(x)" => 1,
        "Laws | gone | stale lift" => 1,
        "Laws | Knob | Val first" => 1,
        "Laws | law | Val first" => 1,
        "Laws | pick | Val(x)" => 1,
        "Uses | pick | Val(x)" => 1)

    clean, _ = probe("tags/clean", "TagsProbeClean")
    clean_methods = Base.invokelatest(Guards.package_methods, clean)
    clean_names = Guards.TOML.parsefile(joinpath(clean.directory, "names.toml"))
    @test Base.invokelatest(Guards.val_tags, clean_methods, clean, clean_names) ==
          Dict{String, Int}()

    planted,
    _ = probe("register/planted", "RegisterProbePlanted";
        extension = "RegisterProbePlantedExt")
    planted_names = Guards.TOML.parsefile(joinpath(planted.directory, "names.toml"))
    planted_files = Guards.source_files(planted.directory)
    @test Base.invokelatest(Guards.register_differences, planted, planted_names) ==
          Dict{String, Int}(
        "Alpha | @twice | unregistered" => 1,
        "Alpha | Shapes | unregistered" => 1,
        "Alpha | _ABSENT | unregistered" => 1,
        "Alpha | _Absent | unregistered" => 1,
        "Alpha | _half | unregistered" => 1,
        "Alpha | area | twice" => 1,
        "Alpha | depth | private, registered public" => 1,
        "Alpha | gone | stale" => 1,
        "Alpha | half | unregistered" => 1,
        "Alpha | outer | meaning" => 1,
        "Alpha | radius | public, registered private" => 1,
        "Alpha | total | unregistered" => 1,
        "Commons | _pad | unregistered" => 1,
        "Missing | stale module" => 1,
        "RegisterProbePlantedExt | _mesh | unregistered" => 1)
    @test Base.invokelatest(Guards.local_functions, planted_files, planted.directory) ==
          Dict{String, Int}(
        "ext/HarmonicPowerModelsLineCableModelsExt/alpha/Alpha.jl | double" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/alpha/Alpha.jl | inner" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/alpha/Alpha.jl | shift" => 1,
        "ext/HarmonicPowerModelsLineCableModelsExt/alpha/Alpha.jl | shifted" => 1)
    @test Base.invokelatest(Guards.forwarders, planted) == Dict{String, Int}(
        "Alpha | alias" => 1,
        "RegisterProbePlantedExt | draw" => 1)

    clean,
    _ = probe("register/clean", "RegisterProbeClean";
        extension = "RegisterProbeCleanExt")
    clean_names = Guards.TOML.parsefile(joinpath(clean.directory, "names.toml"))
    clean_files = Guards.source_files(clean.directory)
    @test Base.invokelatest(Guards.register_differences, clean, clean_names) ==
          Dict{String, Int}()
    @test Base.invokelatest(Guards.local_functions, clean_files, clean.directory) ==
          Dict{String, Int}()
    @test Base.invokelatest(Guards.forwarders, clean) == Dict{String, Int}()

    include("port_probes.jl")
end
