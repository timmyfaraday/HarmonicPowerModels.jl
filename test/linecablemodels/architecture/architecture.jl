include("guards.jl")
using Aqua: Aqua

@testset "architecture" begin
    directory = normpath(joinpath(@__DIR__, "../../.."))
    tree = Guards.PackageTree(Ext, Module[], directory)
    files = Guards.source_files(directory)
    territory = Guards.territory_files(directory)
    methods = Guards.package_methods(tree)
    homes = Guards.home_directories(tree)
    names = Guards.TOML.parsefile(joinpath(@__DIR__, "names.toml"))
    checks = (
        ("ownership", "adapted", () -> Guards.ownership(methods, tree)),
        ("placement", "kept", () -> Guards.placement(methods, tree, homes)),
        ("names", "adapted", () -> Guards.shared_names(tree)),
        ("explicit imports", "kept", () -> Guards.explicit_imports(tree, Module[])),
        ("validate returns its subject", "kept",
            () -> Guards.validate_returns(files, directory)),
        ("reserved verbs", "kept", () -> Guards.reserved_verbs(files, directory)),
        ("storage vocabulary", "adapted",
            () -> Guards.storage_vocabulary(files, directory)),
        ("repeated access names", "kept",
            () -> Guards.repeated_access_names(files, directory)),
        ("symbol switches and probes", "adapted",
            () -> Guards.symbol_switches(files, directory)),
        ("decisions by absence", "kept", () -> Guards.absences(files, directory)),
        ("file hygiene", "kept", () -> Guards.file_hygiene(territory, directory)),
        ("narrow floating types", "kept", () -> Guards.narrow_floats(directory)),
        ("reserved vocabulary", "adapted",
            () -> Guards.reserved_vocabulary(files, directory, Guards.VOCABULARY)),
        ("literal fingerprints", "kept", () -> Guards.fingerprints(files, directory)),
        ("clones", "adapted", () -> Guards.clones(files, methods, tree)),
        ("register", "kept", () -> Guards.register_differences(tree, names)),
        ("local functions", "kept", () -> Guards.local_functions(files, directory)),
        ("forwarders", "kept", () -> Guards.forwarders(tree)),
        ("boxed captures", "kept", () -> Guards.boxed_captures(methods, tree)),
        ("typed values", "kept", () -> Guards.typed_values(tree)),
        ("Val tags", "adapted", () -> Guards.val_tags(methods, tree, names)),
        ("naming rule", "kept", () -> Guards.naming_rule(methods, tree, names)))
    for (name, status, check) in checks
        found = check()
        println(name, " (", status, "): ", sum(values(found); init = 0), " ", found)
        isempty(found) || error("Architecture guard blocked: $name")
        @testset "$name" begin
            @test found == Dict{String, Int}()
        end
    end
    @testset "Aqua piracy" begin
        Aqua.test_piracies(Ext;
            treat_as_own = [HarmonicPowerModels.ExactPiModel,
                HarmonicPowerModels.ParameterBox, LineCableModels.AbstractFormulation])
    end
end

include("probes.jl")
