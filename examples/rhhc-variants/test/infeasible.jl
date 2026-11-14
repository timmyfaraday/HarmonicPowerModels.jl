@testset "Infeasible nominal reports" begin
    directory = mktempdir()
    label = basename(directory)
    paths = [joinpath(directory, label * suffix * ".toml")
             for suffix in ("-infeasible", "-feasible")]
    reports = [joinpath(directory, splitext(basename(path))[1] * ".txt")
               for path in paths]
    try
        config = TOML.parsefile(joinpath(STUDIES, "cases", "ieee39-physical.toml"))
        config["network"]["file"] = joinpath(STUDIES, "networks", "ieee39.json")
        open(last(paths), "w") do io
            TOML.print(io, config)
        end
        # Restore the fixed Cardinal ceiling that produced the infeasible nominal.
        config["branches"]["31"]["rating"] = "Cardinal"
        open(first(paths), "w") do io
            TOML.print(io, config)
        end
        s = study(first(paths))
        io = IOBuffer()
        result = runstudy(s, io)
        text = String(take!(io))
        @test string(result.nominal["termination_status"]) == "INFEASIBLE"
        @test string(result.robust["termination_status"]) == "INFEASIBLE"
        @test !result.nominal["signed_off"]
        @test !result.robust["signed_off"]
        @test isempty(result.nominal["solution"])
        @test isempty(result.robust["solution"])
        @test occursin(r"nominal solution: Dict[^\n]*\(\)", text)
        @test occursin(r"solution: Dict[^\n]*\(\)", text)
        @test occursin("resonances at the nominal point:", text)
        @test occursin("Posterior check not evaluated: the solution is empty.", text)
        @test occursin("LineCableModels solves", text)
        comparison(io, reports, paths...)
        text = String(take!(io))
        @test occursin("INFEASIBLE", text)
        @test occursin("OPTIMAL", text)
        @test occursin("not evaluated: empty solution (INFEASIBLE)", text)
        infeasible = splitext(basename(first(paths)))[1]
        feasible = splitext(basename(last(paths)))[1]
        @test occursin(infeasible * ": nominal solution = Dict", text)
        @test occursin(infeasible * ": solution = Dict", text)
        @test count(line -> startswith(line, feasible * " | " * infeasible * " | "),
            split(text, '\n')) == 2
        @test occursin("nominal solution: Dict", read(first(reports), String))
        @test occursin("a posteriori (1024 random points)", read(last(reports), String))
    finally
        for path in reports
            isfile(path) && rm(path)
        end
        rm(directory; recursive = true)
    end
end
