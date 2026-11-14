@testset "Manual case records" begin
    page = joinpath(@__DIR__, "../../../docs/src/linecablemodels_studies.md")
    blocks = collect(eachmatch(r"(?ms)^```toml\n(.*?)^```", read(page, String)))
    @test !isempty(blocks)
    for block in blocks
        code = block.captures[1]
        parsed = TOML.parse(code)
        @test !isempty(parsed)
        if haskey(parsed, "base") ||
           all(haskey(parsed, key) for key in ("model", "network", "box", "box_names"))
            mktempdir() do directory
                file = joinpath(directory, "manual.toml")
                write(file, code)
                s = study(file)
                @test s.model isa ExactPiModel
                @test Set(uncertain_branches(s.model)) == Set(keys(s.details))
                @test all(isfinite(branch_data(s.model, b, nominal_point(s.box), 5,
                              :negative).br_r) for b in uncertain_branches(s.model))
            end
        end
    end
end
