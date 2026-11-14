@testset "Manual session" begin
    page = joinpath(@__DIR__, "../../docs/src/linecablemodels.md")
    session = Module(:ManualSession)
    Base.include_string(session, "using Test: @test", page)
    blocks = collect(eachmatch(r"(?ms)^```julia\n(.*?)^```", read(page, String)))
    @test !isempty(blocks)
    for block in blocks
        code = block.captures[1]
        startswith(code, "# manual-skip") && continue
        @testset "Block $(block.offset)" begin
            Base.include_string(session, replace(code, r"@assert\b" => "@test"), page)
        end
    end
end
