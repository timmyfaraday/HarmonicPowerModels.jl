
@testset "harmonic transformer model" begin 
    @testset "YN - yn - 0" begin
        data = _PMs.parse_file("data/case_xfmr_YNyn0.m")
        data = replicate(data,harmonics = [1,3,5,9])
    end
end