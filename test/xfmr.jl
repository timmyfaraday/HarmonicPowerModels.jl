
@testset "harmonic transformer model" begin 
    @testset "YN - yn - 0" begin
        path  = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_xfmr_YNyn0.m")
        data  = _PMs.parse_file(path)
        hdata = _HPM.replicate(data)

        run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
    end

    @testset "D - yn - 11" begin
        path  = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_xfmr_Dyn11.m")
        data  = _PMs.parse_file(path)
        hdata = _HPM.replicate(data)

        run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
    end

    @testset "Y - d - 11" begin
        path  = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_xfmr_Yd11.m")
        data  = _PMs.parse_file(path)
        hdata = _HPM.replicate(data)

        run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
    end

    @testset "YN - d - 11" begin
        path  = joinpath(_HPM.BASE_DIR,"test/data/matpower/case_xfmr_YNd11.m")
        data  = _PMs.parse_file(path)
        hdata = _HPM.replicate(data)

        run_hopf_iv(hdata, _PMs.IVRPowerModel, solver)
    end
end