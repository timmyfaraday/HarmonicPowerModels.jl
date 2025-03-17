################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

@testset "Harmonic Power Flow - xfmr" begin

    @testset "Two bus system, single branch - IVR vs fund. HarmonicPowerModel" begin
        # Example comparing a two bus system with a single branch for 
        # PMs.IVRPowerModel and HarmonicPowerModel, for the fundamental harmonic
        path = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_branch_hpf.m")
        data = PMs.parse_file(path)

        # define the set of considered harmonics
        H = [1]

        # build the harmonic data
        hdata = HPM.build_hdata_from_matpower_file(data, H=H)

        # power flow
        results_pf  = PMs.solve_pf_iv(data, PMs.IVRPowerModel, solver_nlp)
        results_hpf = HPM.solve_hpf(hdata, HPM.HarmonicPowerModel, solver_nlp)

        # feasibility tests 
        @test results_pf["termination_status"] == LOCALLY_SOLVED
        @test results_hpf["termination_status"] == LOCALLY_SOLVED

        # equality tests
        sol_pf  = results_pf["solution"]
        sol_hpf = results_hpf["solution"]["nw"]["1"]

        @test sol_pf["bus"]["1"]["vr"] ≈ sol_hpf["bus"]["1"]["vbr"]
        @test sol_pf["bus"]["1"]["vi"] ≈ sol_hpf["bus"]["1"]["vbi"]
        @test sol_pf["bus"]["2"]["vr"] ≈ sol_hpf["bus"]["2"]["vbr"]
        @test sol_pf["bus"]["2"]["vi"] ≈ sol_hpf["bus"]["2"]["vbi"]
     
        @test sol_pf["branch"]["1"]["cr_fr"]    ≈ sol_hpf["branch"]["1"]["cbr_fr"]
        @test sol_pf["branch"]["1"]["ci_fr"]    ≈ sol_hpf["branch"]["1"]["cbi_fr"]
        @test sol_pf["branch"]["1"]["csr_fr"]   ≈ sol_hpf["branch"]["1"]["cbsr_fr"]
        @test sol_pf["branch"]["1"]["csi_fr"]   ≈ sol_hpf["branch"]["1"]["cbsi_fr"]
        @test sol_pf["branch"]["1"]["cr_to"]    ≈ sol_hpf["branch"]["1"]["cbr_to"]
        @test sol_pf["branch"]["1"]["ci_to"]    ≈ sol_hpf["branch"]["1"]["cbi_to"]

        @test sol_pf["gen"]["1"]["crg"] ≈ sol_hpf["gen"]["1"]["cgr"]
        @test sol_pf["gen"]["1"]["cig"] ≈ sol_hpf["gen"]["1"]["cgi"]
    end

    @testset "Two bus system, single branch vs single xfmr (YNyn0)" begin
        # Example comparing a two bus system with a single branch with one with
        # a single xfmr (YNyn0), for the fundamental harmonic
        path_b  = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_branch_hpf.m")
        data_b  = PMs.parse_file(path_b)
        path_x  = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_xfmr_YNyn0_hpf.m")
        data_x  = PMs.parse_file(path_x)

        # define the set of considered harmonics
        H = [1]

        # build the harmonic data
        hdata_b = HPM.build_hdata_from_matpower_file(data_b, H=H)
        hdata_x = HPM.build_hdata_from_matpower_file(data_x, H=H)

        # power flow
        results_hpf_b   = HPM.solve_hpf(hdata_b, HPM.HarmonicPowerModel, solver_nlp)
        results_hpf_x   = HPM.solve_hpf(hdata_x, HPM.HarmonicPowerModel, solver_nlp)

        # feasibility tests 
        @test results_hpf_b["termination_status"] == LOCALLY_SOLVED
        @test results_hpf_x["termination_status"] == LOCALLY_SOLVED

        # equality tests
        sol_hpf_b = results_hpf_b["solution"]["nw"]["1"]
        sol_hpf_x = results_hpf_x["solution"]["nw"]["1"]

        @test sol_hpf_b["bus"]["1"]["vbr"] ≈ sol_hpf_x["bus"]["1"]["vbr"]
        @test sol_hpf_b["bus"]["1"]["vbi"] ≈ sol_hpf_x["bus"]["1"]["vbi"]
        @test sol_hpf_b["bus"]["2"]["vbr"] ≈ sol_hpf_x["bus"]["2"]["vbr"]
        @test sol_hpf_b["bus"]["2"]["vbi"] ≈ sol_hpf_x["bus"]["2"]["vbi"]
     
        @test sol_hpf_b["branch"]["1"]["cbr_fr"] ≈ sol_hpf_x["xfmr"]["1"]["cxr_1"]
        @test sol_hpf_b["branch"]["1"]["cbi_fr"] ≈ sol_hpf_x["xfmr"]["1"]["cxi_1"]
        @test sol_hpf_b["branch"]["1"]["cbr_to"] ≈ sol_hpf_x["xfmr"]["1"]["cxr_2"]
        @test sol_hpf_b["branch"]["1"]["cbi_to"] ≈ sol_hpf_x["xfmr"]["1"]["cxi_2"]
        
        @test sol_hpf_b["gen"]["1"]["cgr"] ≈ sol_hpf_x["gen"]["1"]["cgr"]
        @test sol_hpf_b["gen"]["1"]["cgi"] ≈ sol_hpf_x["gen"]["1"]["cgi"]
    end

end