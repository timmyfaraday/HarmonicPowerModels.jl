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
     
        @test sol_pf["branch"]["1"]["cr_fr"]  ≈ sol_hpf["branch"]["1"]["cbr_fr"]
        @test sol_pf["branch"]["1"]["ci_fr"]  ≈ sol_hpf["branch"]["1"]["cbi_fr"]
        @test sol_pf["branch"]["1"]["csr_fr"] ≈ sol_hpf["branch"]["1"]["cbsr_fr"]
        @test sol_pf["branch"]["1"]["csi_fr"] ≈ sol_hpf["branch"]["1"]["cbsi_fr"]
        @test sol_pf["branch"]["1"]["cr_to"]  ≈ sol_hpf["branch"]["1"]["cbr_to"]
        @test sol_pf["branch"]["1"]["ci_to"]  ≈ sol_hpf["branch"]["1"]["cbi_to"]

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

    @testset "Two bus system, single xfmr (YNyn0), w/wo no load active losses" begin
        # Example comparing a two bus system with a single xfmr (YNyn0) with and
        # without no load losses, for the fundamental harmonic
        path_wo = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_xfmr_YNyn0_hpf.m")
        data_wo = PMs.parse_file(path_wo)
        path_w  = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_xfmr_YNyn0_w_P0_hpf.m")
        data_w  = PMs.parse_file(path_w)

        # define the set of considered harmonics
        H = [1]

        # build the harmonic data
        hdata_wo = HPM.build_hdata_from_matpower_file(data_wo, H=H)
        hdata_w  = HPM.build_hdata_from_matpower_file(data_w, H=H)

        # power flow
        results_hpf_wo = HPM.solve_hpf(hdata_wo, HPM.HarmonicPowerModel, solver_nlp)
        results_hpf_w  = HPM.solve_hpf(hdata_w, HPM.HarmonicPowerModel, solver_nlp)

        # feasibility tests 
        @test results_hpf_wo["termination_status"] == LOCALLY_SOLVED
        @test results_hpf_w["termination_status"]  == LOCALLY_SOLVED

        # equality tests
        sol_hpf_wo = results_hpf_wo["solution"]["nw"]["1"]
        sol_hpf_w  = results_hpf_w["solution"]["nw"]["1"]

        pg_wo   =   sol_hpf_wo["bus"]["1"]["vbr"] * sol_hpf_wo["gen"]["1"]["cgr"] 
                  + sol_hpf_wo["bus"]["1"]["vbi"] * sol_hpf_wo["gen"]["1"]["cgi"]
        qg_wo   =   sol_hpf_wo["bus"]["1"]["vbi"] * sol_hpf_wo["gen"]["1"]["cgr"] 
                  - sol_hpf_wo["bus"]["1"]["vbr"] * sol_hpf_wo["gen"]["1"]["cgi"]

        pg_w    =   sol_hpf_w["bus"]["1"]["vbr"] * sol_hpf_w["gen"]["1"]["cgr"] 
                  + sol_hpf_w["bus"]["1"]["vbi"] * sol_hpf_w["gen"]["1"]["cgi"]
        qg_w    =   sol_hpf_w["bus"]["1"]["vbi"] * sol_hpf_w["gen"]["1"]["cgr"] 
                  - sol_hpf_w["bus"]["1"]["vbr"] * sol_hpf_w["gen"]["1"]["cgi"]

        pl_w    = data_w["load"]["1"]["pd"]
        p0_w    = (sol_hpf_w["xfmr"]["1"]["exr"]^2 + sol_hpf_w["xfmr"]["1"]["exi"]^2) * data_w["xfmr"]["1"]["gsh"]
        p1_w    = (sol_hpf_w["xfmr"]["1"]["cxr_1"]^2 + sol_hpf_w["xfmr"]["1"]["cxi_1"]^2) * data_w["xfmr"]["1"]["r1"]
        p2_w    = (sol_hpf_w["xfmr"]["1"]["cxr_2"]^2 + sol_hpf_w["xfmr"]["1"]["cxi_2"]^2) * data_w["xfmr"]["1"]["r2"]

        pl_wo   = data_wo["load"]["1"]["pd"]
        p1_wo   = (sol_hpf_wo["xfmr"]["1"]["cxr_1"]^2 + sol_hpf_wo["xfmr"]["1"]["cxi_1"]^2) * data_wo["xfmr"]["1"]["r1"]
        p2_wo   = (sol_hpf_wo["xfmr"]["1"]["cxr_2"]^2 + sol_hpf_wo["xfmr"]["1"]["cxi_2"]^2) * data_wo["xfmr"]["1"]["r2"]

        @test pg_wo < pg_w
        @test qg_wo ≈ qg_w

        @test pg_w  ≈ pl_w + p0_w + p1_w + p2_w
        @test pg_wo ≈ pl_wo + p1_wo + p2_wo
    end

    @testset "Two bus system, single xfmr, YNyn0 vs YNyn3" begin
        # Example comparing a two bus system with a single xfmr with vector 
        # group YNyn0 and YNyn3, active and reactive power consumption should 
        # remain unchanged.
        path_0 = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_xfmr_YNyn0_hpf.m")
        data_0 = PMs.parse_file(path_0)
        path_3 = joinpath(HPM.BASE_DIR, "test/data/matpower/xfmr/two_bus_xfmr_YNyn3_hpf.m")
        data_3 = PMs.parse_file(path_3)

        # define the set of considered harmonics
        H = [1,2,3]

        # build the harmonic data
        hdata_0 = HPM.build_hdata_from_matpower_file(data_0, H=H)
        hdata_3 = HPM.build_hdata_from_matpower_file(data_3, H=H)

        # power flow
        results_hpf_0 = HPM.solve_hpf(hdata_0, HPM.HarmonicPowerModel, solver_nlp)
        results_hpf_3 = HPM.solve_hpf(hdata_3, HPM.HarmonicPowerModel, solver_nlp)

        # feasibility tests 
        @test results_hpf_0["termination_status"] == LOCALLY_SOLVED
        @test results_hpf_3["termination_status"]  == LOCALLY_SOLVED

        # equality tests
        sol_hpf_0 = results_hpf_0["solution"]["nw"]["1"]
        sol_hpf_3 = results_hpf_3["solution"]["nw"]["1"]

        pg_0    =   sol_hpf_0["bus"]["1"]["vbr"] * sol_hpf_0["gen"]["1"]["cgr"] 
                  + sol_hpf_0["bus"]["1"]["vbi"] * sol_hpf_0["gen"]["1"]["cgi"]
        qg_0    =   sol_hpf_0["bus"]["1"]["vbi"] * sol_hpf_0["gen"]["1"]["cgr"] 
                  - sol_hpf_0["bus"]["1"]["vbr"] * sol_hpf_0["gen"]["1"]["cgi"]

        pg_3    =   sol_hpf_3["bus"]["1"]["vbr"] * sol_hpf_3["gen"]["1"]["cgr"] 
                  + sol_hpf_3["bus"]["1"]["vbi"] * sol_hpf_3["gen"]["1"]["cgi"]
        qg_3    =   sol_hpf_3["bus"]["1"]["vbi"] * sol_hpf_3["gen"]["1"]["cgr"] 
                  - sol_hpf_3["bus"]["1"]["vbr"] * sol_hpf_3["gen"]["1"]["cgi"]

        vba_0   = atand(sol_hpf_0["bus"]["2"]["vbi"],sol_hpf_0["bus"]["2"]["vbr"])
        vba_3   = atand(sol_hpf_3["bus"]["2"]["vbi"],sol_hpf_3["bus"]["2"]["vbr"])

        @test pg_0 ≈ pg_3
        @test qg_0 ≈ qg_3

        # The secondary winding is vg * 30° lagging with reference to the 
        # primary, positive sense of rotation is counter-clockwise, i.e., 
        # 1 o'clock is one hour behind 12 o'clock.
        @test vba_0 ≈ vba_3 + 3*30
    end

end