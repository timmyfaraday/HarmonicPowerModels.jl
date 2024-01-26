################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

@testset "Harmonic Power Flow" begin

    @testset "Two-Bus Example" begin
        # Example considering a harmonic power flow for a two-bus example 
        # network taken from:
        # > Harmonic Optimal Power Flow with Transformer Excitation by F. Geth 
        #   and T. Van Acker, pg. 7, § IV.A.

        # read-in data
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/two_bus_example_hpf.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = PMs.IVRPowerModel

        # define the set of considered harmonics
        H=[1, 3]

        # build the harmonic data
        hdata = HPM.replicate(data, H=H)

        # power flow
        results_pf = PowerModels.solve_pf_iv(data, form, solver_nlp)
        @test results_pf["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_pf["objective"], 0; atol = 1e0)
        @test isapprox(results_pf["solution"]["bus"]["2"]["vr"],  0.96628746; atol = 1e-3)
        @test isapprox(results_pf["solution"]["bus"]["2"]["vi"], -0.02400000; atol = 1e-3)
        
        results_hpf = HPM.solve_hpf(hdata, form, solver_nlp)

        @test results_hpf["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_hpf["objective"], 0; atol = 1e0)
        @test isapprox(results_hpf["solution"]["nw"]["1"]["bus"]["2"]["vr"],  0.96628746; atol = 1e-3)
        @test isapprox(results_hpf["solution"]["nw"]["1"]["bus"]["2"]["vi"], -0.02400000; atol = 1e-3)
        @test isapprox(results_hpf["solution"]["nw"]["3"]["bus"]["2"]["vr"],  0.03391677; atol = 1e-3)
        @test isapprox(results_hpf["solution"]["nw"]["3"]["bus"]["2"]["vi"], -0.01660243; atol = 1e-3)

    end
end
