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

        # build the harmonic data
        hdata = HPM.replicate(data)

        # harmonic power flow
        results_fund = PowerModels.solve_pf_iv(data, form, solver)
        @test results_fund["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_fund["objective"], 0; atol = 1e0)
        @test isapprox(results_fund["solution"]["bus"]["2"]["vr"],  0.966287; atol = 1e-3)
        @test isapprox(results_fund["solution"]["bus"]["2"]["vi"], -0.024; atol = 1e-3)
        
        
        results_harm = HPM.solve_hpf(hdata, form, solver)

        @test results_harm["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_harm["objective"], 0; atol = 1e0)
        @test isapprox(results_harm["solution"]["nw"]["1"]["bus"]["2"]["vr"], 0.966287; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["1"]["bus"]["2"]["vi"], -0.024; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["3"]["bus"]["2"]["vr"], -0.0160832; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["3"]["bus"]["2"]["vi"], -0.01660243; atol = 1e-3)

    end
end
