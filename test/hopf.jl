################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
################################################################################

"""
Example considering harmonic optimal power flow for a two-bus example.
    This example extends the PF example. 
"""

@testset "Harmonic Optimal Power Flow" begin

    @testset "Two-Bus Example" begin
        # read-in data
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/two_bus_example_hpf.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = PMs.IVRPowerModel

        # define the set of considered harmonics
        H=[1, 3]

        # build the harmonic data
        hdata = HPM.replicate(data, H=H)

        # solve HOPF problem
        results_hopf = HPM.solve_hopf(hdata, form, solver_nlp)

        @test results_hopf["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_hopf["objective"], 0; atol = 1e0)
        @test isapprox(results_hopf["solution"]["nw"]["1"]["bus"]["2"]["vr"],  0.96628746; atol = 1e-3)
        @test isapprox(results_hopf["solution"]["nw"]["1"]["bus"]["2"]["vi"], -0.02400000; atol = 1e-3)
        @test isapprox(results_hopf["solution"]["nw"]["3"]["bus"]["2"]["vr"],  0.03391677; atol = 1e-3)
        @test isapprox(results_hopf["solution"]["nw"]["3"]["bus"]["2"]["vi"], -0.01660243; atol = 1e-3)
    end

    # @testset "Industrial Example" begin
    #     # read-in data
    #     path = joinpath(HarmonicPowerModels.BASE_DIR,"test/data/matpower/industrial_network_hopf.m")
    #     data = PMs.parse_file(path)

    #     # set the formulation
    #     form = PMs.IVRPowerModel

    #     # define the set of considered harmonics
    #     H=[1, 3, 5, 7, 9, 13]

    #     # build the harmonic data
    #     hdata = HarmonicPowerModels.replicate(data, H=H)

    #     # solve HOPF problem
    #     results_hopf = HarmonicPowerModels.solve_hopf(hdata, form, solver_nlp)

    #     @test results_hopf["termination_status"] == LOCALLY_SOLVED
    #     @test isapprox(results_hopf["objective"], 0; atol = 1e0)
    #     @test isapprox(results_hopf["solution"]["nw"]["1"]["bus"]["2"]["vr"], 0.7951823447530401; atol = 1e-3)
    #     @test isapprox(results_hopf["solution"]["nw"]["1"]["bus"]["2"]["vi"], -0.568440001886; atol = 1e-3)
    #     # @test isapprox(results_harm["solution"]["nw"]["5"]["bus"]["8"]["vr"], -0.03252192960993166; atol = 1e-3)
    #     # @test isapprox(results_harm["solution"]["nw"]["5"]["bus"]["8"]["vi"], 0.022870556302375344; atol = 1e-3)
    # end

end
