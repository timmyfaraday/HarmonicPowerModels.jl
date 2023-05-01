################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
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

        # build the harmonic data
        hdata = HPM.replicate(data)

        ihdmax = Dict("1" => 1.10, "3" => 0.05)
        for (nw,ntw) in hdata["nw"], (nb,bus) in ntw["bus"]
            bus["ihdmax"] = ihdmax[nw]
        end 
        # harmonic power flow
        results_harm = HPM.solve_hopf(hdata, form, solver)

        @test results_harm["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_harm["objective"], 0; atol = 1e0)
        @test isapprox(results_harm["solution"]["nw"]["1"]["bus"]["2"]["vr"], 0.966287; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["1"]["bus"]["2"]["vi"], -0.024; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["3"]["bus"]["2"]["vr"], -0.0160832; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["3"]["bus"]["2"]["vi"], -0.01660243; atol = 1e-3)

    end

    @testset "Industrial Example" begin
        # read-in data
        path = joinpath(HarmonicPowerModels.BASE_DIR,"test/data/matpower/industrial_network_hopf.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = PMs.IVRPowerModel

        # build the harmonic data
        hdata = HarmonicPowerModels.replicate(data)
        ihdmax = Dict("1" => 1.10, "3" => 0.05, "5" => 0.06, "7" => 0.05, "9" => 0.015, "13" => 0.03)
        for (nw,ntw) in hdata["nw"], (nb,bus) in ntw["bus"]
            bus["ihdmax"] = ihdmax[nw]
        end 
        results_harm = HarmonicPowerModels.solve_hopf(hdata, form, solver)

        @test results_harm["termination_status"] == LOCALLY_SOLVED
        @test isapprox(results_harm["objective"], 0; atol = 1e0)

        @test isapprox(results_harm["solution"]["nw"]["1"]["bus"]["2"]["vr"], 0.7951823447530401; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["1"]["bus"]["2"]["vi"], -0.568440001886; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["5"]["bus"]["8"]["vr"], -0.0241749759187737; atol = 1e-3)
        @test isapprox(results_harm["solution"]["nw"]["5"]["bus"]["8"]["vi"], 0.028023283097463714; atol = 1e-3)
    end

end
