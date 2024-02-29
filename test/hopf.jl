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

    @testset "Industrial Example" begin
        # read-in data
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hopf.m")
        data = PowerModels.parse_file(path)

        # set the formulation
        form = PMs.IVRPowerModel

        # define the set of considered harmonics
        H=[1, 3, 5, 7, 9, 13]

        # build xfmr magnetization data
        B⁺  = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
        H⁺  = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
        Bᵗ  = vcat(reverse(-B⁺),0.0,B⁺)
        Hᵗ  = vcat(reverse(-H⁺),0.0,H⁺) 
        BH_powercore_h100_23 = Dierckx.Spline1D(Bᵗ, Hᵗ; k=3, bc="nearest")
        magn = Dict("Hᴱ"    => [1, 5], 
                    "Hᴵ"    => [1, 3, 5, 7, 9, 13],
                    "Emax"  => 1.1,
                    "IDH"   => [1.0, 0.06],
                    "pcs"   => [5, 5],
                    "xfmr"  => Dict(1 => Dict(  "l"     => 11.4,
                                                "A"     => 0.5,
                                                "N"     => 500,
                                                "BH"    => BH_powercore_h100_23,
                                                "Vbase" => 150000),
                                    2 => Dict(  "l"     => 8.0,
                                                "A"     => 0.2,
                                                "N"     => 300,
                                                "BH"    => BH_powercore_h100_23,
                                                "Vbase" => 36000),
                                    3 => Dict(  "l"     => 3.1,
                                                "A"     => 0.07,
                                                "N"     => 240,
                                                "BH"    => BH_powercore_h100_23,
                                                "Vbase" => 10000),
                                    4 => Dict(  "l"     => 1.0,
                                                "A"     => 0.001,
                                                "N"     => 240,
                                                "BH"    => BH_powercore_h100_23,
                                                "Vbase" => 690)
                                    )
                    )

        # solve HOPF problem w/o xfmr magnitization
        hdata_wo = HPM.replicate(data, H=H, bus_id=6)
        results_wo = HPM.solve_hopf(hdata_wo, PMs.IVRPowerModel, solver_nlp)

        # solve HOPF problem w. xfmr magnitization
        hdata_w = HPM.replicate(data, H=H, xfmr_magn=magn, bus_id=6)
        results_w = HPM.solve_hopf(hdata_wo, PMs.IVRPowerModel, solver_nlp)

        @testset "HOPF w/o Xfmr Magnitization" begin
            # solved to optimality
            @test results_wo["termination_status"] == LOCALLY_SOLVED
            # objective value
            @test results_wo["objective"] ≈ 0.0
        end

        @testset "HOPF w. Xfmr Magnitization" begin
            # solved to optimality
            @test results_w["termination_status"] == LOCALLY_SOLVED
            # objective value
            @test results_w["objective"] ≈ 0.0
        end
    end

end
