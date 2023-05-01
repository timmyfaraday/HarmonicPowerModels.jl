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

        ihdmax = Dict("1" => 1.10, "3" => 0.05)
        for (nw,ntw) in hdata["nw"], (nb,bus) in ntw["bus"]
            bus["ihdmax"] = ihdmax[nw]
        end 


        # harmonic power flow
        results_fund = PowerModels.solve_pf_iv(data, form, solver)
        # println("Results for the fundamental power flow:")
        # print_summary(results_fund["solution"]) 
        
        
        results_harm = HPM.solve_hpf(hdata, form, solver)

        # tests 
        # println("Results for the harmonic power flow")
        # println("Fundamental harmonic:")
        # print_summary(results_harm["solution"]["nw"]["1"])
        # println("Third harmonic:")
        # print_summary(results_harm["solution"]["nw"]["3"])
    end

end
