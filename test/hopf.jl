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

        # harmonic power flow
        results_harm = HPM.solve_hopf(hdata, form, solver)

    end

end
