################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

@testset "Transformer Model" begin
    @testset "Yy0" begin 
        # read-in data
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/xfmr/xfmr_yy0.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = PMs.IVRPowerModel

        # define the set of considered harmonics
        H⁺ = [1, 7, 13]
        H⁻ = [5]
        H⁰ = [3, 9]

        # pos., neg. and zero-sequence harmonics -- infeasible
        hdata = HPM.replicate(data, H=sort(union(H⁺, H⁻, H⁰)))
        results = HPM.solve_hpf(hdata, form, solver_nlp)

        @test results["termination_status"] == LOCALLY_INFEASIBLE

        # pos. and neg.-sequence harmonics -- feasible
        hdata = HPM.replicate(data, H=sort(union(H⁺, H⁻)))
        results = HPM.solve_hpf(hdata, form, solver_nlp)

        @test results["termination_status"] == LOCALLY_SOLVED

        @testset "Transformer Phase Shifts" begin
            # The phase shift for the harmonic voltage over a transformer, i.e.,
            # between the excitation Eₓₕ and internal voltage of the second 
            # winding V₂ₕ in a balanced (optimal) power flow should relate as 
            # follows:
            # - 'positive sequence' harmonics, i.e., h = 1 + 3n = 1, 4, 7...
            #       θⱽ - θᵉ ≈ vg * π/6
            # - 'negative sequence' harmonics, i.e., h = 2 + 3n = 2, 5, 8...
            #       θⱽ - θᵉ ≈ -vg * π/6
            # - 'zero sequence' harmonics, i.e., h = 3 + 3n = 3, 6, 9...
            #       θⱽ ≈ θᵉ
            for nh ∈ sort(union(H⁺, H⁻)), (nx,xfmr) in data["xfmr"]
                vg = parse(Int, xfmr["vg"][3:end])
                
                sol_xfmr = results["solution"]["nw"]["$nh"]["xfmr"][nx]

                θe = angle(sol_xfmr["ert"] + im * sol_xfmr["eit"])
                θv = angle(sol_xfmr["vrt_to"] + im * sol_xfmr["vit_to"])

                if HPM.is_pos_sequence(nh)
                    @test isapprox(rem(θv - θe - vg * π/6, 2π), 0.0, atol=1e-5)
                elseif HPM.is_neg_sequence(nh)
                    @test isapprox(rem(θv - θe + vg * π/6, 2π), 0.0, atol=1e-5)
                elseif HPM.is_zero_sequence(nh)
                    @test isapprox(rem(θv - θe, 2π), 0.0, atol=1e-5)
                end
            end
        end
    end
end