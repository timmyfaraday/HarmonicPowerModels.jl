
@testset "harmonic transformer model" begin 

    @testset "magnetizing current" begin
        # BH-curve - Thyssenkrupp - PowerCore H 100-23 50Hz
        B⁺ = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
        H⁺ = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
        B = vcat(reverse(-B⁺),0.0,B⁺)
        H = vcat(reverse(-H⁺),0.0,H⁺) 
        BH_powercore_h100_23 = Dierckx.Spline1D(B, H; k=3, bc="nearest")

        # xfmr magnetizing data
        Hᴱ = [1,3]
        Hᴵ = collect(1:2:25)
        Fᴱ = :rectangular
        Fᴵ = :rectangular
        l  = 10.0
        A  = 0.475
        N  = 45
        Vp = 12470
        Vs = 4160
        Ip = 100e6 / 12470

        # excitation voltage - rectangular - H ∈ [1,3]
        Ere = [1.00, 0.30]
        Eim = [0.20, 0.05]
        # excitation voltage - polar - H ∈ [1,3]
        E = hypot.(Ere, Eim)
        θ = atan.(Ere, Eim)

        # manual calculation of the magnitizing current
        f  = HPM.freq
        ω  = 2 * π * f
        dt = 1 / (100 * f * maximum(Hᴵ))
        t  = 0.0:dt:(5.0 / f)
        B  = sum(Vp .* E[ni] ./ (ω * nh * A * N) .* cos.(ω .* nh .* t .+ θ[ni]) for (ni,nh) in enumerate(Hᴱ))
        H  = BH_powercore_h100_23(B)
        Im = l ./ (N .* Ip)  .* H

        # decomposition of the manual calculation
        fq = SDC.Sinusoidal(f .* Hᴵ)
        SDC.decompose(t, Im, fq)

        # 
        HPM.sample_xfmr_magnetizing_current(;xfmr_magn=magn)
    end



    @testset "YN - yn - 0" begin
        path  = joinpath(HPM.BASE_DIR,"test/data/matpower/case_xfmr_YNyn0.m")
        data  = PMs.parse_file(path)
        hdata = HPM.replicate(data)

        solve_hopf(hdata, PMs.IVRPowerModel, solver)
    end

    @testset "D - yn - 11" begin
        path  = joinpath(HPM.BASE_DIR,"test/data/matpower/case_xfmr_Dyn11.m")
        data  = PMs.parse_file(path)
        hdata = HPM.replicate(data)

        solve_hopf(hdata, PMs.IVRPowerModel, solver)
    end

    @testset "Y - d - 11" begin
        path  = joinpath(HPM.BASE_DIR,"test/data/matpower/case_xfmr_Yd11.m")
        data  = PMs.parse_file(path)
        hdata = HPM.replicate(data)

        solve_hopf(hdata, PMs.IVRPowerModel, solver)
    end

    @testset "YN - d - 11" begin
        path  = joinpath(HPM.BASE_DIR,"test/data/matpower/case_xfmr_YNd11.m")
        data  = PMs.parse_file(path)
        hdata = HPM.replicate(data)

        solve_hopf(hdata, PMs.IVRPowerModel, solver)
    end
end