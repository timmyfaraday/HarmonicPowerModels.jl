################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
################################################################################

@testset "Harmonic Hosting Capacity" begin

    @testset "NLP - Industrial Network" begin
        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = dHHC_NLP

        # define the set of considered harmonics
        H   = [1, 3, 5, 7, 9, 13]
        H⁰  = H[(H.%3).==0]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)

        # solve HHC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver_nlp)

        @testset "Feasibility and Objective" begin
            # solved to optimality
            @test results_hhc["termination_status"] == LOCALLY_SOLVED
            # objective value depending on fairness principle
            if !haskey(hdata, "principle")
                @test isapprox(results_hhc["objective"], 0.160236; atol = 1e-4)
            elseif hdata["principle"] == "equality"
                @test isapprox(results_hhc["objective"], 0.019055; atol = 1e-4)
            end
        end

        @testset "Root Mean Square Voltage" begin 
            # Uminᵢ ≤ RMSᵢ = √(∑ₕ(|Uᵢₕ|²)) ≤ Umaxᵢ, ∀ i ∈ I
            for (nb, bus) ∈ data["bus"]
                vmin    = bus["vmin"]
                vmax    = bus["vmax"]
                
                vm      = [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H]

                @test vmin ⪅ sqrt(sum(vm.^2))
                @test sqrt(sum(vm.^2)) ⪅ vmax
            end
        end
        
        @testset "Total Harmonic Voltage Distortion" begin
            # THDᵢ = √(∑ₕ(|Uᵢₕ|²) / |Uᵢ₁|²) ≤ THDmaxᵢ, ∀ i ∈ I
            for (nb,bus) ∈ data["bus"]
                thdmax  = bus["thdmax"]

                vm_fund = results_hhc["solution"]["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H if nh ≠ 1]

                @test sqrt(sum(vm_harm.^2) / vm_fund^2) ⪅ thdmax
            end
        end
        
        @testset "Individual Harmonic Voltage Distortion" begin
            # IHDᵢₕ = √(|Uᵢₕ|² / |Uᵢ₁|²) ≤ IHDmaxᵢₕ, ∀ i ∈ I, h ∈ H
            for nh ∈ H, (nb,bus) in hdata["nw"]["$nh"]["bus"] if nh ≠ 1
                ihdmax  = bus["ihdmax"]

                vm_fund = results_hhc["solution"]["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"]

                @test sqrt(vm_harm^2 / vm_fund^2) ⪅ ihdmax
            end end
        end

        @testset "Branch Currents" begin
            for nh ∈ H, (nb, branch) in hdata["nw"]["$nh"]["branch"]
                f_bus   = branch["f_bus"]
                t_bus   = branch["t_bus"]

                y_fr    = branch["g_fr"] + im * branch["b_fr"]
                y_to    = branch["g_to"] + im * branch["b_to"]

                z       = branch["br_r"] + im * branch["br_x"]

                v_fr    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vr"] +
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vi"] * im
                v_to    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vi"] * im
                
                cs      = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["csr_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["csi_fr"] * im
                c_fr    = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["cr_fr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["ci_fr"] * im
                c_to    = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["cr_to"] +
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["ci_to"] * im

                # branch total current 
                # Īᵦᵢⱼₕ = Ūᵢ * yˢʰᵦᵢⱼₕ ± Īˢᵦᵢⱼₕ, ∀ βij in Tᵇ, h in H
                @test c_fr ≈ v_fr * y_fr + cs
                @test c_to ≈ v_to * y_to - cs
                # Ohm's law 
                # Īˢᵦᵢⱼₕ = (Ūᵢₕ - Ūⱼₕ) * zᵦ, ∀ βij in Tᵇ⁻ᶠʳ, h in H
                @test v_fr - v_to ≈ z * cs
            end
        end

        @testset "Transformer Voltages and Currents" begin
            for nh ∈ H, (nx, xfmr) in hdata["nw"]["$nh"]["xfmr"]
                f_bus   = xfmr["f_bus"]
                t_bus   = xfmr["t_bus"]

                vg      = parse(Int, xfmr["vg"][3:end])

                cnf1    = xfmr["cnf1"]
                cnf2    = xfmr["cnf2"]

                gnd1    = xfmr["gnd1"]
                gnd2    = xfmr["gnd2"]

                t_vg    = xfmr["tr"] + im * xfmr["ti"]

                u_fr    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vr"] +
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vi"] * im
                u_to    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vi"] * im

                rsh     = xfmr["rsh"]
                z       = xfmr["xsc"] * im

                r_fr    = xfmr["r1"]
                z0_fr   = xfmr["re1"] +
                          xfmr["xe1"] * im
                ysh_fr  = ifelse(cnf1 == 'D' && nh ∈ H⁰, 1 / xfmr["r1"], 0.0)

                r_to    = xfmr["r2"]
                z0_to   = xfmr["re1"] +
                          xfmr["re2"] * im
                ysh_to  = ifelse(cnf2 == 'D' && nh ∈ H⁰, 1 / xfmr["r2"], 0.0)

                e       = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["ert"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["eit"] * im
                θe      = angle(e)

                vt_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vit_fr"] * im
                vt_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vrt_to"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vit_to"] * im
                θv      = angle(vt_to)

                cmt     = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cmrt"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cmit"] * im

                ct_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_fr"] * im
                ct_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_to"] * im 

                cst_fr  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_fr"] * im
                cst_to  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_to"] * im

                ## Core
                # Kirchhoff's current law
                # tᵛᵍₓₕ * (Iˢₓᵢⱼₕ - Īᵐₓₕ + Ēₓₕ/rˢʰₓₕ) + Iₓⱼᵢₕ = 0, ∀ xij ∈ Tˣ⁻ᶠʳ, h ∈ H 
                @test conj(t_vg) * (cst_fr - cmt - e / rsh) + cst_to ≈ 0.0
                # Ohm's law
                # Vₓᵢⱼₕ - Ēₓₕ = zˢ * Iˢₓᵢⱼₕ
                @test vt_fr - e ≈ z * cst_fr
                # transformer phase shift
                # Eₓₕ = tᵛᵍₓₕ * Vₓⱼᵢₕ, ∀ xji ∈ Tˣ⁻ᵗᵒ, h ∈ H 
                @test e ≈ t_vg * vt_to

                ## Winding
                # Current balance
                # Iˢₓᵢⱼₕ = Iₓᵢⱼₕ + yˢʰₓᵢⱼₕ * Vₓᵢⱼₕ, ∀ xij ∈ Tˣ, h ∈ H
                @test cst_fr ≈ ct_fr + ysh_fr * vt_fr 
                @test cst_to ≈ ct_to + ysh_to * vt_to
                # Ohm's law (pos./neg. seq.)
                # Uᵢₕ - Vₓᵢⱼₕ = rₓᵢⱼₕ Iₓᵢⱼₕ, ∀ xij ∈ Tˣ, h ∈ H⁺∪H⁻ 
                if nh ∉ H⁰ @test u_fr - vt_fr ≈ r_fr * ct_fr end 
                if nh ∉ H⁰ @test u_to - vt_to ≈ r_to * ct_to end
                # Ohm's law (zero seq.)
                # Uᵢₕ - Vₓᵢⱼₕ = (rₓᵢⱼₕ + z⁰ₓᵢⱼₕ) Iₓᵢⱼₕ, ∀ xij ∈ Tˣ⁻ᵉᵃʳᵗʰᵉᵈ, h ∈ H⁰
                if gnd1 == 1 && nh ∈ H⁰ @test u_fr - vt_fr ≈ (r_fr + z0_fr) * ct_fr end
                if gnd2 == 1 && nh ∈ H⁰ @test u_to - vt_to ≈ (r_to + z0_to) * ct_to end
                # Zero-seq. current blocking
                # |I|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ᵈᵉˡᵗᵃ, h ∈ H⁰
                if cnf1 == 'D' && nh ∈ H⁰ @test abs(ct_fr) ≈ 0.0 end
                if cnf2 == 'D' && nh ∈ H⁰ @test abs(ct_to) ≈ 0.0 end
                # |I|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ⁿᵉ, h ∈ H⁰
                # |Iˢ|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ⁿᵉ, h ∈ H⁰
                if cnf1 ≠ 'D' && gnd1 == 0 && nh ∈ H⁰ @test abs(ct_fr) ≈ 0.0 end
                if cnf1 ≠ 'D' && gnd1 == 0 && nh ∈ H⁰ @test abs(cst_fr) ≈ 0.0 end
                if cnf2 ≠ 'D' && gnd2 == 0 && nh ∈ H⁰ @test abs(ct_to) ≈ 0.0 end
                if cnf2 ≠ 'D' && gnd2 == 0 && nh ∈ H⁰ @test abs(cst_to) ≈ 0.0 end
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
                if HPM.is_pos_sequence(nh)
                    Δθ = ceil(abs(θv - θe - vg * π/6), digits=10)
                elseif HPM.is_neg_sequence(nh)
                    Δθ = ceil(abs(θv - θe + vg * π/6), digits=10)
                elseif HPM.is_zero_sequence(nh)
                    Δθ = ceil(abs(θv - θe), digits=10)
                end
                @test Δθ % 2pi ≈ 0.0
            end
        end
    end

    @testset "SOC - Industrial Network" begin
        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = dHHC_SOC

        # define the set of considered harmonics
        H   = [1, 3, 5, 7, 9, 13]
        H⁰  = H[(H.%3).==0]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)
        for (n, nw) in hdata["nw"]
            for (b, bus) in nw["bus"]
                bus["angle_range"] = 0.0
            end
        end

        # solve HHC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver_soc, solver_nlp)

        @testset "Feasibility and Objective" begin
            # solved to optimality
            @test results_hhc["termination_status"] == ALMOST_OPTIMAL
            # objective value depending on fairness principle
            if !haskey(hdata, "principle")
                @test isapprox(results_hhc["objective"], 0.160236; atol = 1e-4)
            elseif hdata["principle"] == "equality"
                @test isapprox(results_hhc["objective"], 0.019055; atol = 1e-4)
            end
        end

        @testset "Root Mean Square" begin 
            # Uminᵢ ≤ RMSᵢ = √(∑ₕ(|Uᵢₕ|²)) ≤ Umaxᵢ, ∀ i ∈ I
            for (nb, bus) ∈ data["bus"]
                vmin    = bus["vmin"]
                vmax    = bus["vmax"]
                
                vm      = vcat( hdata["nw"]["1"]["bus"][nb]["vm"], 
                                [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                    for nh ∈ H if nh ≠ 1])

                @test vmin ⪅ sqrt(sum(vm.^2))
                @test sqrt(sum(vm.^2)) ⪅ vmax
            end
        end

        @testset "Total Harmonic Distortion" begin
            # THDᵢ = √(∑ₕ(|Uᵢₕ|²) / |Uᵢ₁|²) ≤ THDmaxᵢ, ∀ i ∈ I
            for (nb,bus) ∈ data["bus"]
                thdmax  = bus["thdmax"]

                vm_fund = hdata["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H if nh ≠ 1]

                @test sqrt(sum(vm_harm.^2) / vm_fund^2) ⪅ thdmax
            end
        end

        @testset "Individual Harmonic Distortion" begin
            # IHDᵢₕ = √(|Uᵢₕ|² / |Uᵢ₁|²) ≤ IHDmaxᵢₕ, ∀ i ∈ I, h ∈ H / 1
            for nh ∈ setdiff(H,1), (nb,bus) in hdata["nw"]["$nh"]["bus"]
                ihdmax  = bus["ihdmax"]

                vm_fund = hdata["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"]

                @test sqrt(vm_harm^2 / vm_fund^2) ⪅ ihdmax
            end
        end

        @testset "Branch Currents" begin
            for nh ∈ H, (nb, branch) in hdata["nw"]["$nh"]["branch"] if nh ≠ 1
                f_bus   = branch["f_bus"]
                t_bus   = branch["t_bus"]

                y_fr    = branch["g_fr"] + im * branch["b_fr"]
                y_to    = branch["g_to"] + im * branch["b_to"]

                z       = branch["br_r"] + im * branch["br_x"]

                v_fr    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vr"] +
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vi"] * im
                v_to    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vi"] * im
                
                cs      = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["csr_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["csi_fr"] * im
                c_fr    = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["cr_fr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["ci_fr"] * im
                c_to    = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["cr_to"] +
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["ci_to"] * im

                # branch total current 
                # Īᵦᵢⱼₕ = Ūᵢ * yˢʰᵦᵢⱼₕ ± Īˢᵦᵢⱼₕ, ∀ βij in Tᵇ, h in H
                @test c_fr ≈ v_fr * y_fr + cs
                @test c_to ≈ v_to * y_to - cs
                # Ohm's law 
                # Īˢᵦᵢⱼₕ = (Ūᵢₕ - Ūⱼₕ) * zᵦ, ∀ βij in Tᵇ⁻ᶠʳ, h in H
                @test v_fr - v_to ≈ z * cs
            end end
        end

        @testset "Transformer Voltages and Currents" begin
            for nh ∈ H, (nx, xfmr) in hdata["nw"]["$nh"]["xfmr"] if nh ≠ 1
                f_bus   = xfmr["f_bus"]
                t_bus   = xfmr["t_bus"]

                vg      = parse(Int, xfmr["vg"][3:end])

                cnf1    = xfmr["cnf1"]
                cnf2    = xfmr["cnf2"]

                gnd1    = xfmr["gnd1"]
                gnd2    = xfmr["gnd2"]

                t_vg    = xfmr["tr"] + im * xfmr["ti"]

                u_fr    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vr"] +
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vi"] * im
                u_to    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vi"] * im

                rsh     = xfmr["rsh"]
                z       = xfmr["xsc"] * im

                r_fr    = xfmr["r1"]
                z0_fr   = xfmr["re1"] +
                          xfmr["xe1"] * im
                ysh_fr  = ifelse(cnf1 == 'D' && nh ∈ H⁰, 1 / xfmr["r1"], 0.0)

                r_to    = xfmr["r2"]
                z0_to   = xfmr["re1"] +
                          xfmr["re2"] * im
                ysh_to  = ifelse(cnf2 == 'D' && nh ∈ H⁰, 1 / xfmr["r2"], 0.0)

                e       = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["ert"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["eit"] * im
                θe      = angle(e)

                vt_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vit_fr"] * im
                vt_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vrt_to"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vit_to"] * im
                θv      = angle(vt_to)

                cmt     = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cmrt"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cmit"] * im

                ct_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_fr"] * im
                ct_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_to"] * im 

                cst_fr  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_fr"] * im
                cst_to  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_to"] * im

                ## Core
                # Kirchhoff's current law
                # tᵛᵍₓₕ * (Iˢₓᵢⱼₕ - Īᵐₓₕ + Ēₓₕ/rˢʰₓₕ) + Iₓⱼᵢₕ = 0, ∀ xij ∈ Tˣ⁻ᶠʳ, h ∈ H 
                @test conj(t_vg) * (cst_fr - cmt - e / rsh) + cst_to ≈ 0.0
                # Ohm's law
                # Vₓᵢⱼₕ - Ēₓₕ = zˢ * Iˢₓᵢⱼₕ
                @test vt_fr - e ≈ z * cst_fr
                # transformer phase shift
                # Eₓₕ = tᵛᵍₓₕ * Vₓⱼᵢₕ, ∀ xji ∈ Tˣ⁻ᵗᵒ, h ∈ H 
                @test e ≈ t_vg * vt_to

                ## Winding
                # Current balance
                # Iˢₓᵢⱼₕ = Iₓᵢⱼₕ + yˢʰₓᵢⱼₕ * Vₓᵢⱼₕ, ∀ xij ∈ Tˣ, h ∈ H
                @test cst_fr ≈ ct_fr + ysh_fr * vt_fr 
                @test cst_to ≈ ct_to + ysh_to * vt_to
                # Ohm's law (pos./neg. seq.)
                # Uᵢₕ - Vₓᵢⱼₕ = rₓᵢⱼₕ Iₓᵢⱼₕ, ∀ xij ∈ Tˣ, h ∈ H⁺∪H⁻ 
                if nh ∉ H⁰ @test u_fr - vt_fr ≈ r_fr * ct_fr end 
                if nh ∉ H⁰ @test u_to - vt_to ≈ r_to * ct_to end
                # Ohm's law (zero seq.)
                # Uᵢₕ - Vₓᵢⱼₕ = (rₓᵢⱼₕ + z⁰ₓᵢⱼₕ) Iₓᵢⱼₕ, ∀ xij ∈ Tˣ⁻ᵉᵃʳᵗʰᵉᵈ, h ∈ H⁰
                if gnd1 == 1 && nh ∈ H⁰ @test u_fr - vt_fr ≈ (r_fr + z0_fr) * ct_fr end
                if gnd2 == 1 && nh ∈ H⁰ @test u_to - vt_to ≈ (r_to + z0_to) * ct_to end
                # Zero-seq. current blocking
                # |I|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ᵈᵉˡᵗᵃ, h ∈ H⁰
                if cnf1 == 'D' && nh ∈ H⁰ @test abs(ct_fr) ≈ 0.0 end
                if cnf2 == 'D' && nh ∈ H⁰ @test abs(ct_to) ≈ 0.0 end
                # |I|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ⁿᵉ, h ∈ H⁰
                # |Iˢ|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ⁿᵉ, h ∈ H⁰
                if cnf1 ≠ 'D' && gnd1 == 0 && nh ∈ H⁰ @test abs(ct_fr) ≈ 0.0 end
                if cnf1 ≠ 'D' && gnd1 == 0 && nh ∈ H⁰ @test abs(cst_fr) ≈ 0.0 end
                if cnf2 ≠ 'D' && gnd2 == 0 && nh ∈ H⁰ @test abs(ct_to) ≈ 0.0 end
                if cnf2 ≠ 'D' && gnd2 == 0 && nh ∈ H⁰ @test abs(cst_to) ≈ 0.0 end
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
                println("nh = $nh, xfmr = $nx")
                if HPM.is_pos_sequence(nh)
                    Δθ = ceil(abs(θv - θe - vg * π/6), digits=10)
                elseif HPM.is_neg_sequence(nh)
                    Δθ = ceil(abs(θv - θe + vg * π/6), digits=10)
                elseif HPM.is_zero_sequence(nh)
                    Δθ = ceil(abs(θv - θe), digits=10)
                end
                @test Δθ % 2pi ≈ 0.0
            end end
        end
    end

    @testset "Equivalence NLP vs SOC - Industrial Network" begin
        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form_nlp = dHHC_NLP
        form_soc = dHHC_SOC

        # define the set of considered harmonics
        H=[1, 3, 5, 7, 9, 13]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)
        for (n, nw) in hdata["nw"]
            for (b, bus) in nw["bus"]
                bus["angle_range"] = 0.0
            end
        end

        # solve HHC problem
        results_hhc_nlp = HPM.solve_hhc(hdata, form_nlp, solver_nlp)
        results_hhc_soc = HPM.solve_hhc(hdata, form_soc, solver_soc, solver_nlp)

        @testset "Feasibility" begin
            @test results_hhc_nlp["termination_status"] == LOCALLY_SOLVED
            @test results_hhc_soc["termination_status"] == ALMOST_OPTIMAL
            @test isapprox(results_hhc_nlp["objective"], results_hhc_soc["objective"], atol=1e-4)
        end
    end
end