using Test
include(joinpath(@__DIR__, "..", "run.jl"))

@testset "IEEE30 physical-input study" begin
    for name in ("ieee30-tom", "ieee30-u33", "ieee30-physical", "ieee30-full")
        @testset "$name" begin
            s = study("cases/$name.toml")
            result = runstudy(s, devnull)
            nom, res = result.nominal, result.result
            @testset "Zero-width box" begin
                @test nom["iterations"] == 1
                @test nom["signed_off"]
                @test nom["lower_bound"] ≈ nom["upper_bound"] rtol = 1e-3
                @test all(isapprox(nom["solution"][h], nom["master_solution"][h]; rtol = 1e-3)
                for h in s.config.harmonics)
            end
            @testset "50 Hz per-unit magnitudes" begin
                p = nominal_point(s.box)
                for b in circuits(s)
                    d = s.details[b]
                    row = branch_data(s.model, b, p, 1.0, :positive)
                    per_km = d.len / 1e3 / s.model.branches[b].zbase
                    r, x = row.br_r / per_km, row.br_x / per_km
                    c = (row.b_fr + row.b_to) /
                        (2π * s.config.f1 * s.model.branches[b].zbase * d.len) * 1e9
                    println(name, " 50 Hz branch ", b, ": ", row,
                        "; R/X/C per km: ", (r, x, c))
                    if d.kv == 33
                        @test 0.01 < r < 0.2
                        @test 0.05 < x < 0.2
                        @test 0.15 < c < 0.5
                    else
                        @test 0.01 < r < 0.3
                        @test 0.2 < x < 0.6
                        @test 0.001 < c < 0.03
                    end
                end
            end
            @testset "Capacity bounds and posterior" begin
                @test res["lower_bound"] <= nom["upper_bound"]
                for h in s.config.harmonics
                    @test res["solution"][h] / nom["master_solution"][h] <= 1
                end
                @test res["signed_off"]
                @test res["posterior"].violated_fraction == 0
            end
        end
    end
end
