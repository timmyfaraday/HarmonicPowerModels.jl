@testset "Circuit" begin
    p = nominal_point(box)
    circuit = Circuit(fixture, formulation, box)
    for f in (50.0, 650.0)
        result = circuit(p, f)
        for sequence in (:zero, :positive, :negative)
            reference = direct(p, f)[sequence]
            @test isapprox(result[sequence].z, reference.z; rtol = 1e-12)
            @test isapprox(result[sequence].y, reference.y; rtol = 1e-12)
        end
    end
    @test (@inferred circuit(p, 50.0)) == circuit(p, 50.0)

    fixture_calls[] = 0
    counted = Circuit(counted_fixture, formulation, box)
    first_result = counted(p, 50.0)
    @test counted(p, 50.0) == first_result
    @test length(counted.sequences) == 1
    @test fixture_calls[] == 1

    points = [latent(box, fill(-1.0 + 2(i - 1) / 63, 4)) for i in 1:64]
    serial = [circuit(point, 650.0) for point in points]
    threaded = Circuit(fixture, formulation, box)
    results = similar(serial)
    Threads.@threads for i in eachindex(points)
        results[i] = threaded(points[i], 650.0)
    end
    @test results == serial

    untransposed = Circuit(flat_fixture, CoaxialFormulation(transposition = false), box)
    @test_logs (:warn, r"Sequence matrix is not diagonal within tolerance") match_mode = :any untransposed(
        p, 50.0)
    @test_nowarn Circuit(fixture, formulation, box)(p, 50.0)
end
