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

@testset "Declared coordinates" begin
    expanded = ParameterBox(:eps_xlpe => (2.3, 2.0, 2.6, :linear),
        :r_core => (0.015, 0.014, 0.016, :linear),
        :t_xlpe => (0.008, 0.007, 0.009, :linear),
        :rho_earth => (100.0, 50.0, 200.0, :linear),
        :theta => (20.0, 10.0, 30.0, :linear))
    names = (:eps_xlpe, :r_core, :t_xlpe, :rho_earth)
    branches = Dict{Int, NamedTuple}(
        1 => (problem = counted_fixture, len = 100.0,
            zbase = 10.0, coordinates = names),
        2 => (problem = counted_fixture, len = 200.0,
            zbase = 10.0, coordinates = names))
    model = ExactPiModel(branches, formulation, expanded)
    projected = model.branches[1].zy
    full = Circuit(counted_fixture, formulation, expanded)
    @test projected === model.branches[2].zy
    p = nominal_point(expanded)
    for f in (50.0, 650.0)
        result = projected(p, f)
        for value in (-1.0, -0.5, 0.5, 1.0)
            point = merge(p, (theta = latent(expanded, [0.0, 0.0, 0.0, 0.0, value]).theta,))
            @test projected(point, f) === result
            @test projected(point, f) === full(point, f)
        end
    end
    @test length(projected.sequences) == 2
    @test (@inferred projected(p, 50.0)) === projected(p, 50.0)
    @test isconcretetype(keytype(projected.sequences))
    branches[2] = merge(branches[2], (coordinates = (:eps_xlpe,),))
    @test_throws ArgumentError ExactPiModel(branches, formulation, expanded)
    branches[1] = merge(branches[1], (coordinates = (:unknown,),))
    @test_throws ArgumentError ExactPiModel(branches, formulation, expanded)
end
