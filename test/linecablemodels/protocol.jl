@testset "HPM protocol" begin
    branches = Dict(
        1 => (problem = fixture, len = 1000.0, zbase = 10.89),
        2 => (problem = fixture, len = p -> 1500.0 * p.t_xlpe / 8.0e-3, zbase = 10.89),
        3 => (problem = second_fixture, len = 2000.0, zbase = 10.89))
    model = ExactPiModel(branches, formulation, box)
    @test model.branches[1].zy === model.branches[2].zy
    @test model.branches[1].zy !== model.branches[3].zy
    @test length(unique(br.zy for br in values(model.branches))) == 2

    reference = ExactPiModel(Dict(b => (zy = (p, f, s) -> direct(p, f)[s],
                                      len = br.len, zbase = br.zbase)
    for (b, br) in branches))
    p = nominal_point(box)
    for b in keys(branches), h in (1, 5, 7), sequence in (:zero, :positive, :negative)
        actual = branch_data(model, b, p, h, sequence)
        expected = branch_data(reference, b, p, h, sequence)
        @test isapprox(complex(actual.br_r, actual.br_x),
            complex(expected.br_r, expected.br_x); rtol = 1e-12)
        @test isapprox(complex(actual.g_fr, actual.b_fr),
            complex(expected.g_fr, expected.b_fr); rtol = 1e-12)
        @test isapprox(complex(actual.g_to, actual.b_to),
            complex(expected.g_to, expected.b_to); rtol = 1e-12)
    end
    zy = model.branches[1].zy
    expected = direct(p, 50.0).zero
    @test all(isapprox(a, e; rtol = 1e-12) for (a, e) in zip(zy(p, 50.0, :zero), expected))
    @test_throws ArgumentError("sequence :foo is unavailable; Circuit serves :zero, :positive and :negative") zy(
        p, 50.0, :foo)
    # Olex 19/33 kV single-core, 8.0 mm XLPE, 400-800 mm²: 0.267-0.366 µF/km
    # (`local/physical-inputs/inputs.md`, U33 family table).
    @test 0.25e-9 <= imag(zy(p, 50.0, :positive).y) / (2π * 50) <= 0.40e-9
end
