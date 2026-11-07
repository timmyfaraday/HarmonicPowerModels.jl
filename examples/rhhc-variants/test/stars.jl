@testset "IEEE30 ACSR sequence references" begin
    # Frozen table added at HPM 85cdb95 from the run recorded at 28bec0e:
    # `scan.jl cases/ieee30-full.toml`. Current comparison uses LCM 59cc3f232.
    s = study("cases/ieee30-full.toml")
    p = nominal_point(s.box)
    checked = Set{String}()
    for line in eachline(joinpath(STUDIES, "output", "ieee30-full-sequences.csv"))
        fields = split(line, ',')
        fields[1] in ("Bear/O132-D", "Goat/O132-D", "Lion/O132-D", "Panther/O132-D") ||
            continue
        f = parse(Float64, fields[2])
        f in (50.0, 250.0, 350.0, 550.0, 650.0) || continue
        b = first(b for (b, row) in s.details if row.circuit == fields[1])
        actual = s.model.branches[b].zy(p, f)
        values0 = parse.(Float64, fields[4:end])
        @test isequal(actual.zero.z, complex(values0[1], values0[2]))
        @test isequal(actual.positive.z, complex(values0[3], values0[4]))
        @test isequal(actual.zero.y, complex(values0[5], values0[6]))
        @test isequal(actual.positive.y, complex(values0[7], values0[8]))
        push!(checked, fields[1])
    end
    @test length(checked) == 4
end

@testset "Falcon $name physical model" for (name, count, length0, source_r, source_x) in ((
    "star3", 2, 25690.406433487706, 1.0, 10.0),
    ("star5", 4, 513.8081286697541, 0.02, 0.2))
    s = study("cases/$name-physical.toml")
    p = nominal_point(s.box)
    @test length(s.box) == 11
    @test length(s.details) == count
    @test all(b -> s.model.branches[b].zy === s.model.branches[1].zy, 1:count)
    @test all(b -> s.details[b].len == length0, 1:count)

    # Independent course declaration from Trefinasa printed p4, with the registered
    # EN50182 lay midpoints. HPM base 5e7400e and LCM 59cc3f232.
    # Reference run: scan.jl cases/star3-physical.toml, November 7, 2026.
    steel = metal(s.mat.s1a, s.mat.reference;
        rho = s.families.ACSR.rho_steel.nominal, mu_r = s.families.ACSR.mu_steel.nominal)
    aluminum = metal(s.mat.al1, s.mat.reference; rho = s.families.ACSR.rho_aluminium.nominal)
    conductor = @design "Falcon-reference" begin
        @terminal :core begin
            stranded(steel; center = Circle(0.00262 / 2), shape = Circle(0.00262 / 2),
                n = (6, 12), lay = (LayRatio(21), LayRatio(18)),
                boundary = Circle(5 * 0.00262 / 2), tag = :steel)
            wires(aluminum; shape = Circle(0.00436 / 2), n = 12,
                lay = LayRatio(13), tag = :inner)
            wires(aluminum; shape = Circle(0.00436 / 2), n = 18,
                lay = LayRatio(13), tag = :middle)
            wires(aluminum; shape = Circle(0.00436 / 2), n = 24,
                lay = LayRatio(12), tag = :outer)
        end
    end
    spec = specification("Falcon")
    q = coordinates(s.families.ACSR, (;), (;))
    for da in (-1.0, 1.0), ds in (-1.0, 1.0), steel in (16.0, 26.0),
        second in (14.0, 22.0), inner in (10.0, 16.0), middle in (10.0, 16.0),
        outer in (10.0, 14.0)
        point = merge(q,
            (d_aluminium = da, d_steel = ds, lay_steel = steel,
                lay_steel_outer = second, lay_inner = inner, lay_middle = middle,
                lay_outer = outer))
        realized = realize(spec, point, q)
        @test spec.mass.lower <= realized.mass <= spec.mass.upper
        @test spec.diameter.lower <= realized.diameter <= spec.diameter.upper
        recipe = construction(spec, point, realized.inner, realized.outer)
        @test issorted((recipe.steel_lays..., recipe.aluminum_lays...); rev = true)
        actual = design(spec.family, spec, point, s.mat; realized)
        @test 2r_ex(actual.outer) ≈ realized.diameter rtol = 1e-12
    end
    layout = formation("O132-D")
    structure = coordinates(s.families.O132, (;), (;))
    wire = specification("GW-7x3.15")
    shield = design(wire.family, wire, structure, s.mat)
    positions = merge(layout, (phase_sag = structure.phase_sag,
        shield_sag = structure.shield_sag, shield))
    system0 = system(positions.kind, conductor, positions)
    reference_problem = (point, f) -> LineParametersProblem(system0,
        homogeneous(rho = s.families.EARTH.rho.nominal,
            eps_r = s.families.EARTH.eps.nominal); frequencies = [f])
    reference = ExactPiModel(
        Dict(1 => (problem = reference_problem, len = s.details[1].len, zbase = 100.0)),
        CoaxialFormulation(transposition = true), s.box)
    for f in (50.0, 250.0, 350.0, 550.0, 650.0)
        actual = s.model.branches[1].zy(p, f)
        expected = reference.branches[1].zy(p, f)
        for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

            @test actual[sequence][quantity] ≈ expected[sequence][quantity] rtol = 1e-12
        end
    end
    v = s.model.branches[1].zy(p, 50.0)
    @test v.positive.z ≈ 3.7248844831307796e-5 + 0.00036143382355627963im rtol = 1e-12
    @test v.positive.y ≈ -5.217612165824774e-18 + 3.1698808627923176e-9im rtol = 1e-12
    @test abs(real(v.positive.z) * length0 / source_r - 1) < 0.2
    @test abs(imag(v.positive.z) * length0 / source_x - 1) < 0.2
    @test s.hd["nw"]["1"]["branch"]["1"]["rate_a"] ≈
          sqrt(3) * 100 * 1372.8443238788311 / 1e5

    for endpoint in (s.box.lower, s.box.nominal, s.box.upper)
        point = NamedTuple{Tuple(s.box.names)}(Tuple(endpoint))
        for h in (1, s.config.harmonics...), spoke in 2:count

            a = branch_data(s.model, 1, point, h, :positive)
            b = branch_data(s.model, spoke, point, h, :positive)
            @test a == b
        end
    end
    branch = branch_data(s.model, 1, p, 1, :positive)
    z, y, length0 = v.positive.z, v.positive.y, s.details[1].len
    propagation, characteristic = sqrt(z * y), sqrt(z / y)
    @test branch.br_r ≈ real(characteristic * sinh(propagation * length0)) / 100 rtol = 1e-12
    @test branch.br_x ≈ imag(characteristic * sinh(propagation * length0)) / 100 rtol = 1e-12
    @test branch.b_fr ≈ imag(tanh(propagation * length0 / 2) / characteristic) * 100 rtol = 1e-12
    @test branch.b_fr == branch.b_to

    fixed = ParameterBox(
        s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
    prepared = deepcopy(s.hd)
    nominal = solve_rhhc(prepared, s.model, fixed, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt,
        nbackground = 0, use_seeds = false)
    @test nominal["signed_off"]
    @test posterior_check(prepared, s.model, fixed, nominal["solution"]; n = 1).max_violation <=
          0
    physical = solve_rhhc(prepared, s.model, s.box, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt, posterior = s.config.posterior)
    @test physical["signed_off"]
    @test physical["posterior"].violated_fraction == 0
    @test physical["posterior"].max_violation <= 0
end
