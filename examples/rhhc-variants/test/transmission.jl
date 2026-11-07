@testset "Cardinal $name physical model" for (name, count, form, kv) in (
    ("ieee9", 6, "O230-F", 230.0), ("ieee39", 34, "O345-F", 345.0))
    s = study("cases/$name-physical.toml")
    source = study("cases/$name.toml")
    p = nominal_point(s.box)
    branches = sort(collect(keys(s.details)))
    b = first(branches)
    @test length(s.box) == 10
    @test length(s.details) == count
    @test all(k -> s.model.branches[k].zy === s.model.branches[b].zy, branches)
    spec = specification("Cardinal")
    q = coordinates(s.families.ACSR, (;), (;))
    recipe = construction(spec, q, q.lay_inner, q.lay_outer)
    @test recipe.steel_courses == (6,)
    @test recipe.aluminum_lays == (13, 13, 12)
    for da in (-1.0, 1.0), ds in (-1.0, 1.0), inner in (10.0, 16.0),
        middle in (10.0, 16.0), outer in (10.0, 14.0), steel in (16.0, 26.0)
        point = merge(q,
            (d_aluminium = da, d_steel = ds, lay_inner = inner,
                lay_middle = middle, lay_outer = outer, lay_steel = steel))
        result = realize(spec, point, q)
        @test spec.mass.lower <= result.mass <= spec.mass.upper
        @test spec.diameter.lower <= result.diameter <= spec.diameter.upper
        conductor = design(spec.family, spec, point, s.mat; realized = result)
        @test 2r_ex(conductor.outer) ≈ result.diameter rtol = 1e-12
    end
    # Independent course and placement declaration from the candidate study.
    # HPM parent 3566606. LCM 59cc3f232. Reference run November 7, 2026.
    steel = metal(s.mat.s1a, s.mat.reference; rho = q.rho_steel, mu_r = q.mu_steel)
    aluminum = metal(s.mat.al1, s.mat.reference; rho = q.rho_aluminium)
    conductor = @design "Cardinal-reference" begin
        @terminal :core begin
            stranded(steel; center = Circle(0.00338 / 2), shape = Circle(0.00338 / 2),
                n = (6,), lay = (LayRatio(21),),
                boundary = Circle(3 * 0.00338 / 2), tag = :steel)
            wires(aluminum; shape = Circle(0.00338 / 2), n = 12,
                lay = LayRatio(13), tag = :inner)
            wires(aluminum; shape = Circle(0.00338 / 2), n = 18,
                lay = LayRatio(13), tag = :middle)
            wires(aluminum; shape = Circle(0.00338 / 2), n = 24,
                lay = LayRatio(12), tag = :outer)
        end
    end
    wire = specification("GW-7x3.15")
    shield = design(wire.family, wire, coordinates(s.families.O230, (;), (;)), s.mat)
    reference_system = build(CableSystem,
        [conductor, conductor, conductor, shield, shield],
        [Placement(-10.0, 30 - 20 / 3), Placement(0.0, 30 - 20 / 3),
            Placement(10.0, 30 - 20 / 3), Placement(-3.25, 37.5 - 20 / 3),
            Placement(3.25, 37.5 - 20 / 3)];
        connections = [(core = 1,), (core = 2,), (core = 3,), (core = 0,), (core = 0,)],
        id = "transmission-reference", line_length = 1.0)
    reference = ExactPiModel(
        Dict(1 => (
            problem = (point, f) -> LineParametersProblem(reference_system,
                homogeneous(rho = s.families.EARTH.rho.nominal,
                    eps_r = s.families.EARTH.eps.nominal); frequencies = [f]),
            len = 1.0, zbase = 1.0)),
        CoaxialFormulation(transposition = true),
        s.box)
    for f in (50.0, 250.0, 350.0, 550.0, 650.0)
        actual = s.model.branches[b].zy(p, f)
        expected = reference.branches[1].zy(p, f)
        for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

            @test actual[sequence][quantity] ≈ expected[sequence][quantity] rtol = 1e-12
        end
    end
    v = s.model.branches[b].zy(p, 50.0).positive
    # Reference: scan.jl cases/ieee9-physical.toml, November 7, 2026.
    @test v.z ≈ 6.0257307260845926e-5 + 0.0004357681366383324im rtol = 1e-12
    @test v.y ≈ -2.0289789802773e-17 + 2.636462716284721e-9im rtol = 1e-12
    for k in branches
        br = source.hd["nw"]["1"]["branch"][string(k)]
        detail = s.details[k]
        zbase = kv^2 / source.hd["nw"]["1"]["baseMVA"]
        @test detail.len == br["length_m"]
        pi = branch_data(s.model, k, p, 1, :positive)
        propagation, characteristic = sqrt(v.z * v.y), sqrt(v.z / v.y)
        @test pi.br_r ≈ real(characteristic * sinh(propagation * detail.len)) / zbase rtol = 1e-12
        @test pi.br_x ≈ imag(characteristic * sinh(propagation * detail.len)) / zbase rtol = 1e-12
        @test pi.b_fr ≈ imag(tanh(propagation * detail.len / 2) / characteristic) * zbase rtol = 1e-12
        @test pi.b_fr == pi.b_to
        @test abs(real(v.z) * detail.len / (br["br_r"] * zbase) - 1) < 0.2
        @test abs(imag(v.z) * detail.len / (br["br_x"] * zbase) - 1) < 0.2
        expected_rating = name == "ieee39" && k == 31 ? br["rate_a"] :
                          sqrt(3) * kv * s.config.ratings.Cardinal.current_A /
                          (1000source.hd["nw"]["1"]["baseMVA"])
        @test s.hd["nw"]["1"]["branch"][string(k)]["rate_a"] ≈ expected_rating
    end
    for (id, br) in source.hd["nw"]["1"]["branch"]
        parse(Int, id) in branches && continue
        @test s.hd["nw"]["1"]["branch"][id] == br
    end
    fixed = ParameterBox(
        s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
    prepared = deepcopy(s.hd)
    nominal = solve_rhhc(prepared, s.model, fixed, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt, nbackground = 0, use_seeds = false)
    physical = solve_rhhc(prepared, s.model, s.box, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt, posterior = s.config.posterior)
    @test nominal["signed_off"]
    @test posterior_check(prepared, s.model, fixed, nominal["solution"]; n = 1).max_violation <=
          0
    @test physical["signed_off"]
    @test physical["posterior"].violated_fraction == 0
    @test physical["posterior"].max_violation <= 0
    if name == "ieee39"
        @test Set(k
        for k in branches if !haskey(getproperty(s.config.branches, Symbol(k)), :rating)) ==
              Set([31])
        @test Set(k
        for k in branches
        if max(
            prepared["nw"]["1"]["branch"][string(k)]["cm_fr"],
            prepared["nw"]["1"]["branch"][string(k)]["cm_to"]) >
           sqrt(3) * kv * s.config.ratings.Cardinal.current_A /
           (1000prepared["nw"]["1"]["baseMVA"] * 0.9)) == Set([31])
        @test all(
            k -> max(prepared["nw"]["1"]["branch"][string(k)]["cm_fr"],
                prepared["nw"]["1"]["branch"][string(k)]["cm_to"]) <=
                 prepared["nw"]["1"]["branch"][string(k)]["c_rating"],
            branches)
    end
end
