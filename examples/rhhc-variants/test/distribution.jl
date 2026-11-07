@testset "IEEE33 physical feeder" begin
    s = study("cases/ieee33-physical.toml")
    source = study("cases/ieee33.toml")
    p = nominal_point(s.box)
    @test length(s.details) == 32
    @test length(s.lcm) == 6
    @test length(s.box) == 11
    # Candidate study conductor library and strand patterns, HPM parent 587dd89,
    # LCM 59cc3f232. Reference construction run November 7, 2026.
    declarations = (("Rabbit", 1, 0.00335, 0.00335, (6,)),
        ("Weasel", 1, 0.00259, 0.00259, (6,)),
        ("Dog", 7, 0.00157, 0.00472, (6,)),
        ("Tiger", 7, 0.00236, 0.00236, (12, 18)),
        ("Wolf", 7, 0.00259, 0.00259, (12, 18)),
        ("Bison", 7, 0.003, 0.003, (12, 18, 24)))
    q = coordinates(s.families.ACSR, (;), (;))
    steel = metal(s.mat.s1a, s.mat.reference; rho = q.rho_steel, mu_r = q.mu_steel)
    aluminum = metal(s.mat.al1, s.mat.reference; rho = q.rho_aluminium)
    for (name, count, ds, da, courses) in declarations
        spec = specification(name)
        radius = count == 1 ? ds / 2 : 3ds / 2
        first(courses) == 6 && (radius = max(radius, da / 2))
        core = stranded(steel; center = Circle(ds / 2), shape = Circle(ds / 2),
            n = count == 1 ? () : (6,), lay = count == 1 ? () : (LayRatio(21),),
            boundary = Circle(radius))
        layers = Tuple(wires(aluminum; shape = Circle(da / 2), n = n,
                           lay = LayRatio(i == length(courses) ? 12 : 13))
        for (i, n) in enumerate(courses))
        conductor = build(CableDesign, name, Terminal(:core, Stack(core, layers...)))
        positions = [Placement(x, 10.0) for x in (-0.8, 0.0, 0.8)]
        direct = build(CableSystem, fill(conductor, 3), positions;
            connections = [(core = k,) for k in 1:3], line_length = 1.0, id = name)
        layout = merge(formation("O13-F"), (phase_sag = 1.5,))
        unshielded = system(layout.kind, conductor, layout)
        empty_wires = merge(layout,
            (earth_wire = "GW-7x3.15", shield_x = (), shield_height = (),
                shield_sag = 1.5, shield = conductor))
        zero_wires = system(layout.kind, conductor, empty_wires)
        references = [ExactPiModel(
                          Dict(1 => (
                              problem = (point, f) -> LineParametersProblem(circuit,
                                  homogeneous(rho = s.families.EARTH.rho.nominal,
                                      eps_r = s.families.EARTH.eps.nominal); frequencies = [f]),
                              len = 1.0, zbase = 1.0)),
                          CoaxialFormulation(transposition = true), s.box)
                      for circuit in (direct, unshielded, zero_wires)]
        branches = sort([b for (b, detail) in s.details if detail.spec == name])
        @test all(b -> s.model.branches[b].zy === s.model.branches[first(branches)].zy, branches)
        for f in (50.0, 250.0, 350.0, 550.0, 650.0)
            actual = s.model.branches[first(branches)].zy(p, f)
            for reference in references
                expected = reference.branches[1].zy(p, f)
                for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

                    @test actual[sequence][quantity] ≈ expected[sequence][quantity] rtol = 1e-12
                end
            end
        end
        for da0 in (-1.0, 1.0), ds0 in (-1.0, 1.0), lay in (10.0, 14.0)
            point = merge(q, (d_aluminium = da0, d_steel = ds0, lay_outer = lay))
            result = realize(spec, point, q)
            actual = design(spec.family, spec, point, s.mat; realized = result)
            @test 2r_ex(actual.outer) ≈ result.diameter rtol = 1e-12
            @test spec.diameter.lower <= result.diameter <= spec.diameter.upper
            @test spec.mass.lower <= result.mass <= spec.mass.upper
        end
        v = s.model.branches[first(branches)].zy(p, 50.0).positive
        @test imag(v.y) > 0
        for b in branches
            detail = s.details[b]
            route = getproperty(s.config.routes, Symbol(b))
            @test detail.len == route.length_m
            br = source.hd["nw"]["1"]["branch"][string(b)]
            zbase = 12.66^2 / source.hd["nw"]["1"]["baseMVA"]
            pi = branch_data(s.model, b, p, 1, :positive)
            propagation, characteristic = sqrt(v.z * v.y), sqrt(v.z / v.y)
            @test pi.br_r ≈ real(characteristic * sinh(propagation * detail.len)) / zbase rtol = 1e-12
            @test pi.br_x ≈ imag(characteristic * sinh(propagation * detail.len)) / zbase rtol = 1e-12
            @test pi.b_fr ≈
                  imag(tanh(propagation * detail.len / 2) / characteristic) * zbase rtol = 1e-12
            @test pi.b_fr == pi.b_to
            @test s.hd["nw"]["1"]["branch"][string(b)]["rate_a"] ≈
                  sqrt(3) * 12.66 * getproperty(s.config.ratings, Symbol(name)).current_A /
                  (1000source.hd["nw"]["1"]["baseMVA"])
            @test br["b_fr"] == br["b_to"] == 0
        end
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
    @test all(br -> max(br["cm_fr"], br["cm_to"]) <= br["c_rating"],
        values(prepared["nw"]["1"]["branch"]))
end

@testset "Industrial physical cables" begin
    s = study("cases/industrial-physical.toml")
    source = study("cases/industrial.toml")
    p = nominal_point(s.box)
    @test length(s.details) == length(s.lcm) == 3
    @test length(s.box) == 30
    @test s.hd["nw"]["1"]["xfmr"] == source.hd["nw"]["1"]["xfmr"]
    for (b, name, len, rating, terminals) in ((2, "U36-Al400", 1500.0, 24.5, 2),
        (3, "U11-Al240", 200.0, 6.58, 2),
        (4, "NYY-Cu240-120", 70.0, 0.382, 1))
        detail = s.details[b]
        @test detail.spec == name
        @test detail.len == len
        @test s.hd["nw"]["1"]["branch"][string(b)]["rate_a"] ≈ rating / 100
        spec = specification(name)
        fn = problem(spec, s.families, formation("trefoil"), s.mat,
            (;
                (k => (; (n => Symbol(v) for (n, v) in pairs(row))...)
            for (k, row) in pairs(s.config.bindings))...))
        cable = first(fn(p, 50.0).system.designs)
        @test length(cable.terminals) == terminals
        q = coordinates(getproperty(s.families, spec.family_name), (;), (;))
        # Candidate study resistance caps, HPM parent 23f0b26, LCM 59cc3f232.
        # DC reference at 20 degrees Celsius from the industrial cable library.
        cold = merge(p, NamedTuple{(Symbol("theta_", spec.family_name),)}((20.0,)))
        dc = s.model.branches[b].zy(cold, 0.01).positive
        @test real(dc.z) ≈ q.R20 rtol = 1e-4
        for f in (50.0, 250.0, 350.0, 550.0, 650.0)
            v = s.model.branches[b].zy(p, f)
            @test v.negative.z ≈ v.positive.z rtol = 1e-12
            @test v.negative.y ≈ v.positive.y rtol = 1e-12
            @test imag(v.positive.y) > 0
            pi = branch_data(s.model, b, p, f / 50, :positive)
            gamma = sqrt(v.positive.z * v.positive.y)
            characteristic = sqrt(v.positive.z / v.positive.y)
            expected = characteristic * sinh(gamma * len) / s.model.branches[b].zbase
            @test complex(pi.br_r, pi.br_x) ≈ expected rtol = 1e-12
        end
        changed = merge(p, NamedTuple{(Symbol("theta_", spec.family_name),)}((20.0,)))
        @test real(s.model.branches[b].zy(changed, 50.0).positive.z) <
              real(s.model.branches[b].zy(p, 50.0).positive.z)
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
    @test all(br -> max(br["cm_fr"], br["cm_to"]) <= br["c_rating"],
        values(prepared["nw"]["1"]["branch"]))
end
