using Test
include(joinpath(@__DIR__, "..", "reduce.jl"))

@testset "Equivalent designs" begin
    inputs = reduction_inputs("cases/ieee30-full.toml")
    points = source_points(inputs.box)
    formulation = CoaxialFormulation(transposition = true)
    maximum_blueprint = 0.0
    maximum_zy = 0.0
    for row in inputs.definitions
        source = ExactPiModel(Dict(1 => (problem = row.problem, len = 1.0, zbase = 1.0)),
            formulation, inputs.box)
        reduced = ExactPiModel(
            Dict(1 => (problem = equivalent_problem(row.problem),
                len = 1.0, zbase = 1.0)),
            formulation, inputs.box)
        radial = ExactPiModel(
            Dict(1 => (
                problem = (p, f) -> begin
                    original = row.problem(p, f)
                    equivalent(
                        original, Dict(equivalent_values(original, row.layout.kind)),
                        row.layout.kind)
                end,
                len = 1.0,
                zbase = 1.0)),
            formulation,
            inputs.box)
        @testset "$(row.name)" begin
            for p in points
                original = row.problem(p, 50.0)
                replacement = equivalent(original)
                for (a, b) in zip(original.system.designs, replacement.system.designs)
                    difference = blueprint_error(a, b)
                    maximum_blueprint = max(maximum_blueprint, difference)
                    @test difference <= 1e-12
                end
                for f in (50.0, 650.0)
                    a = source.branches[1].zy(p, f)
                    b = reduced.branches[1].zy(p, f)
                    c = radial.branches[1].zy(p, f)
                    for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

                        difference = relative_error(getproperty(a[sequence], quantity),
                            getproperty(b[sequence], quantity))
                        maximum_zy = max(maximum_zy, difference)
                        @test difference <= 1e-12
                        @test relative_error(a[sequence][quantity], c[sequence][quantity]) <=
                              1e-12
                    end
                end
            end
        end
        println(row.name, ": checked ", length(points), " source points at 50 and 650 Hz")
        flush(stdout)
    end
    println("Maximum blueprint relative error: ", maximum_blueprint)
    println("Maximum sequence relative error: ", maximum_zy)
end

@testset "Effective inputs" begin
    inputs = reduction_inputs("cases/ieee30-full.toml")
    points = [latent(inputs.box, x)
              for x in background_points(inputs.box, METHOD.verification_samples; seed = 2)]
    for row in filter(r -> r.spec.family_name == :ACSR, inputs.definitions)
        @testset "$(row.name) feasible image" begin
            for p in points
                mapped = physical_inputs(row.spec, inputs.families, row.layout,
                    inputs.bindings, p)
                @test feasible(row, mapped)
            end
        end
    end
    nominal = nominal_point(inputs.box)
    cold = merge(nominal,
        (theta_ACSR = inputs.families.ACSR.theta.lower,
            ambient_O132 = inputs.families.O132.ambient.lower))
    hot = merge(cold, (theta_ACSR = inputs.families.ACSR.theta.upper,))
    tied = merge(nominal, (lay_outer_ACSR = inputs.families.ACSR.lay_outer.upper,))
    for row in filter(r -> r.spec.family_name == :ACSR, inputs.definitions)
        @testset "$(row.name) intrinsic height sensitivities" begin
            for fraction in METHOD.difference_step .* (2, 1, 0.5)
                result = sensitivities(inputs, row, [cold, hot], :effective, fraction)
                for (j, axis) in enumerate(result.axes)
                    axis.group == :realized || continue
                    for cable in 1:3
                        k = findfirst(q -> first(q) == "cable.$cable.y", result.reference)
                        @test result.values[k, j, 1] * result.values[k, j, 2] < 0
                    end
                end
            end
        end
        @testset "$(row.name) lay-order map" begin
            source = sensitivities(
                inputs, row, [nominal, tied], :source, METHOD.difference_step)
            physical = sensitivities(
                inputs, row, [nominal, tied], :effective, METHOD.difference_step)
            j = findfirst(a -> a.name == :lay_inner_ACSR, source.axes)
            for cable in 1:3
                k = findfirst(q -> first(q) == "cable.$cable.tube.1.mu_r", source.reference)
                @test source.values[k, j, 1] * source.values[k, j, 2] < 0
                @test all(>(0), physical.values[k, j, :])
            end
        end
    end
end

@testset "Subdivision and containment" begin
    inputs = reduction_inputs("cases/ieee30-full.toml")
    result = hull_study(inputs, "/tmp/hpm-equivalent-tests")
    global tested_hulls = result
    for circuit in result.results
        @testset "$(circuit.row.name)" begin
            @test isempty(circuit.checked.escaped)
            for cells in circuit.hulls.cells, cell in cells

                @test cell.consistent
                @test cell.depth <= METHOD.subdivision_depth
            end
            if circuit.row.spec.family_name == :ACSR
                j = findfirst(a -> a.group == :temperature, circuit.data.axes)
                reference = first(circuit.data.physical).structure.reference_theta
                for cable in 1:3
                    k = findfirst(p -> first(p) == "raw.cable.$cable.y", circuit.data.source.reference)
                    cells = circuit.hulls.cells[k]
                    @test length(cells) == 2
                    @test maximum(c -> c.depth, cells) == 1
                    @test abs(first(cells).upper[j] - reference) <= 1e-5
                end
            end
        end
    end
    groups = equivalent_coordinates(result)
    for group in groups, a in group.members, b in group.members
        @test maximum(abs.(a.profile .- b.profile)) <= METHOD.sharing_tolerance
    end
end

@testset "Boxed equivalent combination" begin
    inputs = tested_hulls.inputs
    groups = equivalent_coordinates(tested_hulls)
    a = first(g
    for g in groups
    if g.family == :ACSR && g.quantity == "resistance" &&
           any(m -> m.row == 1 && m.key == "cable.1.tube.1.resistance", g.members))
    b = first(g
    for g in groups
    if g.family == :ACSR && g.quantity == "y" &&
           any(m -> m.row == 1 && m.key == "cable.1.y", g.members))
    point = NamedTuple{(a.name, b.name)}((1.0, 1.0))
    circuit = equivalent_problem(inputs, 1, [a, b], point, 50.0)
    values = Dict(equivalent_values(circuit, first(inputs.definitions).layout.kind))
    for group in (a, b), member in group.members

        member.row == 1 || continue
        @test relative_error(member.upper, values[member.key]) <= 1e-12
    end
end

@testset "Generated equivalent case" begin
    inputs = tested_hulls.inputs
    evidence = "/tmp/hpm-equivalent-tests"
    proposal(tested_hulls, evidence)
    fresh = generate_equivalent(inputs, evidence;
        destination = "/tmp/hpm-equivalent-regenerated")
    for path in keys(fresh.files)
        @test read(joinpath(STUDIES, path)) ==
              read(joinpath("/tmp/hpm-equivalent-regenerated", path))
    end
    equivalent_inputs = reduction_inputs("cases/ieee30-full-equivalent.toml")
    source_case = TOML.parsefile(inputs.file)
    equivalent_case = TOML.parsefile(equivalent_inputs.file)
    @test source_case["network"] == equivalent_case["network"]
    @test source_case["branches"] == equivalent_case["branches"]
    dataset = readrecord(joinpath(STUDIES, equivalent_case["equivalent_data"]))
    formulation = CoaxialFormulation(transposition = true)
    original_model = ExactPiModel(
        Dict(i => (problem = row.problem, len = 1.0, zbase = 1.0)
        for (i, row) in enumerate(inputs.definitions)),
        formulation,
        inputs.box)
    generated_model = ExactPiModel(
        Dict(i => (problem = row.problem, len = 1.0, zbase = 1.0)
        for (i, row) in enumerate(equivalent_inputs.definitions)),
        formulation,
        equivalent_inputs.box)
    for i in eachindex(inputs.definitions), frequency in (50.0, 650.0)

        a = original_model.branches[i].zy(nominal_point(inputs.box), frequency)
        b = generated_model.branches[i].zy(nominal_point(equivalent_inputs.box), frequency)
        for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

            @test relative_error(a[sequence][quantity], b[sequence][quantity]) <= 1e-12
        end
    end
    for row in inputs.definitions
        members = [member for coordinate in values(dataset.coordinates)
                   for member in coordinate.members if member.circuit == row.name]
        escaped = 0
        for point in tested_hulls.verification
            observed = Dict(equivalent_values(row.problem(point, 50.0), row.layout.kind))
            escaped += count(m -> !(m.lower <= observed[m.field] <= m.upper), members)
        end
        @test escaped == 0
    end
    for point in tested_hulls.verification[1:32]
        observed = Dict(row.name =>
                            Dict(equivalent_values(row.problem(point, 50.0), row.layout.kind))
        for row in inputs.definitions)
        encoded = NamedTuple{Tuple(equivalent_inputs.box.names)}(Tuple(begin
                                                                           coordinate = getproperty(
                                                                               dataset.coordinates,
                                                                               name)
                                                                           member = first(coordinate.members)
                                                                           normalized(
                                                                               observed[member.circuit][member.field],
                                                                               member.nominal,
                                                                               member.lower,
                                                                               member.upper)
                                                                       end
        for name in equivalent_inputs.box.names))
        for row in equivalent_inputs.definitions
            actual = Dict(equivalent_values(row.problem(encoded, 50.0), row.layout.kind))
            for (name, coordinate) in pairs(dataset.coordinates),
                member in coordinate.members

                member.circuit == row.name || continue
                expected = interpolate(member, encoded[name])
                @test relative_error(expected, actual[member.field]) <= 1e-12
            end
        end
    end
end

include(joinpath(@__DIR__, "..", "compare.jl"))

@testset "Corner adversary" begin
    s = study("cases/ieee30-tom.toml")
    fixed = ParameterBox(
        s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
    prepared = deepcopy(s.hd)
    nominal = solve_rhhc(prepared, s.model, fixed, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt,
        nbackground = 0, use_seeds = false)
    solution = nominal["solution"]
    result = corner_adversary(prepared, s.model, s.box, solution)
    independently_checked = [posterior_check(prepared, s.model,
                                 ParameterBox(
                                     s.box.names, point, point, point, s.box.scale), solution; n = 1)
                             for point in result.corners]
    worst = argmax([p.max_violation for p in independently_checked])
    @test result.max_violation ≈ independently_checked[worst].max_violation rtol = 1e-12
    @test result.constraint == independently_checked[worst].constraint
    @test result.count > 1
    @test any(!iszero, result.sensitivities.signs)
    stationary = corner_adversary(prepared, s.model, fixed, solution)
    direct = posterior_check(prepared, s.model, fixed, solution; n = 1)
    @test stationary.count == 1
    @test stationary.point == nominal_point(fixed)
    @test stationary.max_violation ≈ direct.max_violation rtol = 1e-12
    doubled = corner_adversary(prepared, s.model, fixed, Dict(h => 2current
    for (h, current) in solution))
    @test doubled.max_violation ≈ 4(stationary.max_violation + 1) - 1 rtol = 1e-12
end

include("surrogates.jl")

include("stars.jl")
