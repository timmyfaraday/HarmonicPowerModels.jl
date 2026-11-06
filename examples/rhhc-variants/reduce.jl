include(joinpath(@__DIR__, "builders.jl"))
using SHA: sha256

"""
    reduction_inputs(path)

Read the case's selected source box and every distinct spec and formation.
Build the same problem functions as the network studies, without solving a network.
"""
function reduction_inputs(path)
    file = isfile(path) ? abspath(path) : joinpath(STUDIES, path)
    config = readrecord(file)
    mat = readrecord(joinpath(STUDIES, "data", "materials.toml"))
    records = readrecord(config)
    families = records.families
    box = parameterbox(config.box, families, config.box_names)
    bindings = (;
        (k => (; (n => Symbol(v) for (n, v) in pairs(b))...)
    for (k, b) in pairs(config.bindings))...)
    rows = sort(unique([(row.spec, row.formation) for row in values(config.branches)]))
    definitions = map(rows) do (name, placement)
        spec = specification(name, placement, records)
        row = formation(placement)
        layout = haskey(row, :earth_wire) ?
                 merge(row, (shield_spec = specification(row.earth_wire),)) : row
        (; name = name * "/" * placement, spec, layout,
            problem = problem(spec, families, row, mat, bindings))
    end
    return (; file, box, definitions, families, mat, bindings, config)
end

"""
    source_points(box)

Return the nominal point, the two endpoints of every source axis with other
coordinates nominal, and 256 Sobol points. Use HPM's seed 1 and box mapping.
"""
function source_points(box)
    p = nominal_point(box)
    points = [p]
    for (i, key) in enumerate(box.names), bound in (box.lower, box.upper)

        push!(points, merge(p, NamedTuple{(key,)}((bound[i],))))
    end
    append!(points, [latent(box, x) for x in background_points(box, 256; seed = 1)])
    return points
end

"""
    equivalent(source::LineParametersProblem)

Replace every design with its plain-tube equivalent. Keep the placements,
connections, earth, temperature and frequencies of the source problem.
"""
function equivalent(source::LineParametersProblem)
    circuit = source.system
    connections = [(;
                       (design.terminals[row.terminal].name => row.connection
                   for row in circuit.terminals if row.cable == k)...)
                   for (k, design) in enumerate(circuit.designs)]
    reduced = build(CableSystem, equivalent.(circuit.designs), circuit.placements;
        connections, id = circuit.id, line_length = circuit.line_length)
    return LineParametersProblem(reduced, source.earth;
        temperature = source.temperature, frequencies = source.frequencies)
end

"""
    quantities(problem)

Read the blueprint quantities and external geometry used by the coaxial computation.
Keys identify cable, terminal or layer, and physical field. Values use SI units. The
relative material properties and the lay correction are dimensionless. The public formulation law evaluates resistivity at the problem temperature.
"""
function quantities(problem)
    rows = Pair{String, Float64}[]
    formulation = CoaxialFormulation(transposition = true)
    for (k, design) in enumerate(problem.system.designs)
        blueprint = LineCableModels.flatten(design)
        for (i, tube) in enumerate(blueprint.conductors)
            prefix = "cable.$k.tube.$i."
            append!(rows,
                [prefix * "r_in" => tube.r_in, prefix * "r_ex" => tube.r_ex,
                    prefix * "rho" =>
                        formulation.resistivity(tube.material, problem.temperature),
                    prefix * "mu_r" => tube.material.mu_r,
                    prefix * "mu_r_ins" => blueprint.mu_r_ins[i]])
        end
        for (i, layer) in enumerate(blueprint.dielectrics)
            prefix = "cable.$k.layer.$i."
            append!(rows,
                [prefix * "r_in" => layer.r_in, prefix * "r_ex" => layer.r_ex,
                    prefix * "eps_r" => layer.material.eps_r,
                    prefix * "rho" => layer.material.rho,
                    prefix * "mu_r" => layer.material.mu_r,
                    prefix * "tan_delta" => layer.material.tan_delta])
        end
        center = problem.system.placements[k] * Placement(only(blueprint.positions)...)
        append!(rows,
            ["cable.$k.x" => center.x, "cable.$k.y" => center.y,
                "cable.$k.radius" => support(design.outer)])
    end
    for (k, layer) in enumerate(problem.earth.layers),
        field in (:rho, :eps_r, :mu_r, :thickness)

        push!(rows, "earth.$k.$field" => getproperty(layer, field))
    end
    return rows
end

"""
    blueprint_error(source, reduced)

Return the largest relative difference among numerical blueprint quantities.
Zero values must agree exactly. Topology and terminal ordering must agree.
"""
function blueprint_error(source, reduced)
    a = LineCableModels.flatten(source)
    b = LineCableModels.flatten(reduced)
    for field in (:dielectric_ranges, :assembly_ranges, :interior_ranges, :positions)
        getproperty(a, field) == getproperty(b, field) || return Inf
    end
    values = Float64[]
    for field in (:conductors, :dielectrics)
        length(getproperty(a, field)) == length(getproperty(b, field)) || return Inf
        for (x, y) in zip(getproperty(a, field), getproperty(b, field))
            for key in (:r_in, :r_ex)
                push!(values, relative_error(getproperty(x, key), getproperty(y, key)))
            end
            for key in fieldnames(typeof(x.material))
                push!(values, relative_error(getproperty(x.material, key), getproperty(y.material, key)))
            end
        end
    end
    append!(values, relative_error.(a.mu_r_ins, b.mu_r_ins))
    return maximum(values; init = 0.0)
end

"""
    relative_error(reference, value)

Return the relative magnitude of a difference, with exact agreement for zeros
and infinite material properties.
"""
function relative_error(reference, value)
    reference == value ? 0.0 :
    (iszero(reference) ? Inf : abs(value - reference) / abs(reference))
end

"""
    equivalent_problem(definition)

Build a callable with the source definition's point and frequency arguments.
The result uses the equivalent design path.
"""
function equivalent_problem(definition)
    return (p, f) -> equivalent(definition(p, f))
end

"""
    effective_axes(inputs, row)

Describe the effective coordinates driven by the selected source box. Actual
wire diameters use meters. Lays and the U33 registered positions are dimensionless.
Temperatures use degrees Celsius. Other coordinates keep their source SI units.
"""
function effective_axes(inputs, row)
    axes = NamedTuple[]
    for (j, name) in enumerate(inputs.box.names)
        binding = getproperty(inputs.config.box, name)
        family = Symbol(binding.family)
        field = Symbol(binding.coordinate)
        group = family == :EARTH ? :ground : :q
        relevant = family == :EARTH || family == row.spec.family_name
        if haskey(row.layout, :structure) && family == row.layout.structure
            group = :structure
            relevant = true
        end
        relevant || (group = :unread)
        lower, upper = inputs.box.lower[j], inputs.box.upper[j]
        if row.spec.family_name == :ACSR && family == :ACSR
            if field == :theta
                group, field = :temperature, :temperature
            elseif field == :lay_inner
                group, field = :realized, :inner
            elseif field == :lay_outer
                group, field = :realized, :outer
            elseif field == :d_aluminium
                group, field = :realized, :da
                lower, upper = row.spec.d_aluminium.lower, row.spec.d_aluminium.upper
            elseif field == :d_steel
                group, field = :realized, :ds
                lower, upper = row.spec.d_steel.lower, row.spec.d_steel.upper
            end
        end
        push!(axes, (; name, group, field, lower, upper))
    end
    return axes
end

"""
    effective_value(physical, axis)

Read one effective coordinate from the mapped physical record.
"""
function effective_value(physical, axis)
    axis.group == :temperature && return physical.temperature
    return getproperty(getproperty(physical, axis.group), axis.field)
end

"""
    effective_point(row, physical, axis, value)

Change one effective coordinate without applying a source map. Recompute mass,
area and diameter from the changed construction, with the other inputs fixed.
"""
function effective_point(row, physical, axis, value)
    if axis.group == :temperature
        return merge(physical, (; temperature = value))
    end
    nested = merge(getproperty(physical, axis.group), NamedTuple{(axis.field,)}((value,)))
    changed = merge(physical, NamedTuple{(axis.group,)}((nested,)))
    if row.spec.family_name == :ACSR
        r = changed.realized
        realized = merge(r, construction(row.spec, changed.q, changed.reference,
            r.da, r.ds, r.inner, r.outer))
        changed = merge(changed, (; realized))
    end
    return changed
end

"""
    feasible(row, physical)

Test the effective order constraints and EN 50182 diameter and mass limits.
Numerical tolerance is 32 machine spacings of each aggregate upper limit.
"""
function feasible(row, physical)
    row.spec.family_name == :ACSR || return true
    r = physical.realized
    r.outer <= r.inner || return false
    haskey(row.layout, :structure) && physical.temperature < physical.structure.ambient &&
        return false
    row.spec.d_aluminium.lower <= r.da <= row.spec.d_aluminium.upper || return false
    row.spec.d_steel.lower <= r.ds <= row.spec.d_steel.upper || return false
    for field in (:diameter, :mass)
        interval = getproperty(row.spec, field)
        value = getproperty(r, field)
        tolerance = 32eps(interval.upper)
        interval.lower - tolerance <= value <= interval.upper + tolerance || return false
    end
    return true
end

"""
    effective_interval(row, physical, axis, fraction)

Find a feasible finite-difference interval on the coordinate line through a mapped
point. Bisect a clipped side to its active constraint. A side shorter than one
millionth of the requested step stays at the original point.
"""
function effective_interval(row, physical, axis, fraction)
    value = effective_value(physical, axis)
    step = fraction * (axis.upper - axis.lower)
    bounds = [max(axis.lower, value - step), min(axis.upper, value + step)]
    for k in eachindex(bounds)
        candidate = bounds[k]
        if !feasible(row, effective_point(row, physical, axis, candidate))
            inside, outside = value, candidate
            for _ in 1:50
                trial = (inside + outside) / 2
                if feasible(row, effective_point(row, physical, axis, trial))
                    inside = trial
                else
                    outside = trial
                end
            end
            bounds[k] = abs(inside - value) > METHOD.zero_band * step ? inside : value
        end
    end
    return bounds
end

"""
    sensitivities(inputs, row, points, domain, fraction)

Compute every finite difference in source or effective coordinates. Effective
perturbations stay in the feasible set and never pass through the source maps.
The returned arrays include normalized derivatives and feasible interval widths.
"""
function sensitivities(inputs, row, points, domain, fraction, measure = quantities)
    axes = effective_axes(inputs, row)
    reference = measure(row.problem(nominal_point(inputs.box), 50.0))
    values = zeros(length(reference), length(axes), length(points))
    spans = zeros(length(axes), length(points))
    for (i, p) in enumerate(points)
        physical = physical_inputs(
            row.spec, inputs.families, row.layout, inputs.bindings, p)
        feasible(row, physical) || error("mapped source point is infeasible")
        for (j, axis) in enumerate(axes)
            if domain == :effective
                axis.group == :unread && continue
                lo, hi = effective_interval(row, physical, axis, fraction)
                spans[j, i] = hi - lo
                hi > lo || continue
                a = measure(effective_problem(row.spec,
                    effective_point(row, physical, axis, lo), row.layout, inputs.mat, 50.0))
                b = measure(effective_problem(row.spec,
                    effective_point(row, physical, axis, hi), row.layout, inputs.mat, 50.0))
                width = axis.upper - axis.lower
            else
                index = findfirst(==(axis.name), inputs.box.names)
                width = inputs.box.upper[index] - inputs.box.lower[index]
                lo = max(inputs.box.lower[index], p[axis.name] - fraction * width)
                hi = min(inputs.box.upper[index], p[axis.name] + fraction * width)
                spans[j, i] = hi - lo
                hi > lo || continue
                a = measure(row.problem(merge(p, NamedTuple{(axis.name,)}((lo,))), 50.0))
                b = measure(row.problem(merge(p, NamedTuple{(axis.name,)}((hi,))), 50.0))
            end
            for k in eachindex(reference)
                nominal = last(reference[k])
                isfinite(nominal) || continue
                scale = iszero(nominal) ? 1.0 : abs(nominal)
                values[k, j, i] = (last(b[k]) - last(a[k])) / (hi - lo) * width / scale
            end
        end
    end
    return (; reference, axes, values, spans)
end

"""
    dependence(inputs, destination)

Analyze every quantity at all 303 points in source and effective coordinates.
Write complete dependence and sign-change tables to the chosen evidence directory.
Classify a source sign change as map-induced when its effective derivative has no
sign change. A sign change in effective coordinates is intrinsic. Build no hull
when an intrinsic sign change remains. The numerical zero band is 0.000001.
"""
function dependence(inputs, destination)
    points = source_points(inputs.box)
    mkpath(destination)
    results = NamedTuple[]
    open(joinpath(destination, "equivalent-ranges.csv"), "w") do io
        println(io, "circuit,quantity,nominal,sampled_minimum,sampled_maximum")
        for row in inputs.definitions
            observed = [quantities(row.problem(p, 50.0)) for p in points]
            for k in eachindex(first(observed))
                values = [last(q[k]) for q in observed]
                println(io,
                    join(
                        (row.name, first(first(observed)[k]), first(values),
                            minimum(values), maximum(values)),
                        ','))
            end
        end
    end
    open(joinpath(destination, "equivalent-sensitivities.csv"), "w") do io
        println(io,
            "circuit,quantity,coordinate,domain,minimum,maximum,negative_sample,positive_sample,blocked_samples")
        for row in inputs.definitions
            println("Analyzing ", row.name)
            flush(stdout)
            source = sensitivities(inputs, row, points, :source, METHOD.difference_step)
            physical = sensitivities(
                inputs, row, points, :effective, METHOD.difference_step)
            for k in eachindex(source.reference), (j, axis) in enumerate(source.axes)

                for (domain, table) in ((:source, source), (:effective, physical))
                    derivatives = table.values[k, j, :]
                    minimum_value, negative = findmin(derivatives)
                    maximum_value, positive = findmax(derivatives)
                    println(io,
                        join(
                            (row.name, first(source.reference[k]), axis.name,
                                domain, minimum_value, maximum_value, negative, positive,
                                count(iszero, table.spans[j, :])),
                            ','))
                end
                source_flip = minimum(source.values[k, j, :]) < -METHOD.zero_band &&
                              maximum(source.values[k, j, :]) > METHOD.zero_band
                effective_flip = minimum(physical.values[k, j, :]) < -METHOD.zero_band &&
                                 maximum(physical.values[k, j, :]) > METHOD.zero_band
                if source_flip || effective_flip
                    derivatives = physical.values[k, j, :]
                    low, high = argmin(derivatives), argmax(derivatives)
                    push!(results,
                        (; circuit = row.name, quantity = first(source.reference[k]),
                            coordinate = axis.name, classification = effective_flip ?
                                                                     "intrinsic" :
                                                                     "map-induced",
                            source_minimum = minimum(source.values[k, j, :]),
                            source_maximum = maximum(source.values[k, j, :]),
                            effective_minimum = minimum(derivatives), effective_maximum = maximum(derivatives),
                            negative_sample = low, positive_sample = high))
                end
            end
            flush(io)
            println("Finished ", row.name)
            flush(stdout)
        end
    end
    open(joinpath(destination, "equivalent-signs.csv"), "w") do io
        println(io,
            "circuit,quantity,coordinate,classification,source_minimum,source_maximum,effective_minimum,effective_maximum,negative_sample,positive_sample")
        for row in results
            println(io, join(values(row), ','))
        end
    end
    open(joinpath(destination, "equivalent-stop.txt"), "w") do io
        intrinsic = filter(r -> r.classification == "intrinsic", results)
        println(io, "Completed every circuit and sample in source and effective coordinates.")
        println(io, "Map-induced changes: ", count(r -> r.classification == "map-induced", results))
        println(io, "Intrinsic changes: ", length(intrinsic))
        for row in intrinsic
            println(io, row.circuit, " ", row.quantity, " ", row.coordinate,
                " ", row.effective_minimum, " ", row.effective_maximum)
        end
        isempty(intrinsic) ||
            println(io, "STOP after full analysis. Hulls and proposal await review.")
    end
    return results
end

"""
    observations(problem, kind)

Combine the coaxial inventory with the equivalent coordinates. Raw keys describe
LCM inputs. Equivalent keys describe the reversible radial parameterization.
"""
function observations(problem, kind)
    return ["raw." * first(p) => last(p) for p in quantities(problem)] ∪
           ["eq." * first(p) => last(p) for p in equivalent_values(problem, kind)]
end

"""
    derivative(inputs, row, physical, axis, metric, fraction)

Evaluate one normalized effective sensitivity while keeping the perturbation
feasible. A blocked direction returns zero and contributes no sign evidence.
"""
function derivative(inputs, row, physical, axis, metric, fraction)
    lo, hi = effective_interval(row, physical, axis, fraction)
    hi > lo || return 0.0
    a = observations(
        effective_problem(row.spec, effective_point(row, physical, axis, lo),
            row.layout, inputs.mat, 50.0),
        row.layout.kind)[metric].second
    b = observations(
        effective_problem(row.spec, effective_point(row, physical, axis, hi),
            row.layout, inputs.mat, 50.0),
        row.layout.kind)[metric].second
    return (b - a) / (hi - lo) * (axis.upper - axis.lower)
end

"""
    sign_conflicts(gradients, indices, drivers)

Return sensitivity columns with positive and negative samples in a cell.
"""
function sign_conflicts(gradients, indices, drivers)
    return [j
            for j in drivers
            if minimum(gradients[j, indices]) < -METHOD.zero_band &&
        maximum(gradients[j, indices]) > METHOD.zero_band]
end

"""
    split_cell(inputs, row, data, metric, indices, lower, upper, depth)

Choose the driving coordinate and threshold that minimize disagreeing derivative
pairs in the children. Refine a bracketed derivative zero on that coordinate.
Choose splits from the measured sensitivities, up to `METHOD.subdivision_depth`.
"""
function split_cell(inputs, row, data, metric, indices, lower, upper, depth)
    gradients = data.effective.values[metric, :, :]
    drivers = findall(j -> maximum(abs, gradients[j, :]) > METHOD.zero_band, eachindex(data.axes))
    conflicts = sign_conflicts(gradients, indices, drivers)
    isempty(conflicts) && return [(; lower, upper, depth, indices, consistent = true)]
    depth == METHOD.subdivision_depth &&
        return [(; lower, upper, depth, indices, consistent = false)]
    best = (score = Inf, axis = 0, threshold = 0.0, left = Int[], right = Int[])
    for j in drivers
        positions = sort(unique(data.positions[j, indices]))
        for n in 1:(length(positions) - 1)
            threshold = (positions[n] + positions[n + 1]) / 2
            left = filter(i -> data.positions[j, i] <= threshold, indices)
            right = filter(i -> data.positions[j, i] > threshold, indices)
            score = sum(count(i -> gradients[k, i] < -METHOD.zero_band, part) *
                        count(i -> gradients[k, i] > METHOD.zero_band, part)
            for k in conflicts
            for part in (left, right))
            if score < best.score
                best = (; score, axis = j, threshold, left, right)
            end
        end
    end
    best.axis == 0 && return [(;
        lower, upper, depth = METHOD.subdivision_depth, indices, consistent = false)]
    j = best.axis
    threshold = best.threshold
    for k in conflicts
        first_sign = sum(sign(gradients[k, i]) for i in best.left)
        second_sign = sum(sign(gradients[k, i]) for i in best.right)
        first_sign * second_sign < 0 || continue
        left_witnesses = filter(
            i -> gradients[k, i] * first_sign > 0 &&
                 abs(gradients[k, i]) > METHOD.zero_band, best.left)
        right_witnesses = filter(
            i -> gradients[k, i] * second_sign > 0 &&
                 abs(gradients[k, i]) > METHOD.zero_band,
            best.right)
        if isempty(left_witnesses) || isempty(right_witnesses)
            continue
        end
        a = left_witnesses[argmin(abs.(data.positions[j, left_witnesses] .- threshold))]
        b = right_witnesses[argmin(abs.(data.positions[j, right_witnesses] .- threshold))]
        anchor = data.physical[a]
        lo, hi = data.positions[j, a], data.positions[j, b]
        left_state = effective_point(row, anchor, data.axes[j], lo)
        right_state = effective_point(row, anchor, data.axes[j], hi)
        feasible(row, left_state) && feasible(row, right_state) || continue
        dl = derivative(
            inputs, row, left_state, data.axes[k], metric, METHOD.difference_step)
        dh = derivative(
            inputs, row, right_state, data.axes[k], metric, METHOD.difference_step)
        dl * dh < 0 || continue
        for _ in 1:40
            middle = (lo + hi) / 2
            dm = derivative(
                inputs, row, effective_point(row, anchor, data.axes[j], middle),
                data.axes[k], metric, METHOD.difference_step)
            if iszero(dm)
                lo = hi = middle
                break
            elseif dl * dm > 0
                lo, dl = middle, dm
            else
                hi = middle
            end
            hi - lo <= 1e-9 * (upper[j] - lower[j]) && break
        end
        threshold = (lo + hi) / 2
        break
    end
    left = filter(i -> data.positions[j, i] <= threshold, indices)
    right = filter(i -> data.positions[j, i] >= threshold, indices)
    if isempty(left) || isempty(right) || threshold <= lower[j] || threshold >= upper[j]
        return [(;
            lower, upper, depth = METHOD.subdivision_depth, indices, consistent = false)]
    end
    lower_right, upper_left = copy(lower), copy(upper)
    lower_right[j], upper_left[j] = threshold, threshold
    return [split_cell(inputs, row, data, metric, left, lower, upper_left, depth + 1);
            split_cell(inputs, row, data, metric, right, lower_right, upper, depth + 1)]
end

"""
    constraint_groups(row, axes)

Describe the effective order and aggregate constraints for the construction.
The numerical search closes its driving set over these coupled coordinates.
"""
function constraint_groups(row, axes)
    row.spec.family_name == :ACSR || return Vector{Int}[]
    metal = findall(a -> a.group == :realized, axes)
    thermal = findall(
        a -> a.group == :temperature ||
             (a.group == :structure && a.field == :ambient), axes)
    return [metal, thermal]
end

function constraint_values(row, physical)
    row.spec.family_name == :ACSR || return Float64[]
    r = physical.realized
    thermal = haskey(row.layout, :structure) ?
              physical.temperature - physical.structure.ambient : 1.0
    return [r.inner - r.outer, thermal,
        r.diameter - row.spec.diameter.lower, row.spec.diameter.upper - r.diameter,
        r.mass - row.spec.mass.lower, row.spec.mass.upper - r.mass]
end

"""
    cell_interval(row, state, axis, lower, upper)

Intersect a coordinate line with a cell and the effective feasible set.
Bisection approaches an active constraint from its feasible side.
"""
function cell_interval(row, state, axis, lower, upper)
    current = effective_value(state, axis)
    result = [lower, upper]
    for k in 1:2
        if !feasible(row, effective_point(row, state, axis, result[k]))
            inside, outside = current, result[k]
            for _ in 1:52
                trial = (inside + outside) / 2
                if feasible(row, effective_point(row, state, axis, trial))
                    inside = trial
                else
                    outside = trial
                end
            end
            result[k] = inside
        end
    end
    return result
end

"""
    golden(function_value, lower, upper)

Minimize a scalar function on an interval with 32 golden-section steps.
"""
function golden(function_value, lower, upper)
    ratio = (sqrt(5.0) - 1) / 2
    a, b = lower, upper
    c, d = b - ratio * (b - a), a + ratio * (b - a)
    fc, fd = function_value(c), function_value(d)
    for _ in 1:32
        if fc <= fd
            b, d, fd = d, c, fc
            c = b - ratio * (b - a)
            fc = function_value(c)
        else
            a, c, fc = c, d, fd
            d = a + ratio * (b - a)
            fd = function_value(d)
        end
    end
    candidates = [
        (lower, function_value(lower)), (upper, function_value(upper)), (c, fc), (d, fd)]
    return candidates[argmin(last.(candidates))]
end

function metric_value(inputs, row, state, metric)
    return observations(effective_problem(row.spec, state, row.layout, inputs.mat, 50.0),
        row.layout.kind)[metric].second
end

"""
    face_point(row, state, axes, first_axis, second_axis, value, face, cell)

Move on an active constraint face by setting one coordinate and solving for a
second coordinate by bisection. Return a feasibility flag with the resulting state.
"""
function face_point(row, state, axes, first_axis, second_axis, value, face, cell)
    changed = effective_point(row, state, axes[first_axis], value)
    lo, hi = cell.lower[second_axis], cell.upper[second_axis]
    a = effective_point(row, changed, axes[second_axis], lo)
    b = effective_point(row, changed, axes[second_axis], hi)
    fa, fb = constraint_values(row, a)[face], constraint_values(row, b)[face]
    fa * fb <= 0 && fa != fb || return false, state
    for _ in 1:48
        middle = (lo + hi) / 2
        trial = effective_point(row, changed, axes[second_axis], middle)
        fm = constraint_values(row, trial)[face]
        if fa * fm > 0
            lo, fa = middle, fm
        else
            hi = middle
        end
    end
    for value in (lo, hi, (lo + hi) / 2)
        candidate = effective_point(row, changed, axes[second_axis], value)
        feasible(row, candidate) && return true, candidate
    end
    return false, state
end

"""
    refine_cell(inputs, row, data, metric, cell, initial, orientation, drivers)

Refine a feasible candidate along coordinate lines and along active faces.
One-dimensional face searches solve the companion coordinate by bisection.
"""
function refine_cell(
        inputs, row, data, metric, cell, initial, orientation, drivers, passes = 4)
    state = initial
    best = orientation * metric_value(inputs, row, state, metric)
    for _ in 1:passes
        previous = best
        for j in drivers
            lo, hi = cell_interval(row, state, data.axes[j], cell.lower[j], cell.upper[j])
            hi > lo || continue
            value, objective = golden(
                x -> orientation * metric_value(inputs, row,
                    effective_point(row, state, data.axes[j], x), metric),
                lo,
                hi)
            if objective < best
                state = effective_point(row, state, data.axes[j], value)
                best = objective
            end
        end
        margins = constraint_values(row, state)
        for face in eachindex(margins)
            abs(margins[face]) <= 1e-8 * max(1.0, abs(margins[face])) || continue
            coupled = [j
                       for j in drivers
                       if abs(
                constraint_values(row, effective_point(row, state, data.axes[j], cell.upper[j]))[face] -
                constraint_values(row, effective_point(row, state, data.axes[j], cell.lower[j]))[face]) >
                          1e-12]
            for (n, j) in enumerate(coupled), k in coupled[(n + 1):end]

                grid = range(cell.lower[j], cell.upper[j]; length = 9)
                feasible_values = [x
                                   for x in grid
                                   if first(face_point(row, state,
                    data.axes, j, k, x, face, cell))]
                length(feasible_values) >= 2 || continue
                lo, hi = extrema(feasible_values)
                value, objective = golden(
                    x -> begin
                        valid, trial = face_point(row, state, data.axes, j, k, x, face, cell)
                        valid ? orientation * metric_value(inputs, row, trial, metric) :
                        Inf
                    end,
                    lo,
                    hi)
                if objective < best
                    valid, trial = face_point(
                        row, state, data.axes, j, k, value, face, cell)
                    valid && (state = trial; best = objective)
                end
            end
        end
        previous - best <= 1e-11 * max(abs(best), eps()) && break
    end
    return state, orientation * best
end

"""
    cell_hull(inputs, row, data, metric, cell)

Enumerate driving-coordinate corners, order ties and radial intersections with
active aggregate limits. Refine the best feasible candidates on their faces.
"""
function cell_hull(inputs, row, data, metric, cell)
    drivers = findall(
        j -> maximum(abs, data.effective.values[metric, j, :]) > METHOD.zero_band,
        eachindex(data.axes))
    constrained = Int[]
    for group in constraint_groups(row, data.axes)
        isempty(intersect(group, drivers)) || (constrained = union(constrained, group))
    end
    search_axes = sort(constrained)
    states = [data.physical[i] for i in cell.indices]
    isempty(states) && error("cell has no feasible sample")
    extrema_states = NamedTuple[]
    extrema_values = Float64[]
    total_candidates = 0
    for orientation in (1, -1)
        signs = [sum(data.effective.values[metric, j, cell.indices]) >= 0 ? 1 : -1
                 for j in eachindex(data.axes)]
        anchor = first(states)
        # Independent monotone axes take their oriented endpoint before enumeration.
        for j in setdiff(drivers, constrained)
            value = signs[j] * orientation > 0 ? cell.lower[j] : cell.upper[j]
            anchor = effective_point(row, anchor, data.axes[j], value)
        end
        candidates = typeof(anchor)[]
        for mask in 0:(2 ^ length(constrained) - 1)
            corner = anchor
            for (n, j) in enumerate(constrained)
                value = iszero(mask & (1 << (n - 1))) ? cell.lower[j] : cell.upper[j]
                corner = effective_point(row, corner, data.axes[j], value)
            end
            if feasible(row, corner)
                push!(candidates, corner)
            else
                inside, outside = 0.0, 1.0
                for _ in 1:48
                    fraction = (inside + outside) / 2
                    trial = anchor
                    for j in constrained
                        value = (1 - fraction) * effective_value(anchor, data.axes[j]) +
                                fraction * effective_value(corner, data.axes[j])
                        trial = effective_point(row, trial, data.axes[j], value)
                    end
                    if feasible(row, trial)
                        inside = fraction
                    else
                        outside = fraction
                    end
                end
                trial = anchor
                for j in constrained
                    value = (1 - inside) * effective_value(anchor, data.axes[j]) +
                            inside * effective_value(corner, data.axes[j])
                    trial = effective_point(row, trial, data.axes[j], value)
                end
                push!(candidates, trial)
            end
        end
        coordinates = [[effective_value(state, data.axes[j]) for j in drivers]
                       for state in candidates]
        keep = [i
                for i in eachindex(candidates)
                if !any(
            l -> l != i &&
                 all(
                     n -> signs[drivers[n]] * orientation *
                          (coordinates[l][n] - coordinates[i][n]) <= 0,
                     eachindex(drivers)) &&
                 any(
                     n -> signs[drivers[n]] * orientation *
                          (coordinates[l][n] - coordinates[i][n]) < 0,
                     eachindex(drivers)),
            eachindex(candidates))]
        candidates = unique(candidates[keep])
        total_candidates += length(candidates)
        best_state = anchor
        best = orientation * metric_value(inputs, row, anchor, metric)
        for candidate in candidates
            refined, value = refine_cell(inputs, row, data, metric, cell, candidate,
                orientation, search_axes, 1)
            if orientation * value < best
                best_state, best = refined, orientation * value
            end
        end
        refined, value = refine_cell(inputs, row, data, metric, cell, best_state,
            orientation, search_axes)
        push!(extrema_states, refined)
        push!(extrema_values, value)
    end
    low, high = extrema_values
    allowance = 64eps(max(abs(low), abs(high)))
    return (; lower = prevfloat(low - allowance), upper = nextfloat(high + allowance),
        low_state = first(extrema_states), high_state = last(extrema_states),
        candidates = total_candidates)
end

"""
    reduction_data(inputs, row)

Measure the coaxial inputs and equivalent coordinates at the source samples,
with their source and effective sensitivities. Keep physical samples for searches.
"""
function reduction_data(inputs, row)
    points = source_points(inputs.box)
    source = sensitivities(inputs, row, points, :source, METHOD.difference_step,
        p -> observations(p, row.layout.kind))
    effective = sensitivities(inputs, row, points, :effective, METHOD.difference_step,
        p -> observations(p, row.layout.kind))
    physical = [physical_inputs(row.spec, inputs.families, row.layout, inputs.bindings, p)
                for p in points]
    positions = [a.group == :unread ? 0.0 : effective_value(state, a)
                 for a in effective.axes, state in physical]
    observed = reduce(hcat, [last.(observations(row.problem(p, 50.0), row.layout.kind))
                             for p in points])
    return (;
        points, source, effective, physical, positions, observed, axes = effective.axes)
end

"""
    metric_classes(data)

Group affine aliases within a circuit after checking their normalized sample
profiles and driving sets. A class shares one search, with measured alias error.
"""
function metric_classes(data)
    classes = NamedTuple[]
    for k in eachindex(data.source.reference)
        drivers = findall(
            j -> maximum(abs, data.effective.values[k, j, :]) > METHOD.zero_band,
            eachindex(data.axes))
        isempty(drivers) && continue
        lo, hi = extrema(data.observed[k, :])
        hi > lo || continue
        profile = (data.observed[k, :] .- lo) ./ (hi - lo)
        assigned = false
        for group in classes
            group.drivers == drivers || continue
            error = maximum(abs.(profile .- group.profile))
            if error <= 1e-11
                push!(group.members, k)
                push!(group.errors, error)
                assigned = true
                break
            end
        end
        assigned || push!(classes, (;
            representative = k, members = [k], errors = [0.0], drivers, profile))
    end
    return classes
end

"""
    quantity_hulls(inputs, row, data)

Subdivide inconsistent scalar quantities and search the feasible cells. Affine
aliases inherit the searched interval with their measured numerical allowance.
Fixed quantities receive an outward roundoff envelope and no independent axis.
"""
function quantity_hulls(inputs, row, data)
    n = length(data.source.reference)
    bounds = [(; lower = minimum(data.observed[k, :]),
                  upper = maximum(data.observed[k, :])) for k in 1:n]
    cells = [NamedTuple[] for _ in 1:n]
    methods = fill("fixed", n)
    lower = [a.lower for a in data.axes]
    upper = [a.upper for a in data.axes]
    for group in metric_classes(data)
        k = group.representative
        partition = split_cell(
            inputs, row, data, k, collect(eachindex(data.points)), lower, upper, 0)
        all(c -> c.consistent, partition) ||
            error("depth six remains inconsistent: $(row.name) $(first(data.source.reference[k]))")
        searched = [cell_hull(inputs, row, data, k, cell) for cell in partition]
        lo = minimum(h -> h.lower, searched)
        hi = maximum(h -> h.upper, searched)
        reference_min, reference_max = extrema(data.observed[k, :])
        for (index, member) in enumerate(group.members)
            member_min, member_max = extrema(data.observed[member, :])
            ratio = (member_max - member_min) / (reference_max - reference_min)
            allowance = group.errors[index] * (member_max - member_min)
            a = member_min + ratio * (lo - reference_min) - allowance
            b = member_min + ratio * (hi - reference_min) + allowance
            rounding = 64eps(max(abs(a), abs(b)))
            bounds[member] = (;
                lower = prevfloat(a - rounding), upper = nextfloat(b + rounding))
            cells[member] = [merge(cell,
                                 (
                                     value_lower = prevfloat(member_min +
                                                             ratio *
                                                             (bound.lower - reference_min) -
                                                             allowance - rounding),
                                     value_upper = nextfloat(member_min +
                                                             ratio *
                                                             (bound.upper - reference_min) +
                                                             allowance + rounding),
                                     candidates = bound.candidates))
                             for (cell, bound) in zip(partition, searched)]
            methods[member] = member == k ?
                              "corners, projected aggregate faces, golden face search" :
                              "verified affine alias"
        end
        println("  ", first(data.source.reference[k]), " ",
            length(partition), " cells ", lo, " ", hi)
        flush(stdout)
    end
    for k in 1:n
        methods[k] == "fixed" || continue
        nominal = last(data.source.reference[k])
        isfinite(nominal) || continue
        rounding = 64eps(max(1.0, abs(nominal)))
        bounds[k] = (; lower = prevfloat(bounds[k].lower - rounding),
            upper = nextfloat(bounds[k].upper + rounding))
        cells[k] = [(; lower, upper, depth = 0,
            indices = collect(eachindex(data.points)), consistent = true,
            value_lower = bounds[k].lower, value_upper = bounds[k].upper, candidates = 0)]
    end
    return (; bounds, cells, methods)
end

"""
    containment(inputs, row, hulls, points)

Check mapped source observations against their computed hulls, with no comparison
tolerance. Return escaped indices and the largest normalized outward margin.
A nonpositive outward margin means containment.
"""
function containment(inputs, row, hulls, points)
    escaped = Tuple{Int, Int}[]
    largest = -Inf
    observed = Matrix{Float64}(undef, length(hulls.bounds), length(points))
    for (i, p) in enumerate(points)
        values = last.(observations(row.problem(p, 50.0), row.layout.kind))
        observed[:, i] = values
        physical = physical_inputs(
            row.spec, inputs.families, row.layout, inputs.bindings, p)
        axes = effective_axes(inputs, row)
        positions = [a.group == :unread ? 0.0 : effective_value(physical, a) for a in axes]
        for (k, (value, interval)) in enumerate(zip(values, hulls.bounds))
            isfinite(value) || continue
            matching = [cell
                        for cell in hulls.cells[k]
                        if all(
                j -> axes[j].group == :unread ||
                     cell.lower[j] <= positions[j] <= cell.upper[j],
                eachindex(axes))]
            scale = max(abs(interval.lower), abs(interval.upper), eps())
            margin = minimum(
                (max(cell.value_lower - value, value - cell.value_upper) / scale
                for cell in matching);
                init = Inf)
            largest = max(largest, margin)
            margin <= 0 || push!(escaped, (k, i))
        end
    end
    return (; escaped, largest, observed)
end

"""
    hull_study(inputs, destination)

Compute the cells and hulls of the full dataset. Verify
`METHOD.verification_samples` independent mapped Sobol points, and write cell
and bound evidence to the requested directory.
"""
function hull_study(inputs, destination)
    mkpath(destination)
    verification = [latent(inputs.box, p)
                    for p in background_points(inputs.box, METHOD.verification_samples; seed = 3)]
    results = NamedTuple[]
    for row in inputs.definitions
        println("Hull study: ", row.name)
        flush(stdout)
        data = reduction_data(inputs, row)
        hulls = quantity_hulls(inputs, row, data)
        checked = containment(inputs, row, hulls, verification)
        push!(results, (; row, data, hulls, checked))
        println("Escapes: ", length(checked.escaped), "; largest margin: ", checked.largest)
        flush(stdout)
    end
    open(joinpath(destination, "equivalent-hulls.csv"), "w") do io
        println(io, "circuit,quantity,nominal,lower,upper,cells,method,escapes")
        for result in results, (k, pair) in enumerate(result.data.source.reference)

            interval = result.hulls.bounds[k]
            println(io,
                join(
                    (result.row.name, first(pair), last(pair),
                        interval.lower, interval.upper,
                        length(result.hulls.cells[k]), replace(result.hulls.methods[k], ',' =>
                            ';'),
                        count(p -> first(p) == k, result.checked.escaped)),
                    ','))
        end
    end
    open(joinpath(destination, "equivalent-cells.csv"), "w") do io
        println(io,
            "circuit,quantity,cell,depth,coordinate,lower,upper,samples,consistent,value_lower,value_upper,candidates")
        for result in results, (k, cells) in enumerate(result.hulls.cells),
            (i, cell) in enumerate(cells)
            for (j, axis) in enumerate(result.data.axes)
                axis.group == :unread && continue
                println(io,
                    join(
                        (result.row.name,
                            first(result.data.source.reference[k]), i, cell.depth,
                            axis.name, cell.lower[j], cell.upper[j],
                            length(cell.indices), cell.consistent,
                            cell.value_lower, cell.value_upper, cell.candidates),
                        ','))
            end
        end
    end
    open(joinpath(destination, "equivalent-containment.csv"), "w") do io
        println(io, "circuit,samples,escapes,largest_outward_margin")
        for result in results
            println(io,
                join(
                    (result.row.name, length(verification),
                        length(result.checked.escaped), result.checked.largest),
                    ','))
        end
    end
    return (; inputs, results, verification)
end

"""
    normalized(value, nominal, lower, upper)

Map an equivalent value to a signed position with nominal at zero and hull
endpoints at minus and plus one. Ignore differences of at most 64 machine spacings.
"""
function normalized(value, nominal, lower, upper)
    abs(value - nominal) <= 64eps(max(abs(value), abs(nominal))) && return 0.0
    return value < nominal ? (value - nominal) / (nominal - lower) :
           (value - nominal) / (upper - nominal)
end

function equivalent_value(position, member)
    return member.nominal +
           position * (position < 0 ? member.nominal - member.lower :
            member.upper - member.nominal)
end

function metric_family(row, key)
    startswith(key, "earth.") && return :EARTH
    chunks = split(key, '.')
    if first(chunks) == "cable" && haskey(row.layout, :shield_spec) &&
       parse(Int, chunks[2]) > 3
        return row.layout.shield_spec.family_name
    end
    return row.spec.family_name
end

"""
    equivalent_coordinates(study)

Group equivalent quantities only when their signed normalized positions agree
at every analysis and verification point within `METHOD.sharing_tolerance`.
Complete-link grouping
checks a new member against each existing member. Source dependencies follow the
computed source sensitivities, including the construction maps.
"""
function equivalent_coordinates(study)
    groups = NamedTuple[]
    for (ri, result) in enumerate(study.results),
        (k, pair) in enumerate(result.data.source.reference)

        startswith(first(pair), "eq.") || continue
        result.hulls.methods[k] == "fixed" && continue
        key = first(pair)[4:end]
        interval = result.hulls.bounds[k]
        nominal = last(pair)
        values = [result.data.observed[k, :]; result.checked.observed[k, :]]
        profile = [normalized(value, nominal, interval.lower, interval.upper)
                   for value in values]
        drivers = [study.inputs.box.names[j]
                   for j in eachindex(study.inputs.box.names)
                   if maximum(abs, result.data.source.values[k, j, :]) > METHOD.zero_band]
        member = (; row = ri, metric = k, key, nominal, lower = interval.lower,
            upper = interval.upper, profile, drivers)
        family = metric_family(result.row, key)
        quantity = last(split(key, '.'))
        added = false
        for group in groups
            group.family == family && group.quantity == quantity || continue
            error = maximum(maximum(abs.(profile .- m.profile)) for m in group.members)
            if error <= METHOD.sharing_tolerance
                push!(group.members, member)
                group.error[1] = max(group.error[1], error)
                added = true
                break
            end
        end
        if !added
            index = 1 + count(g -> g.family == family && g.quantity == quantity, groups)
            name = Symbol(quantity, "_", family, "_", index)
            push!(groups, (; name, family, quantity, members = [member], error = [0.0]))
        end
    end
    return groups
end

"""
    coordinate_problem(inputs, row, group, row_index, point, frequency)

Change one equivalent coordinate group from nominal through the public plain-tube
builder. The HPM point has one normalized field, `q`.
"""
function coordinate_problem(inputs, row, group, row_index, point, frequency)
    source = row.problem(nominal_point(inputs.box), frequency)
    values = Dict(equivalent_values(source, row.layout.kind))
    for member in group.members
        member.row == row_index || continue
        values[member.key] = equivalent_value(point.q, member)
    end
    return equivalent(source, values, row.layout.kind)
end

"""
    equivalent_problem(inputs, row_index, groups, point, frequency)

Build a circuit at a point of the proposed equivalent box. A missing coordinate
keeps its nominal value. The caller supplies signed normalized positions.
"""
function equivalent_problem(inputs, row_index, groups, point, frequency)
    row = inputs.definitions[row_index]
    source = row.problem(nominal_point(inputs.box), frequency)
    values = Dict(equivalent_values(source, row.layout.kind))
    for group in groups
        position = get(point, group.name, 0.0)
        for member in group.members
            member.row == row_index || continue
            values[member.key] = equivalent_value(position, member)
        end
    end
    return equivalent(source, values, row.layout.kind)
end

"""
    admission(study, groups)

Apply the sweep measure at each equivalent hull endpoint with the other coordinates
nominal. Ask the extension for all sequences at 50 and 650 Hz. Admit at 0.1 percent.
"""
function admission(study, groups)
    decisions = NamedTuple[]
    formulation = CoaxialFormulation(transposition = true)
    box = ParameterBox(:q => (0.0, -1.0, 1.0))
    for group in groups
        rows = unique(m.row for m in group.members)
        branches = Dict(ri => (
                            problem = (p, f) -> coordinate_problem(study.inputs,
                                study.results[ri].row, group, ri, p, f),
                            len = 1.0,
                            zbase = 1.0) for ri in rows)
        model = ExactPiModel(branches, formulation, box)
        best = (change = -1.0, circuit = "", sequence = :zero,
            quantity = :z, frequency = 50.0, endpoint = 0.0)
        for ri in rows, frequency in (50.0, 650.0), endpoint in (-1.0, 1.0)
            nominal = model.branches[ri].zy((q = 0.0,), frequency)
            result = model.branches[ri].zy((q = endpoint,), frequency)
            for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

                relative = abs(result[sequence][quantity] - nominal[sequence][quantity]) /
                           abs(nominal[sequence][quantity])
                if relative > best.change
                    best = (; change = relative, circuit = study.results[ri].row.name,
                        sequence, quantity, frequency, endpoint)
                end
            end
        end
        push!(decisions,
            merge(
                (;
                    name = group.name, admitted = best.change >= METHOD.admission_threshold), best))
        println("Coordinate ", group.name, " effect ", best.change)
        flush(stdout)
    end
    return decisions
end

function source_positions(box, point)
    return [box.scale[j] == :log ?
            (2log(point[key]) - log(box.lower[j]) - log(box.upper[j])) /
            (log(box.upper[j]) - log(box.lower[j])) :
            (2point[key] - box.lower[j] - box.upper[j]) / (box.upper[j] - box.lower[j])
            for (j, key) in enumerate(box.names)]
end

function pair_distance(
        study, first_group, second_group, first_target, second_target, position)
    point = latent(study.inputs.box, position)
    a, b = first(first_group.members), first(second_group.members)
    row_a, row_b = study.results[a.row].row, study.results[b.row].row
    values_a = Dict(equivalent_values(row_a.problem(point, 50.0), row_a.layout.kind))
    values_b = a.row == b.row ? values_a :
               Dict(equivalent_values(row_b.problem(point, 50.0), row_b.layout.kind))
    x = normalized(values_a[a.key], a.nominal, a.lower, a.upper)
    y = normalized(values_b[b.key], b.nominal, b.lower, b.upper)
    return max(abs(x - first_target), abs(y - second_target))
end

"""
    pair_extreme(study, a, b, x, y)

Refine the closest source point to one proposed pair of equivalent endpoints.
"""
function pair_extreme(study, a, b, x, y)
    pa, pb = first(a.members).profile, first(b.members).profile
    first_drivers = unique([d for m in a.members for d in m.drivers])
    second_drivers = unique([d for m in b.members for d in m.drivers])
    errors = max.(abs.(pa .- x), abs.(pb .- y))
    index = argmin(errors)
    point = index <= length(study.results[1].data.points) ?
            study.results[1].data.points[index] :
            study.verification[index - length(study.results[1].data.points)]
    position = source_positions(study.inputs.box, point)
    drivers = findall(k -> k in union(first_drivers, second_drivers), study.inputs.box.names)
    objective = pair_distance(study, a, b, x, y, position)
    for _ in 1:3
        previous = objective
        for j in drivers
            grid = collect(range(-1.0, 1.0; length = 9))
            values = [pair_distance(study, a, b, x, y,
                          [i == j ? value : position[i] for i in eachindex(position)])
                      for value in grid]
            best = argmin(values)
            lo, hi = grid[max(1, best - 1)], grid[min(length(grid), best + 1)]
            value, distance = golden(
                t -> pair_distance(study, a, b, x, y,
                    [i == j ? t : position[i] for i in eachindex(position)]),
                lo,
                hi)
            if distance < objective
                position[j] = value
                objective = distance
            end
        end
        previous - objective < 1e-8 && break
    end
    return (; distance = objective, position)
end

"""
    decorrelations(study, groups, decisions)

Report pairs of admitted coordinates with common source inputs. Select the
extreme combination farthest from the sampled source image, then refine the
closest source point by low-dimensional coordinate searches. Try another corner
when refinement reaches the first one. Record joint residuals and endpoint values.
"""
function decorrelations(study, groups, decisions)
    admitted = findall(d -> d.admitted, decisions)
    pairs = NamedTuple[]
    for (n, first_index) in enumerate(admitted), second_index in admitted[(n + 1):end]

        a, b = groups[first_index], groups[second_index]
        first_drivers = unique([d for m in a.members for d in m.drivers])
        second_drivers = unique([d for m in b.members for d in m.drivers])
        shared = intersect(first_drivers, second_drivers)
        isempty(shared) && continue
        pa, pb = first(a.members).profile, first(b.members).profile
        alternatives = [(x, y, minimum(max.(abs.(pa .- x), abs.(pb .- y))))
                        for x in (-1.0, 1.0) for y in (-1.0, 1.0)]
        sort!(alternatives; by = last, rev = true)
        x, y, _ = first(alternatives)
        refined = pair_extreme(study, a, b, x, y)
        for alternative in alternatives[2:end]
            refined.distance > METHOD.zero_band && break
            trial_x, trial_y, _ = alternative
            trial = pair_extreme(study, a, b, trial_x, trial_y)
            if trial.distance > refined.distance
                x, y, refined = trial_x, trial_y, trial
            end
        end
        refined.distance > METHOD.zero_band ||
            error("no excluded corner established for $(a.name) and $(b.name)")
        objective, position = refined.distance, refined.position
        first_member, second_member = first(a.members), first(b.members)
        push!(pairs,
            (; first = a.name, second = b.name, shared = join(string.(shared), ';'),
                first_endpoint = x, second_endpoint = y,
                first_value = equivalent_value(x, first_member), second_value = equivalent_value(
                    y, second_member),
                residual = objective, source = string(latent(study.inputs.box, position))))
        println("Pair ", a.name, " / ", b.name, " residual ", objective)
        flush(stdout)
    end
    return pairs
end

"""
    proposal(study, destination)

Write the coordinate definitions, sharing measurements, admission results and
source correlations, while leaving the reduced case to the review step.
"""
function proposal(study, destination)
    any(r -> !isempty(r.checked.escaped), study.results) && error("containment failed")
    groups = equivalent_coordinates(study)
    decisions = admission(study, groups)
    println("PROPOSAL DIMENSION: ", count(d -> d.admitted, decisions))
    open(joinpath(destination, "equivalent-coordinates.csv"), "w") do io
        println(io,
            "coordinate,family,quantity,circuit,member,nominal,lower,upper,sharing_error,source_drivers,admitted,maximum_change")
        for (i, group) in enumerate(groups), member in group.members

            println(io,
                join(
                    (group.name, group.family, group.quantity,
                        study.results[member.row].row.name,
                        member.key, member.nominal, member.lower, member.upper, group.error[1],
                        join(string.(member.drivers), ';'), decisions[i].admitted, decisions[i].change),
                    ','))
        end
    end
    open(joinpath(destination, "equivalent-admission.csv"), "w") do io
        println(io, "coordinate,admitted,maximum_change,circuit,sequence,quantity,frequency,endpoint")
        for decision in decisions
            println(io, join(values(decision), ','))
        end
    end
    open(joinpath(destination, "equivalent-source-map.csv"), "w") do io
        println(io, "source_coordinate,equivalent_count,equivalent_coordinates")
        for source in study.inputs.box.names
            names = [group.name
                     for (i, group) in enumerate(groups)
                     if decisions[i].admitted &&
                any(m -> source in m.drivers, group.members)]
            println(io, join((source, length(names), join(string.(names), ';')), ','))
        end
    end
    open(joinpath(destination, "equivalent-sharing.csv"), "w") do io
        println(io, "first,second,family,quantity,maximum_normalized_difference")
        for (i, a) in enumerate(groups), b in groups[(i + 1):end]

            a.family == b.family && a.quantity == b.quantity || continue
            difference = maximum(maximum(abs.(x.profile .- y.profile)) for x in a.members
            for y in b.members)
            println(io, join((a.name, b.name, a.family, a.quantity, difference), ','))
        end
    end
    pairs = decorrelations(study, groups, decisions)
    open(joinpath(destination, "equivalent-decorrelations.csv"), "w") do io
        println(io,
            "first,second,shared_sources,first_endpoint,second_endpoint,first_value,second_value,joint_residual,best_source_point")
        for row in pairs
            print(io, join(values(row)[1:8], ','))
            println(io, ",\"", replace(row.source, '"' => "\"\""), "\"")
        end
    end
    println("PROPOSAL DIMENSION: ", count(d -> d.admitted, decisions))
    return (; study, groups, decisions, pairs)
end

"""
    reduction_table(path)

Read the reduction's CSV tables. A final quoted source-point field remains text.
The preceding fields use commas as delimiters and occupy one line.
"""
function reduction_table(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    return [Dict(zip(header, split(line, ','; limit = length(header))))
            for line in lines[2:end]]
end

plain(value) = value
plain(value::Symbol) = string(value)
plain(value::Tuple) = [plain(v) for v in value]
plain(value::AbstractVector) = [plain(v) for v in value]
function plain(value::NamedTuple)
    Dict{String, Any}(string(k) => plain(v) for (k, v) in pairs(value))
end

function toml_text(value)
    io = IOBuffer()
    TOML.print(io, value; sorted = true)
    return String(take!(io))
end

"""
    generated_equivalent(inputs, evidence)

Construct deterministic case and data documents from the computed reduction
and nominal public blueprints. The output preserves the source branch assignments.
"""
function generated_equivalent(inputs, evidence)
    members = reduction_table(joinpath(evidence, "equivalent-coordinates.csv"))
    decisions = reduction_table(joinpath(evidence, "equivalent-admission.csv"))
    admitted = [r["coordinate"] for r in decisions if parse(Bool, r["admitted"])]
    source = relpath(inputs.file, STUDIES)
    commit = readchomp(`git -C $STUDIES rev-parse HEAD`)
    command = "julia --project=examples/rhhc-variants examples/rhhc-variants/reduce.jl " *
              source
    method = "effective coordinates, subdivision, constrained hull search, exact sharing at $(METHOD.sharing_tolerance), $(100METHOD.admission_threshold) percent admission"
    origin = "$source at $commit, LCM public flatten and the computed equivalent hulls"
    data_path = "data/equivalent/" * splitext(basename(inputs.file))[1] * ".toml"
    case_path = "cases/" * splitext(basename(inputs.file))[1] * "-equivalent.toml"
    coordinates_data = Dict{String, Any}()
    for name in admitted
        rows = filter(r -> r["coordinate"] == name, members)
        meanings = Dict(
            "width" => ("Radial width of an equivalent tube or annular layer", "m"),
            "resistance" =>
                ("Operating DC resistance per unit length of an equivalent tube", "ohm/m"),
            "mu_r" => ("Relative permeability of an equivalent tube", "1"),
            "eps_r" => ("Relative permittivity of an annular layer", "1"),
            "y" => ("Vertical placement of a circuit member", "m"),
            "rho" => ("Homogeneous earth resistivity", "ohm*m"))
        definition, unit = meanings[first(rows)["quantity"]]
        coordinates_data[name] = Dict("family" => first(rows)["family"],
            "definition" => definition, "unit" => "1", "physical_unit" => unit,
            "nominal" => 0.0, "lower" => -1.0,
            "upper" => 1.0, "scale" => "linear", "provenance" => origin,
            "sharing_error" => parse(Float64, first(rows)["sharing_error"]),
            "members" => [Dict("circuit" => r["circuit"], "field" => r["member"],
                 "nominal" => parse(Float64, r["nominal"]),
                 "lower" => parse(Float64, r["lower"]), "upper" =>
                     parse(Float64, r["upper"]),
                 "scale" => "linear", "unit" => unit, "provenance" => origin)
             for r in rows])
    end
    circuits_data = Dict{String, Any}()
    p = nominal_point(inputs.box)
    for row in inputs.definitions
        original = row.problem(p, 50.0)
        values = Dict(equivalent_values(original, row.layout.kind))
        entry = plain(specification(first(original.system.designs), 1, values))
        entry["family_name"] = string(row.spec.family_name)
        entry["provenance"] = "$origin, nominal circuit $(row.name)"
        entry["earth"] = Dict{String, Any}("rho" => Dict("nominal" => values["earth.rho"]),
            "eps" => Dict("nominal" => values["earth.eps_r"]))
        if haskey(values, "formation.y")
            entry["center_y"] = Dict("nominal" => values["formation.y"])
        else
            entry["phase_x"] = [v.x for v in original.system.placements[1:3]]
            entry["phase_y"] = Dict{String, Any}[Dict("nominal" => values["cable.$k.y"])
                                                 for k in 1:3]
            if haskey(row.layout, :earth_wire)
                entry["shield"] = plain(specification(last(original.system.designs), 4, values))
                entry["shield_x"] = last(original.system.placements).x
                entry["shield_y"] = Dict("nominal" => values["cable.4.y"])
            end
        end
        for name in admitted, member in coordinates_data[name]["members"]

            member["circuit"] == row.name || continue
            physical = Dict("nominal" => member["nominal"], "lower" => member["lower"],
                "upper" => member["upper"], "coordinate" => name, "provenance" => origin)
            chunks = split(member["field"], '.')
            if first(chunks) == "earth"
                entry["earth"][chunks[2] == "eps_r" ? "eps" : chunks[2]] = physical
            elseif first(chunks) == "cable"
                cable = parse(Int, chunks[2])
                if length(chunks) == 3
                    if cable <= 3
                        entry["phase_y"][cable] = physical
                    else
                        entry["shield_y"] = physical
                    end
                elseif cable == 1 || cable == 4
                    owner = cable == 1 ? entry : entry["shield"]
                    key = join(chunks[3:4], '.')
                    index = findfirst(part -> part["key"] == key, owner["parts"])
                    owner["parts"][index][last(chunks)] = physical
                end
            end
        end
        circuits_data[row.name] = entry
    end
    lcm_commit = readchomp(`git -C $(pkgdir(LineCableModels)) rev-parse HEAD`)
    lcm_changes = read(`git -C $(pkgdir(LineCableModels)) diff HEAD -- src Project.toml`, String)
    lcm_fingerprint = bytes2hex(sha256(lcm_changes))
    data = Dict(
        "source_case" => source, "source_commit" => commit, "lcm_commit" => lcm_commit,
        "lcm_source_diff_sha256" => lcm_fingerprint,
        "method" => method, "regenerate" => command,
        "coordinates" => coordinates_data, "circuits" => circuits_data)
    config = TOML.parsefile(inputs.file)
    config["box_names"] = admitted
    config["equivalent_data"] = data_path
    config["box"] = Dict(name => Dict("family" => coordinates_data[name]["family"],
                             "coordinate" => name, "origin" => "$data_path, coordinates.$name")
    for name in admitted)
    config["bindings"] = Dict{String, Any}()
    for name in admitted
        family = coordinates_data[name]["family"]
        bindings = get!(config["bindings"], family, Dict{String, Any}())
        bindings[name] = name
    end
    config["note"] = "Generated equivalent intervals with solid grounding and ideal transposition. " *
                     "Source correlations are dropped between independent coordinates. " *
                     "Five screened radial widths keep their nominal values."
    header = "# Generated. Source: $source\n# Source commit: $commit\n# LCM commit: $lcm_commit\n# Method: $method\n# Regenerate: `$command`\n\n"
    return (; case_path, data_path, case_text = header * toml_text(config),
        data_text = header * toml_text(data), data, admitted, decisions, source, commit, command)
end

"""
    reduction_report(generated, evidence)

Write the compact reduction evidence from the computed numerical tables.
"""
function reduction_report(generated, evidence)
    io = IOBuffer()
    digest = bytes2hex(sha256(generated.data_text))
    println(io, "Equivalent uncertainty reduction")
    println(io, "Source: ", generated.source, "; source commit: ", generated.commit)
    println(io, "LCM commit: ", generated.data["lcm_commit"])
    println(io, "Method: ", generated.data["method"])
    println(io, "Regenerate: ", generated.command)
    println(io, "Data: ", generated.data_path, "; SHA-256: ", digest)
    source_dimension = length(TOML.parsefile(joinpath(STUDIES, generated.source))["box_names"])
    println(io, "Source dimension: ", source_dimension, ". Equivalent dimension: ",
        length(generated.admitted), ". Vertex threshold: 10.")
    println(io, "Normalized coordinates have nominal 0, bounds [-1, 1], and a linear scale.")
    println(io, "Physical values interpolate separately below and above nominal.")
    println(io, "Five screened widths stay nominal. Containment applies to the admitted physical members.")
    println(io,
        "coordinate | definition | members | maximum sharing error | maximum endpoint change | physical bounds")
    for name in generated.admitted
        coordinate = generated.data["coordinates"][name]
        decision = only(filter(r -> r["coordinate"] == name, generated.decisions))
        bounds = unique("$(m["circuit"]): $(m["field"]) [$(m["lower"]), $(m["upper"])]"
        for m in coordinate["members"]
        if !startswith(m["field"], "cable.2.") &&
            !startswith(m["field"], "cable.3."))
        println(
            io, name, " | ", coordinate["definition"], " | ", length(coordinate["members"]),
            " | ", coordinate["sharing_error"], " | ",
            decision["maximum_change"], " | ", join(bounds, "; "))
    end
    cells = reduction_table(joinpath(evidence, "equivalent-cells.csv"))
    println(io, "Subdivision cells (other finite quantities have one cell):")
    for cell in cells
        cell["quantity"] == "raw.cable.1.y" && cell["coordinate"] == "theta_ACSR" ||
            continue
        println(
            io, cell["circuit"], ": cell ", cell["cell"], ", effective phase temperature [",
            cell["lower"], ", ", cell["upper"], "] degrees C, depth ", cell["depth"],
            ", ", cell["samples"], " samples. Each phase uses this split.")
    end
    println(io, "Containment:")
    for row in reduction_table(joinpath(evidence, "equivalent-containment.csv"))
        println(
            io, row["circuit"], ": ", row["samples"], " mapped samples, ", row["escapes"],
            " escapes, largest normalized outward margin ", row["largest_outward_margin"])
    end
    pairs = reduction_table(joinpath(evidence, "equivalent-decorrelations.csv"))
    sort!(pairs; by = r -> (-parse(Float64, r["joint_residual"]), r["first"], r["second"]))
    println(io, "Decorrelated pairs: ", length(pairs),
        ". Five largest added extremes by joint-search residual:")
    for row in pairs[1:min(5, length(pairs))]
        println(io, row["first"], " = ", row["first_value"], " at ", row["first_endpoint"],
            "; ", row["second"], " = ", row["second_value"], " at ", row["second_endpoint"],
            "; shared sources ", row["shared_sources"], "; joint-search residual ", row["joint_residual"])
    end
    return String(take!(io))
end

"""
    generate_equivalent(inputs, evidence; destination = STUDIES)

Write the generated case, coordinate data and compact reduction report.
Raw numerical tables remain in the supplied temporary evidence directory.
"""
function generate_equivalent(inputs, evidence; destination = STUDIES)
    generated = generated_equivalent(inputs, evidence)
    case_name = splitext(basename(generated.case_path))[1]
    report_path = "output/" * case_name * "-reduction.txt"
    files = Dict(generated.case_path => generated.case_text,
        generated.data_path => generated.data_text,
        report_path => reduction_report(generated, evidence))
    for (path, contents) in files
        file = joinpath(destination, path)
        mkpath(dirname(file))
        write(file, contents)
    end
    return (; generated, files)
end

if abspath(PROGRAM_FILE) == @__FILE__
    inputs = reduction_inputs(isempty(ARGS) ? "cases/ieee30-full.toml" : first(ARGS))
    destination = length(ARGS) > 1 ? ARGS[2] : "/tmp/hpm-equivalent-step3"
    study_result = hull_study(inputs, destination)
    proposal(study_result, destination)
    generate_equivalent(inputs, destination)
end
