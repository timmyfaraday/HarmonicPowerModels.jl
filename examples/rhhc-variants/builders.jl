using HarmonicPowerModels, LineCableModels, PowerModels, JuMP, Ipopt, Clarabel
using Printf, LinearAlgebra, TOML

const HPM = HarmonicPowerModels
const STUDIES = @__DIR__
const METHOD = (difference_step = 1e-4, zero_band = 1e-6,
    subdivision_depth = 6, sharing_tolerance = 1e-10,
    verification_samples = 4096, admission_threshold = 1e-3,
    interior_starts = (corners = true, center = true),
    interior_tolerance = 1e-10, interior_iterations = 32)
include(joinpath(STUDIES, "..", "shhc", "pce_params.jl"))
include(joinpath(STUDIES, "..", "shhc", "ieee30_scenarios.jl"))
PowerModels.silence()
BLAS.set_num_threads(1)

# Decode each TOML record once. The files keep the quantities and their citations.
record(x) = x
record(x::AbstractVector) = Tuple(record(v) for v in x)
function record(x::AbstractDict)
    haskey(x, "value") && return record(x["value"])
    return (; (Symbol(k) => record(v) for (k, v) in x)...)
end
readrecord(path) = record(TOML.parsefile(path))

function metal(row, reference; rho = row.rho,
        mu_r = get(row, :mu_r, reference.mu_r), temperature = reference.temperature)
    return Material(:conductor,
        rho * (1 + row.alpha * (temperature - reference.temperature)),
        get(row, :eps_r, reference.eps_r), mu_r, reference.temperature, row.alpha)
end
function dielectric(row; eps_r = row.eps_r)
    Material(:insulator, row.rho, eps_r; tan_delta = get(row, :tan_delta, 0.0))
end
semiconductor(row) = Material(:semicon, row.rho, row.eps_r)

"""
    equivalent(source::CableDesign)

Build plain tubes and annular layers from the public LCM blueprint of `source`.
Keep terminal names and materials. Radii are in meters. This path serves designs
with one concentric assembly at the origin and no helical dielectric correction.
"""
function equivalent(source::CableDesign)
    blueprint = LineCableModels.flatten(source)
    length(blueprint.positions) == 1 && only(blueprint.positions) == (0, 0) ||
        error("equivalent design needs one assembly at the origin")
    parts = AbstractCablePart[]
    for (k, tube) in enumerate(blueprint.conductors)
        shape = iszero(tube.r_in) ? Circle(tube.r_ex) : Annulus(tube.r_in, tube.r_ex)
        push!(parts, Terminal(source.terminals[k].name,
            Region(Symbol("tube_", k), shape, tube.material)))
        for l in blueprint.dielectric_ranges[k]
            layer = blueprint.dielectrics[l]
            push!(parts, Region(Symbol("layer_", l),
                Annulus(layer.r_in, layer.r_ex), layer.material))
        end
    end
    return build(CableDesign, source.id, parts...)
end

function specification(name)
    row = readrecord(joinpath(STUDIES, "data", "specs", name * ".toml"))
    return merge(row, (
        family = Val(Symbol(row.family)), family_name = Symbol(row.family), name = name))
end

function formation(name)
    row = readrecord(joinpath(STUDIES, "data", "formations", name * ".toml"))
    haskey(row, :structure) && (row = merge(row, (structure = Symbol(row.structure),)))
    return merge(row,
        (name = name, kind = Val(Symbol(row.kind)),
            reference = Val(Symbol(get(row, :reference, "support")))))
end

function coordinates(family, bindings, p)
    return (;
        (k => (haskey(bindings, k) ? get(p, getproperty(bindings, k), row.nominal) :
          row.nominal)
    for (k, row) in pairs(family))...)
end

# A signed normalized coordinate selects each product's registered interval.
function interpolate(row, coordinate)
    width = coordinate < zero(coordinate) ? row.nominal - row.lower :
            row.upper - row.nominal
    return row.nominal + coordinate * width
end

function design(::Val{:TOM}, spec, q, mat; temperature = q.theta)
    aluminum = metal(mat.tom_aluminium, mat.reference; temperature)
    copper = metal(mat.tom_copper, mat.reference; temperature)
    semicon = semiconductor(mat.tom_semicon)
    xlpe = dielectric(mat.tom_xlpe; eps_r = q.eps_xlpe)
    pe = dielectric(mat.tom_pe)
    return @design "TOM-630" begin
        @terminal :core begin
            core(aluminum; r = q.r_core)
        end
        screen(semicon; t = spec.t_sc1, tag = :semicon_inner)
        insulation(xlpe; t = q.t_xlpe)
        screen(semicon; t = spec.t_sc2, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = spec.t_tube)
        end
        jacket(pe; t = spec.t_jkt)
    end
end

function design(::Val{:U33}, spec, q, mat; temperature = q.theta)
    radius = interpolate(spec.d_c, q.d_c) / 2
    aluminum = metal(mat.aluminium, mat.reference;
        rho = interpolate(spec.R20, q.R20) * π * radius^2, temperature)
    copper = metal(mat.copper, mat.reference; rho = q.rho_copper, temperature)
    inner = semiconductor(merge(mat.semicon_inner, (
        rho = q.rho_inner, eps_r = q.eps_semicon)))
    outer = semiconductor(merge(mat.semicon_outer, (
        rho = q.rho_outer, eps_r = q.eps_semicon)))
    xlpe = dielectric(merge(mat.xlpe, (rho = q.rho_xlpe, tan_delta = q.tand_xlpe)); eps_r = q.eps)
    pvc = dielectric(merge(mat.pvc, (rho = q.rho_jacket, eps_r = q.eps_jacket)))
    inner_thickness = interpolate(spec.t_sc1, q.t_sc1)
    screen_radius = radius + inner_thickness + q.t_ins + q.t_sc2
    tube = sqrt(screen_radius^2 + q.rho_copper / (π * q.Rscreen)) - screen_radius
    thickness = interpolate(spec.t_jkt, q.t_jkt) + q.d_wire - tube
    return @design "U33" begin
        @terminal :core begin
            core(aluminum; r = radius)
        end
        screen(inner; t = inner_thickness, tag = :semicon_inner)
        insulation(xlpe; t = q.t_ins)
        screen(outer; t = q.t_sc2, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = tube)
        end
        jacket(pvc; t = thickness)
    end
end

# A concentric strand has 1 + 3k(k + 1) wires and spans 2k + 1 diameters.
strandspan(count) = 2 * ((isqrt(12count - 3) - 3) ÷ 6) + 1

"""
    construction(spec, q, inner, outer)

Resolve the steel course counts and the dimensionless steel and aluminum lay ratios.
Keep the supplied effective inner and outer aluminum lays unchanged.
"""
function construction(spec, q, inner, outer)
    steel_courses = Tuple(6i for i in 1:((strandspan(spec.n_steel) - 1) ÷ 2))
    steel_lays = accumulate(min,
        Tuple(i == 1 ? q.lay_steel : q.lay_steel_outer for i in eachindex(steel_courses)))
    aluminum_lays = Tuple(i == length(spec.courses) ? outer :
                          i == 1 ? inner : min(q.lay_middle, inner)
    for i in eachindex(spec.courses))
    return (; steel_courses, steel_lays, aluminum_lays)
end

"""
    construction(spec, q, reference, da, ds, inner, outer)

Compute mass in kg/m, metal area in square meters and outside
diameter in meters from the supplied wire diameters. Each component mass scales
with wire area and the ratio of its summed helical wire lengths to their nominal values.
The caller supplies the effective inner and outer aluminum lays unchanged.
"""
function construction(spec, q, reference, da, ds, inner, outer)
    recipe = construction(spec, q, inner, outer)
    nominal_lays = (i == length(spec.courses) ? reference.lay_outer :
                    i == 1 ? reference.lay_inner : reference.lay_middle
    for i in eachindex(spec.courses))
    aluminum = spec.mass_al *
               sum(n * sqrt(1 + (π / lay)^2)
    for (n, lay) in zip(spec.courses, recipe.aluminum_lays)) /
               sum(n * sqrt(1 + (π / lay)^2)
    for (n, lay) in zip(spec.courses, nominal_lays))
    steel = spec.mass_st *
            (1 + sum(
        (n * sqrt(1 + (π / lay)^2)
        for (n, lay) in zip(recipe.steel_courses, recipe.steel_lays)); init = 0.0)) /
            (1 + sum(
        (n * sqrt(1 + (π / (i == 1 ? reference.lay_steel : reference.lay_steel_outer))^2)
        for (i, n) in enumerate(recipe.steel_courses));
        init = 0.0))
    diameter = strandspan(spec.n_steel) * ds + 2length(spec.courses) * da
    mass = aluminum * (da / spec.d_aluminium.nominal)^2 +
           steel * (ds / spec.d_steel.nominal)^2
    area = π * (spec.n_aluminium * da^2 + spec.n_steel * ds^2) / 4
    return (;
        da, ds, inner = first(recipe.aluminum_lays), outer = last(recipe.aluminum_lays),
        diameter, mass, area)
end

# EN 50182:2001 5.4 and 5.7.3 constrain diameter and mass jointly.
function realize(spec, q, reference)
    inner = min(q.lay_inner, spec.n_steel > 7 ?
                             min(q.lay_steel, q.lay_steel_outer) : q.lay_steel)
    outer = min(q.lay_outer, inner,
        length(spec.courses) > 2 ? q.lay_middle : inner)
    da = interpolate(spec.d_aluminium, q.d_aluminium)
    ds = interpolate(spec.d_steel, q.d_steel)
    physical = construction(spec, q, reference, da, ds, inner, outer)
    factor = one(da + ds)
    if !(spec.diameter.lower <= physical.diameter <= spec.diameter.upper &&
         spec.mass.lower <= physical.mass <= spec.mass.upper)
        low = zero(factor)
        high = factor
        for _ in 1:100
            trial = (low + high) / 2
            a = spec.d_aluminium.nominal + trial * (da - spec.d_aluminium.nominal)
            d = spec.d_steel.nominal + trial * (ds - spec.d_steel.nominal)
            physical = construction(spec, q, reference, a, d, inner, outer)
            if spec.diameter.lower <= physical.diameter <= spec.diameter.upper &&
               spec.mass.lower <= physical.mass <= spec.mass.upper
                low = trial
            else
                high = trial
            end
        end
        factor = low
        da = spec.d_aluminium.nominal + factor * (da - spec.d_aluminium.nominal)
        ds = spec.d_steel.nominal + factor * (ds - spec.d_steel.nominal)
    end
    return merge(construction(spec, q, reference, da, ds, inner, outer), (; factor))
end

"""
    effective(kind, spec, q, reference, structure, row)

Map family inputs to the values consumed by a circuit design. For ACSR, keep the
realized construction and phase temperature in degrees Celsius as separate values.
"""
effective(kind, spec, q, reference, structure, row) = (; q, reference, structure)
function effective(::Val{:ACSR}, spec, q, reference, structure, row)
    temperature = haskey(row, :structure) ? max(q.theta, structure.ambient) : q.theta
    return (; q, reference, structure, temperature, realized = realize(spec, q, reference))
end

"""
    physical_inputs(spec, families, layout, bindings, p)

Resolve a source point and apply the construction and temperature maps once.
Return a record of effective values and earth inputs.
"""
function physical_inputs(spec, families, layout, bindings, p)
    family = getproperty(families, spec.family_name)
    q = coordinates(family, get(bindings, spec.family_name, (;)), p)
    reference = coordinates(family, (;), (;))
    structure = haskey(layout, :structure) ?
                coordinates(getproperty(families, layout.structure), get(bindings, layout.structure, (;)), p) :
                layout
    ground = coordinates(families.EARTH, get(bindings, :EARTH, (;)), p)
    return physical_inputs(spec.family, spec,
        effective(spec.family, spec, q, reference, structure, layout), ground)
end

"""
    effective_problem(spec, physical, layout, mat, f)

Build a problem from effective construction values at frequency `f` in hertz.
The design path consumes the realized construction and phase temperature.
"""
function effective_problem(spec, physical, layout, mat, f)
    conductor, positions = arrangement(spec.family, spec, physical, layout, mat)
    circuit = system(positions.kind, conductor, positions)
    soil = homogeneous(rho = physical.ground.rho, eps_r = physical.ground.eps)
    return LineParametersProblem(circuit, soil;
        temperature = mat.reference.temperature, frequencies = [f])
end

function design(::Val{:ACSR}, spec, q, mat; temperature = q.theta, realized)
    steel = metal(mat.s1a, mat.reference; rho = q.rho_steel, mu_r = q.mu_steel, temperature)
    aluminum = metal(mat.al1, mat.reference; rho = q.rho_aluminium, temperature)
    da, ds = realized.da, realized.ds
    recipe = construction(spec, q, realized.inner, realized.outer)
    core = stranded(steel; center = Circle(ds / 2), shape = Circle(ds / 2),
        n = recipe.steel_courses, lay = LayRatio.(Tuple(recipe.steel_lays)),
        boundary = Circle(strandspan(spec.n_steel) * ds / 2), tag = :steel)
    courses = Tuple(wires(aluminum; shape = Circle(da / 2), n = count,
                        lay = LayRatio(lay), tag = Symbol("aluminium_", i))
    for (i, (count, lay)) in enumerate(zip(spec.courses, recipe.aluminum_lays)))
    return build(CableDesign, "ACSR", Terminal(:core, Stack(core, courses...)))
end

function design(::Val{:O132}, spec, q, mat; temperature = q.ambient)
    steel = metal(mat.s1a, mat.reference; rho = q.rho_steel, mu_r = q.mu_steel, temperature)
    return @design spec.name begin
        @terminal :core begin
            stranded(steel; center = Circle(q.diameter / 2),
                shape = Circle(q.diameter / 2), n = (spec.count - 1,),
                lay = LayRatio(q.lay), boundary = Circle(strandspan(spec.count) *
                                                         q.diameter / 2))
        end
    end
end

centre(::Val{:centre}, row, spacing) = Placement(0.0, -row.depth)
centre(::Val{:lower}, row, spacing) = Placement(0.0, -row.depth + spacing / (2sqrt(3)))

function system(::Val{:trefoil}, conductor, row)
    spacing = 2r_ex(conductor.outer)
    positions = trefoil(conductor; spacing, at = centre(row.reference, row, spacing))
    return build(CableSystem, fill(conductor, 3), positions;
        connections = [(core = k, sheath = 0) for k in 1:3], id = "trefoil", line_length = 1.0)
end

function system(::Val{:overhead}, conductor, row)
    # The parabolic span's mean height is support height minus two thirds of the sag.
    positions = [Placement(x, h - (2 / 3) * row.phase_sag)
                 for (x, h) in zip(row.phase_x, row.phase_height)]
    if !haskey(row, :earth_wire)
        return build(CableSystem, fill(conductor, 3), positions;
            connections = [(core = k,) for k in 1:3], id = row.name, line_length = 1.0)
    end
    push!(positions, Placement(row.shield_x, row.shield_height - (2 / 3) * row.shield_sag))
    return build(CableSystem, [conductor, conductor, conductor, row.shield], positions;
        connections = [(core = 1,), (core = 2,), (core = 3,), (core = 0,)],
        id = row.name, line_length = 1.0)
end

function problem(spec, families, row, mat, bindings)
    layout = haskey(row, :earth_wire) ?
             merge(row, (shield_spec = specification(row.earth_wire),)) : row
    return (p, f) -> effective_problem(spec,
        physical_inputs(spec, families, layout, bindings, p), layout, mat, f)
end

# The register uses a level parabolic span and small elastic and thermal strains.
function sag(span, mass, area, modulus, expansion, reference, theta, theta0, gravity)
    load = mass * gravity
    tension0 = load * span^2 / (8reference)
    low = zero(theta + modulus + reference)
    high = zero(low) + span
    for _ in 1:100
        value = (low + high) / 2
        tension = load * span^2 / (8value)
        residual = 8(value^2 - reference^2) / (3span^2) -
                   (tension - tension0) / (modulus * area) - expansion * (theta - theta0)
        if residual > 0
            high = value
        else
            low = value
        end
    end
    return (low + high) / 2
end

function arrangement(::Val{:TOM}, spec, q, reference, structure, row, mat)
    return design(spec.family, spec, q, mat), row
end

function arrangement(::Val{:U33}, spec, q, reference, structure, row, mat)
    return design(spec.family, spec, q, mat), merge(row, (depth = q.depth,))
end

function arrangement(kind, spec, physical, row, mat)
    return arrangement(
        kind, spec, physical.q, physical.reference, physical.structure, row, mat)
end

function arrangement(::Val{:ACSR}, spec, physical, row, mat)
    (; q, structure, temperature, realized) = physical
    conductor = design(spec.family, spec, q, mat; temperature, realized)
    haskey(row, :structure) || return conductor, row
    phase_sag = sag(structure.span, realized.mass, realized.area,
        get(spec, :E, q.E) * 1e9, get(spec, :expansion, q.expansion),
        structure.phase_sag, temperature, structure.reference_theta, structure.gravity)
    wire = row.shield_spec
    shield = design(wire.family, wire, structure, mat)
    wire_area = π * structure.diameter^2 / 4
    wire_mass = structure.density * wire_area *
                (1 + (wire.count - 1) * sqrt(1 + (π / structure.lay)^2))
    shield_sag = sag(structure.span, wire_mass, wire.count * wire_area, structure.E * 1e9,
        structure.expansion, structure.shield_sag, structure.ambient,
        structure.reference_theta, structure.gravity)
    return conductor, merge(row, (; phase_sag, shield_sag, shield))
end

# The scan uses three earthwire copies as phases to check their DC resistance.
function arrangement(::Val{:O132}, spec, q, reference, structure, row, mat)
    conductor = design(spec.family, spec, q, mat)
    return conductor,
    merge(row, (phase_sag = structure.phase_sag,
        shield_sag = structure.shield_sag, shield = conductor))
end

function parameterbox(rows, families, names)
    entries = map(names) do text
        name = Symbol(text)
        row = getproperty(rows, name)
        values = haskey(row, :family) ?
                 getproperty(getproperty(families, Symbol(row.family)),
            Symbol(row.coordinate)) : row
        return name =>
            (values.nominal, values.lower, values.upper,
                Symbol(get(values, :scale, "linear") == "registered" ? "linear" :
                       get(values, :scale, "linear")))
    end
    return ParameterBox(entries...)
end

function isxfmr(b)
    get(b, "transformer", false) || !(get(b, "tap", 1.0) in (0.0, 1.0)) ||
        get(b, "shift", 0.0) != 0.0
end

function transformers(rows, data)
    for (name, row) in pairs(rows)
        id = string(name)
        haskey(data["branch"], id) || error("declared transformer branch $id is absent")
        br = data["branch"][id]
        isxfmr(br) || error("declared branch $id is not a transformer")
        (br["f_bus"], br["t_bus"]) == (row.f_bus, row.t_bus) ||
            error("transformer $id does not join declared buses $(row.f_bus)–$(row.t_bus)")
    end
    return rows
end

# MATPOWER route comments remain the source of lengths, in meters after decoding.
function network(config, ipopt)
    source = config.network
    source.kind == "scenario" && return ieee30_scenario(source.scenario, ipopt)
    file = normpath(joinpath(STUDIES, source.file))
    data = HPM.parse_file(file)
    if source.kind == "matpower"
        text = read(file, String)
        table = match(r"mpc\.branch\s*=\s*\[(.*?)\];"s, text)
        for (id, line) in enumerate(filter(
            l -> !isempty(strip(l)) && !startswith(strip(l), "%"),
            split(table.captures[1], '\n')))
            route = match(r"(?:[lL]\s*=\s*|[,;]\s*)([0-9.]+)\s*km", line)
            if route !== nothing
                data["branch"][string(id)]["length_m"] = 1e3parse(Float64, route.captures[1])
            end
        end
    end
    return data
end

function route_length(br, kv, zbase, config)
    if haskey(config, :routes)
        row = getproperty(config.branches, Symbol(br["index"]))
        return getproperty(config.routes, Symbol(row.route)).length_m
    end
    haskey(br, "length_m") && return br["length_m"]
    haskey(br, "cable_km") && return 1e3br["cable_km"]
    config.network.kind == "scenario" && return 1e3br["br_x"] * zbase / OHL_X[kv]
    error("LCM branch $(br["index"]) has no route length in its source data")
end

function harmonic_units(count::Integer, data)
    sort(parse.(Int, collect(keys(data["load"]))))[1:count]
end
harmonic_units(ids::Tuple, data) = collect(ids)

function study(path; all_coordinates = false)
    file = isfile(path) ? abspath(path) : joinpath(STUDIES, path)
    config = readrecord(file)
    haskey(config, :factor_data) && return study(file, config,
        study(joinpath(STUDIES, config.source_case)))
    mat = readrecord(joinpath(STUDIES, "data", "materials.toml"))
    records = readrecord(config)
    families = records.families
    box = parameterbox(config.box, families,
        all_coordinates ? Tuple(string.(keys(config.box))) : config.box_names)
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    conic = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    data = network(config, ipopt)
    haskey(config, :principle) && (data["principle"] = config.principle)
    transformers(get(config, :transformers, (;)), data)
    bindings = (;
        (k => (; (n => Symbol(v) for (n, v) in pairs(b))...)
    for (k, b) in pairs(get(config, :bindings, (;))))...)
    definitions = Dict{Tuple{String, String}, Function}()
    details = Dict{Int, NamedTuple}()
    branches = Dict{Int, NamedTuple}()
    for (name, row) in pairs(get(config, :branches, (;)))
        b = parse(Int, string(name))
        br = data["branch"][string(b)]
        kv = data["bus"][string(br["f_bus"])]["base_kv"]
        zbase = kv^2 / data["baseMVA"]
        len = route_length(br, kv, zbase, config)
        if haskey(row, :rating)
            rating = getproperty(config.ratings, Symbol(row.rating))
            br["rate_a"] = sqrt(3) * kv * rating.current_A / (1000data["baseMVA"])
        end
        key = (row.spec, row.formation)
        if !haskey(definitions, key)
            spec = specification(row.spec, row.formation, records)
            definitions[key] = problem(
                spec, families, formation(row.formation), mat, bindings)
        end
        branches[b] = (problem = definitions[key], len = len, zbase = zbase)
        details[b] = (branch = b, f_bus = br["f_bus"], t_bus = br["t_bus"], kv = kv,
            spec = row.spec, formation = row.formation, len = len,
            circuit = row.spec * "/" * row.formation,
            source_x = br["br_x"] * zbase / (len / 1e3),
            source_c = (br["b_fr"] + br["b_to"]) / (2π * config.f1 * zbase * (len / 1e3)))
    end
    if config.model == "exact"
        formulation = CoaxialFormulation(transposition = true)
        model = ExactPiModel(branches, formulation, box; f1 = config.f1)
        p = nominal_point(box)
        for b in keys(branches), (k, v) in pairs(branch_data(model, b, p, 1.0, :positive))

            data["branch"][string(b)][string(k)] = v
        end
    end
    units = harmonic_units(config.units, data)
    all(id -> haskey(data["load"], string(id)), units) ||
        error("harmonic unit is absent from the source load table")
    hd = HPM.replicate(deepcopy(data); H = [1; collect(config.harmonics)])
    if config.model == "scaled"
        classes = Dict(b => Symbol(name) for (name, ids) in pairs(config.classes)
        for b in ids)
        model = ScaledBranchModel(hd, classes)
    end
    hd["angle_reference"] = "local"
    for h in config.harmonics, k in collect(keys(hd["nw"][string(h)]["load"]))

        parse(Int, k) in units || delete!(hd["nw"][string(h)]["load"], k)
    end
    pce = Dict{String, Any}(deepcopy(PCE_PARAMS)...)
    pce["harmonics"] = Dict{Int, Any}(h => (PCE_PARAMS["harmonics"][h].v_min,
                                          PCE_PARAMS["harmonics"][h].v_max)
    for h in config.harmonics)
    pce["harmonic_units"] = units
    pce["bound"] = config.bound
    pce["epsilon"] = config.epsilon
    ids = unique(b -> details[b].circuit, sort(collect(keys(details))))
    lcm = config.model == "exact" ? Tuple(model.branches[b].zy for b in ids) : ()
    return (; file, config, mat, families, box, model, details, hd, pce, ipopt, conic, lcm)
end

function circuits(s)
    ids = sort(collect(keys(s.details)))
    return unique(b -> s.details[b].circuit, ids)
end

function header(io, s)
    println(io, "git commit: ", readchomp(`git -C $STUDIES rev-parse HEAD`))
    println(io, "LCM commit: ", readchomp(`git -C $(pkgdir(LineCableModels)) rev-parse HEAD`))
    println(io, "case: ", relpath(s.file, dirname(dirname(STUDIES))))
    println(io, "threads: ", Threads.nthreads(), "; solver: Clarabel")
end

function declarations(io, s)
    isempty(get(s.config, :transformers, (;))) && return
    println(io, "Transformer declarations. Not applied to the network: HPM reads vector groups")
    println(io, "from xfmr, which has no tap. IEEE30 transformers are tapped branches.")
    for (id, row) in sort(collect(pairs(s.config.transformers)); by = x -> parse(Int, string(first(x))))
        br = s.hd["nw"]["1"]["branch"][string(id)]
        kv1 = s.hd["nw"]["1"]["bus"][string(br["f_bus"])]["base_kv"]
        kv2 = s.hd["nw"]["1"]["bus"][string(br["t_bus"])]["base_kv"]
        zb1, zb2 = kv1^2 / s.hd["nw"]["1"]["baseMVA"], kv2^2 / s.hd["nw"]["1"]["baseMVA"]
        println(io, "branch ", id, " (", row.f_bus, "–", row.t_bus, "): vg=", row.vg,
            " gnd1=", row.gnd1, " gnd2=", row.gnd2,
            " re1=", row.re1_ohm / zb1, " xe1=", row.xe1_ohm / zb1,
            " re2=", row.re2_ohm / zb2, " xe2=", row.xe2_ohm / zb2, " [pu]")
        println(io, row.note)
        haskey(row, :resistance_note) && println(io, row.resistance_note)
    end
end

"""
    equivalent_values(source, kind)

Describe plain-tube equivalents by positive radial widths in meters, operating
DC resistances in ohms per meter, and material properties. The widths reconstruct
the blueprint radii. Contact formations retain their geometric spacing law.
"""
function equivalent_values(source::LineParametersProblem, kind)
    rows = Pair{String, Float64}[]
    formulation = CoaxialFormulation(transposition = true)
    for (k, design) in enumerate(source.system.designs)
        blueprint = LineCableModels.flatten(design)
        for (i, tube) in enumerate(blueprint.conductors)
            prefix = "cable.$k.tube.$i."
            i == 1 && push!(rows, prefix * "hole" => tube.r_in)
            push!(rows, prefix * "width" => tube.r_ex - tube.r_in)
            rho = formulation.resistivity(tube.material, source.temperature)
            push!(rows, prefix * "resistance" => rho / (π * (tube.r_ex^2 - tube.r_in^2)))
            push!(rows, prefix * "mu_r" => tube.material.mu_r)
        end
        for (i, layer) in enumerate(blueprint.dielectrics)
            prefix = "cable.$k.layer.$i."
            push!(rows, prefix * "width" => layer.r_ex - layer.r_in)
            for field in (:rho, :eps_r, :mu_r, :tan_delta)
                push!(rows, prefix * string(field) => getproperty(layer.material, field))
            end
        end
    end
    append!(rows, placement_values(kind, source))
    length(source.earth.layers) == 2 || error("equivalent values need homogeneous earth")
    for field in (:rho, :eps_r, :mu_r)
        push!(rows, "earth.$field" => getproperty(last(source.earth.layers), field))
    end
    return rows
end

function placement_values(::Val{:trefoil}, source)
    conductor = first(source.system.designs)
    relative = trefoil(conductor; spacing = 2support(conductor.outer))
    center = first(source.system.placements) * inv(first(relative))
    return ["formation.x" => center.x, "formation.y" => center.y]
end
function placement_values(::Val{:overhead}, source)
    return ["cable.$k.$field" => getproperty(p, field)
            for (k, p) in enumerate(source.system.placements) for field in (:x, :y)]
end

"""
    equivalent(source, values, kind)

Reconstruct plain tubes and annular layers from `equivalent_values`. Tube
resistivity equals DC resistance times metal area at the operating temperature.
Keep the source connections and formation, with its equivalent placement inputs.
"""
function equivalent(source::LineParametersProblem, values::AbstractDict, kind)
    mat = (reference = (temperature = source.temperature,),)
    designs = [design(Val(:equivalent), specification(conductor, k, values), (;), mat)
               for (k, conductor) in enumerate(source.system.designs)]
    positions = equivalent_positions(kind, source, designs, values)
    connections = [(;
                       (design.terminals[row.terminal].name => row.connection
                   for row in source.system.terminals if row.cable == k)...)
                   for (k, design) in enumerate(source.system.designs)]
    circuit = build(CableSystem, designs, positions; connections,
        id = source.system.id, line_length = source.system.line_length)
    soil = homogeneous(rho = values["earth.rho"], eps_r = values["earth.eps_r"],
        mu_r = values["earth.mu_r"])
    return LineParametersProblem(circuit, soil;
        temperature = source.temperature, frequencies = source.frequencies)
end

function equivalent_positions(::Val{:trefoil}, source, designs, values)
    conductor = first(designs)
    return trefoil(conductor; spacing = 2support(conductor.outer),
        at = Placement(values["formation.x"], values["formation.y"]))
end
function equivalent_positions(::Val{:overhead}, source, designs, values)
    return [Placement(values["cable.$k.x"], values["cable.$k.y"])
            for k in eachindex(designs)]
end

"""
    readrecord(config::NamedTuple)

Read the physical families and, for a generated equivalent case, its coordinate
records. Keep the physical families available to the unchanged catalogue checks.
"""
function readrecord(config::NamedTuple)
    families = (;
        (Symbol(splitext(name)[1]) =>
             readrecord(joinpath(STUDIES, "data", "families", name))
    for name in sort(readdir(joinpath(STUDIES, "data", "families")))
    if endswith(name, ".toml"))...)
    haskey(config, :equivalent_data) || return (; families)
    dataset = readrecord(joinpath(STUDIES, config.equivalent_data))
    for (name, coordinate) in pairs(dataset.coordinates)
        family = Symbol(coordinate.family)
        declaration = (nominal = coordinate.nominal, lower = coordinate.lower,
            upper = coordinate.upper, scale = coordinate.scale)
        rows = merge(getproperty(families, family), NamedTuple{(name,)}((declaration,)))
        families = merge(families, NamedTuple{(family,)}((rows,)))
    end
    return (; families, dataset)
end

function specification(name, placement, records)
    haskey(records, :dataset) || return specification(name)
    row = getproperty(records.dataset.circuits, Symbol(name * "/" * placement))
    return merge(row, (;
        family = Val(:equivalent), family_name = Symbol(row.family_name), name))
end

"""
    interpolate(row, coordinates::NamedTuple)

Resolve a generated physical member from its signed coordinate, or keep its
nominal value when the member has no admitted coordinate.
"""
function interpolate(row, coordinates::NamedTuple)
    return haskey(row, :coordinate) ?
           interpolate(row, getproperty(coordinates, Symbol(row.coordinate))) : row.nominal
end

physical_inputs(kind, spec, physical, ground) = merge(physical, (; ground))
function physical_inputs(::Val{:equivalent}, spec, physical, ground)
    return merge(physical,
        (ground = (rho = interpolate(spec.earth.rho, ground),
            eps = interpolate(spec.earth.eps, ground)),))
end

function specification(original::Material{T, K}) where {T, K}
    return (; kind = string(K), T0 = original.T0, alpha = original.alpha,
        rho_thermal = original.rho_thermal, theta_max = original.theta_max,
        sigma_solar = original.sigma_solar)
end

"""
    specification(source::CableDesign, cable, values)

Describe the public blueprint as tube and layer declarations. Physical members
are nominal records until the generator attaches their admitted coordinates.
"""
function specification(source::CableDesign, cable, values)
    blueprint = LineCableModels.flatten(source)
    parts = NamedTuple[]
    for (i, tube) in enumerate(blueprint.conductors)
        prefix = "cable.$cable.tube.$i."
        part = merge(specification(tube.material),
            (; key = "tube.$i", terminal = string(source.terminals[i].name),
                width = (nominal = values[prefix * "width"],),
                resistance = (nominal = values[prefix * "resistance"],),
                mu_r = (nominal = values[prefix * "mu_r"],),
                eps_r = (nominal = tube.material.eps_r,),
                tan_delta = (nominal = tube.material.tan_delta,)))
        i == 1 && (part = merge(part, (hole = (nominal = values[prefix * "hole"],),)))
        push!(parts, part)
        for l in blueprint.dielectric_ranges[i]
            layer = blueprint.dielectrics[l]
            prefix = "cable.$cable.layer.$l."
            push!(parts,
                merge(specification(layer.material),
                    (; key = "layer.$l", width = (nominal = values[prefix * "width"],),
                        rho = (nominal = values[prefix * "rho"],),
                        eps_r = (nominal = values[prefix * "eps_r"],),
                        mu_r = (nominal = values[prefix * "mu_r"],),
                        tan_delta = (nominal = values[prefix * "tan_delta"],))))
        end
    end
    return (; id = source.id, parts = Tuple(parts))
end

"""
    design(::Val{:equivalent}, spec, coordinates, materials)

Build plain tubes and annular layers through LCM's public design API. Resolve
widths in meters and tube DC resistance in ohms per meter through `interpolate`.
Convert operating resistivity to the declared material reference temperature.
"""
function design(::Val{:equivalent}, spec, q, mat)
    parts = AbstractCablePart[]
    radius = 0.0
    for part in spec.parts
        haskey(part, :hole) && (radius = interpolate(part.hole, q))
        outer = radius + interpolate(part.width, q)
        rho = haskey(part, :terminal) ?
              interpolate(part.resistance, q) * π * (outer^2 - radius^2) /
              (1 + part.alpha * (mat.reference.temperature - part.T0)) :
              interpolate(part.rho, q)
        medium = Material(Symbol(part.kind), rho, interpolate(part.eps_r, q),
            interpolate(part.mu_r, q), part.T0, part.alpha;
            rho_thermal = part.rho_thermal, theta_max = part.theta_max,
            tan_delta = interpolate(part.tan_delta, q), sigma_solar = part.sigma_solar)
        shape = iszero(radius) ? Circle(outer) : Annulus(radius, outer)
        region = Region(Symbol(replace(part.key, '.' => '_')), shape, medium)
        push!(parts, haskey(part, :terminal) ? Terminal(Symbol(part.terminal), region) :
                     region)
        radius = outer
    end
    return build(CableDesign, spec.id, parts...)
end

function arrangement(::Val{:equivalent}, spec, physical, row, mat)
    return arrangement(row.kind, spec, physical.q, physical.structure, row, mat)
end
function arrangement(::Val{:trefoil}, spec, q, structure, row, mat)
    conductor = design(Val(:equivalent), spec, q, mat)
    spacing = 2r_ex(conductor.outer)
    offset = centre(row.reference, merge(row, (depth = 0.0,)), spacing).y
    return conductor, merge(row, (depth = offset - interpolate(spec.center_y, q),))
end
function arrangement(::Val{:overhead}, spec, q, structure, row, mat)
    conductor = design(Val(:equivalent), spec, q, mat)
    positions = merge(row,
        (phase_sag = 0.0, phase_x = spec.phase_x,
            phase_height = Tuple(interpolate(member, q) for member in spec.phase_y)))
    haskey(spec, :shield) || return conductor, positions
    shield = design(Val(:equivalent), spec.shield, structure, mat)
    return conductor,
    merge(positions,
        (shield_sag = 0.0,
            shield_x = spec.shield_x, shield_height = interpolate(spec.shield_y, structure), shield))
end

function factor_values(values, p, family)
    kR, kL, kC = (get(p, Symbol("k", name, "_", family), 1.0) for name in ("R", "L", "C"))
    return (z = complex(kR * real(values.z), kL * imag(values.z)),
        y = complex(real(values.y), kC * imag(values.y)))
end

function factor_callable(original, nominal, family)
    return (p, f, sequences...) -> begin
        values = original(nominal, f)
        if isempty(sequences)
            return (;
                (name => factor_values(value, p, family)
            for (name, value) in pairs(values))...)
        end
        return factor_values(getproperty(values, only(sequences)), p, family)
    end
end

function study(file, config, source)
    dataset = readrecord(joinpath(STUDIES, config.factor_data))
    declarations = (;
        (name => (nominal = row.nominal, lower = row.lower,
             upper = row.upper, scale = row.scale)
    for (name, row) in pairs(dataset.coordinates))...)
    box = parameterbox(declarations, source.families, config.box_names)
    return study(file, config, source, box)
end

function study(file, config, source, box)
    classes = Dict(b => Symbol(name) for (name, ids) in pairs(config.classes) for b in ids)
    hd = deepcopy(source.hd)
    if config.model == "exact"
        p = nominal_point(source.box)
        branches = Dict(b => (zy = factor_callable(row.zy, p, classes[b]),
                            len = row.len, zbase = row.zbase)
        for (b, row) in source.model.branches)
        model = ExactPiModel(branches; f1 = config.f1)
    else
        model = ScaledBranchModel(hd, classes)
    end
    return merge(source, (; file, config, box, model, hd))
end

lcm_solves(s) = sum((length(c.sequences) for c in s.lcm); init = 0)
