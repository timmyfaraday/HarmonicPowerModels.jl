using HarmonicPowerModels, LineCableModels, PowerModels, JuMP, Ipopt, Clarabel
using Printf, LinearAlgebra, TOML

const HPM = HarmonicPowerModels
const STUDIES = @__DIR__
const METHOD = (admission_threshold = 1e-3,)
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
    construction(spec, q, reference, da, ds, inner, outer)

Compute mass, metal area and outside diameter from effective wire diameters in
meters and effective lay ratios. The caller supplies already ordered and projected values.
"""
function construction(spec, q, reference, da, ds, inner, outer)
    aluminum = spec.mass_al *
               (spec.courses[1] * sqrt(1 + (π / inner)^2) +
                spec.courses[2] * sqrt(1 + (π / outer)^2)) /
               (spec.courses[1] * sqrt(1 + (π / reference.lay_inner)^2) +
                spec.courses[2] * sqrt(1 + (π / reference.lay_outer)^2))
    steel = spec.mass_st * (1 + (spec.n_steel - 1) * sqrt(1 + (π / q.lay_steel)^2)) /
            (1 + (spec.n_steel - 1) * sqrt(1 + (π / reference.lay_steel)^2))
    diameter = strandspan(spec.n_steel) * ds + 2length(spec.courses) * da
    mass = aluminum * (da / spec.d_aluminium.nominal)^2 +
           steel * (ds / spec.d_steel.nominal)^2
    area = π * (spec.n_aluminium * da^2 + spec.n_steel * ds^2) / 4
    return (; da, ds, inner, outer, diameter, mass, area)
end

# EN 50182:2001 5.4 and 5.7.3 constrain diameter and mass jointly.
function realize(spec, q, reference)
    inner = q.lay_inner
    outer = min(q.lay_outer, inner)
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
    return @design "ACSR" begin
        @terminal :core begin
            stranded(steel; center = Circle(ds / 2),
                shape = Circle(ds / 2), n = (spec.n_steel - 1,),
                lay = LayRatio(q.lay_steel), boundary = Circle(strandspan(spec.n_steel) *
                                                               ds / 2),
                tag = :steel)
            wires(aluminum; shape = Circle(da / 2), n = spec.courses[1],
                lay = LayRatio(realized.inner), tag = :aluminium_inner)
            wires(aluminum; shape = Circle(da / 2), n = spec.courses[2],
                lay = LayRatio(realized.outer), tag = :aluminium_outer)
        end
    end
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
    phase_sag = sag(structure.span, realized.mass, realized.area, q.E * 1e9, q.expansion,
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
        br["vg"] = row.vg
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
        key = (row.spec, row.formation)
        if !haskey(definitions, key)
            spec = specification(row.spec)
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
    println(io,
        "LCM commit: ",
        TOML.parsefile(joinpath(STUDIES, "Project.toml"))["sources"]["LineCableModels"]["rev"])
    println(io, "case: ", relpath(s.file, dirname(dirname(STUDIES))))
    println(io, "threads: ", Threads.nthreads(), "; solver: Clarabel")
end

function declarations(io, s)
    isempty(get(s.config, :transformers, (;))) && return
    println(io, "Transformer vector groups are applied to tapped branches.")
    println(io, "Winding earthing is declared below; HPM applies it only in its IV formulation.")
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

function readrecord(config::NamedTuple)
    families = (;
        (Symbol(splitext(name)[1]) =>
             readrecord(joinpath(STUDIES, "data", "families", name))
    for name in sort(readdir(joinpath(STUDIES, "data", "families")))
    if endswith(name, ".toml"))...)
    return (; families)
end

physical_inputs(kind, spec, physical, ground) = merge(physical, (; ground))

lcm_solves(s) = sum((length(c.sequences) for c in s.lcm); init = 0)
