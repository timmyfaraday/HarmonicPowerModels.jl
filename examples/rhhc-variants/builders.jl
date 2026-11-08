using HarmonicPowerModels, LineCableModels, PowerModels, JuMP, Ipopt, Clarabel
using Printf, LinearAlgebra, TOML
import LineCableModels.Commons: nominal

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

# Candidate study MV cables use a solid core and an equal-area copper screen tube.
# The jacket fills the difference between the wire envelope and the tube.
function design(::Union{Val{:U36}, Val{:U11}}, spec, q, mat; temperature = q.theta)
    radius = q.d_c / 2
    aluminum = metal(mat.aluminium, mat.reference;
        rho = q.R20 * π * radius^2, temperature)
    copper = metal(mat.copper, mat.reference; rho = q.rho_copper, temperature)
    inner = Material(:semicon, q.rho_inner, q.eps_semicon)
    outer = Material(:semicon, q.rho_outer, q.eps_semicon)
    xlpe = Material(:insulator, q.rho_ins, q.eps; tan_delta = q.tand)
    protection = Material(:insulator, q.rho_jacket, q.eps_jacket)
    screen_radius = radius + q.t_sc1 + q.t_ins + q.t_sc2
    area = spec.screen_count * π * q.d_wire^2 / 4
    tube = sqrt(screen_radius^2 + area / π) - screen_radius
    return @design spec.name begin
        @terminal :core begin
            core(aluminum; r = radius)
        end
        screen(inner; t = q.t_sc1, tag = :semicon_inner)
        insulation(xlpe; t = q.t_ins)
        screen(outer; t = q.t_sc2, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = tube)
        end
        jacket(protection; t = q.t_jkt + q.d_wire - tube)
    end
end

# The coaxial formulation needs one concentric assembly per cable.
# Variant: three unscreened single cores replace the NYY multicore construction.
function design(::Val{:NYY}, spec, q, mat; temperature = q.theta)
    radius = q.d_c / 2
    copper = metal(mat.copper, mat.reference;
        rho = q.R20 * π * radius^2, temperature)
    pvc = Material(:insulator, q.rho_ins, q.eps; tan_delta = q.tand)
    protection = Material(:insulator, q.rho_jacket, q.eps_jacket)
    return @design spec.name begin
        @terminal :core begin
            core(copper; r = radius)
        end
        insulation(pvc; t = q.t_ins)
        jacket(protection; t = q.t_jkt)
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
    # The six-wire aluminum course needs a center circle at least one wire diameter across.
    core_diameter = first(spec.courses) == 6 ? max(strandspan(spec.n_steel) * ds, da) :
                    strandspan(spec.n_steel) * ds
    diameter = core_diameter + 2length(spec.courses) * da
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
        # The converged feasible point lies on one diameter or mass limit.
        a0 = spec.d_aluminium.nominal + factor * (da - spec.d_aluminium.nominal)
        d0 = spec.d_steel.nominal + factor * (ds - spec.d_steel.nominal)
        physical = construction(spec, q, reference, a0, d0, inner, outer)
        limits = ((:diameter, spec.diameter.lower), (:diameter, spec.diameter.upper),
            (:mass, spec.mass.lower), (:mass, spec.mass.upper))
        active = argmin(map(limits) do (quantity, bound)
            abs(nominal(getproperty(physical, quantity)) - bound) / abs(bound)
        end)
        quantity, bound = limits[active]
        delta_a = nominal(da) - spec.d_aluminium.nominal
        delta_d = nominal(ds) - spec.d_steel.nominal
        if quantity === :diameter
            slope = 2length(spec.courses) * delta_a + strandspan(spec.n_steel) * delta_d
        else
            # Mass is A*a^2 + S*d^2 at fixed lays. Evaluate its two
            # coefficients on plain inputs, then differentiate analytically.
            q0, ref0 = map(nominal, q), map(nominal, reference)
            i0, o0 = nominal(inner), nominal(outer)
            A = construction(spec, q0, ref0, 1.0, 0.0, i0, o0).mass
            S = construction(spec, q0, ref0, 0.0, 1.0, i0, o0).mass
            slope = 2A * nominal(a0) * delta_a + 2S * nominal(d0) * delta_d
        end
        correction = (getproperty(physical, quantity) - bound) / slope
        factor -= correction - nominal(correction)
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
    return merge(effective(spec.family, spec, q, reference, structure, layout), (; ground))
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
    core_diameter = first(spec.courses) == 6 ? max(strandspan(spec.n_steel) * ds, da) :
                    strandspan(spec.n_steel) * ds
    core = stranded(steel; center = Circle(ds / 2), shape = Circle(ds / 2),
        n = recipe.steel_courses, lay = LayRatio.(Tuple(recipe.steel_lays)),
        boundary = Circle(core_diameter / 2), tag = :steel)
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

function system(::Val{:unscreened}, conductor, row)
    spacing = 2r_ex(conductor.outer)
    positions = trefoil(conductor; spacing, at = centre(row.reference, row, spacing))
    return build(CableSystem, fill(conductor, 3), positions;
        connections = [(core = k,) for k in 1:3], id = "trefoil", line_length = 1.0)
end

function system(::Val{:overhead}, conductor, row)
    # The parabolic span's mean height is support height minus two thirds of the sag.
    positions = [Placement(x, h - (2 / 3) * row.phase_sag)
                 for (x, h) in zip(row.phase_x, row.phase_height)]
    if !haskey(row, :earth_wire)
        return build(CableSystem, fill(conductor, 3), positions;
            connections = [(core = k,) for k in 1:3], id = row.name, line_length = 1.0)
    end
    append!(positions,
        [Placement(x, h - (2 / 3) * row.shield_sag)
         for (x, h) in zip(row.shield_x, row.shield_height)])
    count = length(positions) - 3
    return build(CableSystem, vcat(fill(conductor, 3), fill(row.shield, count)), positions;
        connections = vcat([(core = k,) for k in 1:3], fill((core = 0,), count)),
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
    x0 = (low + high) / 2
    tension = load * span^2 / (8x0)
    residual = 8(x0^2 - reference^2) / (3span^2) -
               (tension - tension0) / (modulus * area) - expansion * (theta - theta0)
    residual_x = 16x0 / (3span^2) + load * span^2 / (8modulus * area * x0^2)
    correction = residual / residual_x
    return x0 - (correction - nominal(correction))
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

function arrangement(::Union{Val{:U36}, Val{:U11}}, spec, q, reference, structure, row, mat)
    return design(spec.family, spec, q, mat), merge(row, (depth = q.depth,))
end

function arrangement(::Val{:NYY}, spec, q, reference, structure, row, mat)
    return design(spec.family, spec, q, mat),
    merge(row, (depth = q.depth, kind = Val(:unscreened)))
end

function arrangement(::Val{:ACSR}, spec, physical, row, mat)
    (; q, structure, temperature, realized) = physical
    conductor = design(spec.family, spec, q, mat; temperature, realized)
    haskey(row, :structure) || return conductor, row
    phase_sag = sag(structure.span, realized.mass, realized.area,
        get(spec, :E, q.E) * 1e9, get(spec, :expansion, q.expansion),
        structure.phase_sag, temperature, structure.reference_theta, structure.gravity)
    haskey(row, :earth_wire) || return conductor, merge(row, (; phase_sag))
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
        br["gnd1"], br["gnd2"] = row.gnd1, row.gnd2
        zb1 = data["bus"][string(row.f_bus)]["base_kv"]^2 / data["baseMVA"]
        zb2 = data["bus"][string(row.t_bus)]["base_kv"]^2 / data["baseMVA"]
        br["re1"], br["xe1"] = row.re1_ohm / zb1, row.xe1_ohm / zb1
        br["re2"], br["xe2"] = row.re2_ohm / zb2, row.xe2_ohm / zb2
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
    mat = readrecord(joinpath(STUDIES, "data", "materials.toml"))
    families = (;
        (Symbol(splitext(name)[1]) =>
             readrecord(joinpath(STUDIES, "data", "families", name))
    for name in sort(readdir(joinpath(STUDIES, "data", "families")))
    if endswith(name, ".toml"))...)
    box = parameterbox(config.box, families,
        all_coordinates ? Tuple(string.(keys(config.box))) : config.box_names)
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    conic = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    data = network(config, ipopt)
    haskey(config, :principle) && (data["principle"] = config.principle)
    transformers(get(config, :transformers, (;)), data)
    for bus in get(config, :star_points, ())
        haskey(data["bus"], string(bus)) || error("declared star point $bus is absent")
        data["bus"][string(bus)]["star_point"] = true
    end
    for (name, row) in pairs(get(config, :generators, (;)))
        id = string(name)
        haskey(data["gen"], id) || error("declared generator $id is absent")
        generator = data["gen"][id]
        if get(row, :zero_sequence_open, false)
            generator["zero_sequence_open"] = true
        else
            generator["rsc0"], generator["xsc0"] = row.rsc0, row.xsc0
        end
    end
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
    recombination = get(config, :recombination, true)
    pce["recombination"] = recombination
    pce["harmonics"] = Dict{Int, Any}(h => (PCE_PARAMS["harmonics"][h].v_min,
                                          PCE_PARAMS["harmonics"][h].v_max)
    for h in config.harmonics if recombination || haskey(PCE_PARAMS["harmonics"], h))
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
    println(io, "Transformer vector groups are applied to tapped branches.")
    println(io, "HPM applies winding earthing in the zero sequence.")
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

lcm_solves(s) = sum((length(c.sequences) for c in s.lcm); init = 0)
