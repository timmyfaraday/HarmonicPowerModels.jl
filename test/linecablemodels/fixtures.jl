# Physical inputs: LCM `examples/cable_30kv_630mm2_trefoil.jl`.
# Parameter bounds: HPM `examples/rhhc/ieee30_lcm.jl`.
function fixture(p, f)
    aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.0, 20.0, 0.00403)
    copper = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    semicon = Material(:semicon, 1000.0, 1000.0)
    xlpe = Material(:insulator, 1.0e14, p.eps_xlpe)
    polyethylene = Material(:insulator, 1.0e14, 2.3)
    design = @design "xlpe_18_30kv_1x630al" begin
        @terminal :core begin
            core(aluminum; r = p.r_core, tag = :core)
        end
        screen(semicon; t = 0.8e-3, tag = :semicon_inner)
        insulation(xlpe; t = p.t_xlpe, tag = :insulation)
        screen(semicon; t = 0.8e-3, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = 0.3e-3, tag = :screen)
        end
        jacket(polyethylene; t = 2.5e-3, tag = :jacket)
    end
    spacing = 2r_ex(design.outer) # [m]
    placements = trefoil(design;
        spacing, at = Placement(0.0, -1.0 + spacing / (2sqrt(3))))
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = k, sheath = 0) for k in 1:3],
        id = "cable_30kv_630mm2_trefoil", line_length = 1.0)
    earth = homogeneous(rho = p.rho_earth, eps_r = 10.0) # [Ω·m]
    return LineParametersProblem(system, earth; temperature = 20.0, frequencies = [f])
end

# Flat formation at twice the cable diameter, specified for the warning test.
function flat_fixture(p, f)
    aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.0, 20.0, 0.00403)
    copper = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    semicon = Material(:semicon, 1000.0, 1000.0)
    xlpe = Material(:insulator, 1.0e14, p.eps_xlpe)
    polyethylene = Material(:insulator, 1.0e14, 2.3)
    design = @design "xlpe_18_30kv_1x630al" begin
        @terminal :core begin
            core(aluminum; r = p.r_core, tag = :core)
        end
        screen(semicon; t = 0.8e-3, tag = :semicon_inner)
        insulation(xlpe; t = p.t_xlpe, tag = :insulation)
        screen(semicon; t = 0.8e-3, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = 0.3e-3, tag = :screen)
        end
        jacket(polyethylene; t = 2.5e-3, tag = :jacket)
    end
    spacing = 2 * 2r_ex(design.outer) # [m]
    placements = flat(design; n = 3, spacing, at = Placement(0.0, -1.0))
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = k, sheath = 0) for k in 1:3],
        id = "cable_30kv_630mm2_flat", line_length = 1.0)
    earth = homogeneous(rho = p.rho_earth, eps_r = 10.0) # [Ω·m]
    return LineParametersProblem(system, earth; temperature = 20.0, frequencies = [f])
end

# A distinct function with the same physical input, to test shared branches.
second_fixture(p, f) = fixture(p, f)

const fixture_calls = Threads.Atomic{Int}(0)
function counted_fixture(p, f)
    Threads.atomic_add!(fixture_calls, 1)
    return fixture(p, f)
end

const box = ParameterBox(
    :eps_xlpe => (2.3, 2.25, 2.5),
    :t_xlpe => (8.0e-3, 0.95 * 8.0e-3, 1.10 * 8.0e-3),
    :r_core => (14.16e-3, 0.99 * 14.16e-3, 1.01 * 14.16e-3),
    :rho_earth => (100.0, 10.0, 1000.0, :log))
const formulation = CoaxialFormulation(transposition = true)

# Independent sequence values from LCM's ideally transposed phase matrices.
function direct(p, f)
    result = compute(fixture(p, f), formulation)
    zs = result.Z[1, 1, 1]
    zm = result.Z[1, 2, 1]
    ys = result.Y[1, 1, 1]
    ym = result.Y[1, 2, 1]
    return (zero = (z = zs + 2zm, y = ys + 2ym),
        positive = (z = zs - zm, y = ys - ym),
        negative = (z = zs - zm, y = ys - ym))
end
