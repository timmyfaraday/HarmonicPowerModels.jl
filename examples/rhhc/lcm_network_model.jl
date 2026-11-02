################################################################################
# examples/rhhc/lcm_network_model.jl                                           #
# LineCableModels.jl as the parametric network model of the robust HHC.        #
#                                                                               #
# A 33 kV network of two XLPE cable sections and one overhead line            #
# (data/hello_world.m). LineCableModels computes the positive-sequence z1(f),  #
# y1(f) per meter of the cable and the line; `ExactPiModel` turns them into    #
# branch data at any harmonic order and any point of a parameter box on the   #
# cable inputs. The script scans the driving-point impedance at bus 3 over the #
# box vertices: the parallel resonance moves across the 7th harmonic.         #
#                                                                               #
# Requires LineCableModels.jl `release/v0.2.0` (not yet registered):          #
#   julia> ] develop <path to LineCableModels.jl checkout>                    #
################################################################################

using HarmonicPowerModels, PowerModels, LineCableModels, LinearAlgebra, Printf
const HPM = HarmonicPowerModels
PowerModels.silence()

# ── LineCableModels side ──────────────────────────────────────────────────────
# One 18/30 kV XLPE cable type (1x630 Al, trefoil at 1 m depth, sheaths grounded at both
# ends) and one 33 kV overhead line, each built from a named parameter tuple. The phases
# are ideally transposed, so the 3×3 matrices are circulant and z1 = zs − zm.

const CABLE = (rho_core = 2.8264e-8, r_core = 14.16e-3, t_xlpe = 8.0e-3, eps_xlpe = 2.3,
               rho_earth = 100.0)
const LINE  = (rho_wire = 2.8264e-8, r_wire = 10.0e-3, height = 10.0, spacing = 1.5,
               rho_earth = 100.0)

positive_sequence(r) = (r.Z[1, 1, 1] - r.Z[1, 2, 1], r.Y[1, 1, 1] - r.Y[1, 2, 1])

function cable_zy(p, f)
    core_metal = Material(:conductor, p.rho_core, 1.0, 1.0, 20.0, 0.00403)
    semicon    = Material(:semicon, 1_000.0, 1_000.0)
    xlpe       = Material(:insulator, 1.0e14, p.eps_xlpe)
    copper     = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    pe         = Material(:insulator, 1.0e14, 2.3)
    design = @cable "xlpe_18_30kv_1x630al" begin
        @terminal :core begin
            core(core_metal; r = p.r_core, tag = :core)
        end
        screen(semicon; t = 0.8e-3, tag = :semicon_inner)
        insulation(xlpe; t = p.t_xlpe, tag = :insulation)
        screen(semicon; t = 0.8e-3, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = 0.3e-3, tag = :screen)
        end
        jacket(pe; t = 2.5e-3, tag = :jacket)
    end
    d = 2 * r_ex(design.outer)
    placements = [Placement(-d / 2, -1.0), Placement(d / 2, -1.0),
                  Placement(0.0, -1.0 + d * sqrt(3) / 2)]
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = i, sheath = 0) for i in 1:3], id = "cable", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = p.rho_earth, eps_r = 10.0);
        temperature = 20.0, frequencies = [f])
    return positive_sequence(compute(problem, CoaxialFormulation(transposition = true)))
end

function line_zy(p, f)
    wire   = Material(:conductor, p.rho_wire, 1.0, 1.0, 20.0, 0.00403)
    design = @cable "ohl_33kv_al" begin
        @terminal :core begin
            core(wire; r = p.r_wire, tag = :wire)
        end
    end
    placements = [Placement(x, p.height) for x in (-p.spacing, 0.0, p.spacing)]
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = i,) for i in 1:3], id = "line", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = p.rho_earth, eps_r = 10.0);
        temperature = 20.0, frequencies = [f])
    return positive_sequence(compute(problem, CoaxialFormulation(transposition = true)))
end

# ── HarmonicPowerModels side ──────────────────────────────────────────────────

data  = HPM.parse_file(joinpath(@__DIR__, "data", "hello_world.m"))
hdata = HPM.replicate(data; H = [1, 5, 7, 11, 13])
zbase = 33.0^2 / 100.0                                   # Ω, 33 kV on 100 MVA

# the uncertain inputs of the cable type, shared by both cable sections
box = ParameterBox(:eps_xlpe  => (2.3, 2.25, 2.5),
                   :t_xlpe    => (8.0e-3, 0.95 * 8.0e-3, 1.10 * 8.0e-3),
                   :r_core    => (14.16e-3, 0.99 * 14.16e-3, 1.01 * 14.16e-3),
                   :rho_earth => (100.0, 10.0, 1000.0, :log))

cable(p, f, seq) = (seq === :zero && error("no zero-sequence data"); cable_zy(merge(CABLE, p), f))
line(p, f, seq)  = (seq === :zero && error("no zero-sequence data"); line_zy(LINE, f))

# branch order of hello_world.m: 1 source, 2 cable 15 km, 3 cable 10 km, 4 line 10 km
model = ExactPiModel(Dict(2 => (zy = cable, len = 15_000.0, zbase = zbase),
                          3 => (zy = cable, len = 10_000.0, zbase = zbase),
                          4 => (zy = line,  len = 10_000.0, zbase = zbase)))

"Driving-point impedance |Z_nn| at harmonic order `h`, clean buses eliminated."
function driving_point(hdata, model, p, h, n)
    Y, idx = admittance_at(hdata, model, p, h)
    keep   = [idx[b["bus_i"]] for b in values(hdata["nw"]["1"]["bus"]) if !get(b, "is_clean", false)]
    sort!(keep)
    Z = inv(Y[keep, keep])
    return abs(Z[findfirst(==(idx[n]), keep), findfirst(==(idx[n]), keep)])
end

hs = collect(1.0:0.02:25.0)
scan(p) = [driving_point(hdata, model, p, h, 3) for h in hs]
peak(z) = (k = argmax(z); (h = hs[k], z = z[k]))

z_nom = scan(nominal_point(box))
p_nom = peak(z_nom)
@printf("nominal: resonance at h = %.2f, |Z_33| = %.1f pu\n", p_nom.h, p_nom.z)

vtx  = Iterators.product(ntuple(_ -> (-1.0, 1.0), length(box))...)
pks  = [peak(scan(latent(box, collect(v)))) for v in vtx]
@printf("over the %d vertices: resonance between h = %.2f and %.2f\n",
        length(pks), minimum(pk.h for pk in pks), maximum(pk.h for pk in pks))
