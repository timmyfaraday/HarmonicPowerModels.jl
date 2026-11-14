# LineCableModels

HarmonicPowerModels (HPM) joins a harmonic network to physical cable and overhead-line
inputs through its LineCableModels (LCM) package extension. Loading both packages loads
`HarmonicPowerModelsLineCableModelsExt`. HPM implements the network-model protocol and the exact
π-model. LCM computes series impedance and shunt admittance per meter from materials,
geometry, connections, earth and frequency. The extension reduces these parameters to
sequences and supplies them to HPM.

## Install and load

Use Julia 1.12. In an environment containing this version of HPM, add LCM's dependency
and source records to `Project.toml`, then run `julia --project -e 'using Pkg; Pkg.instantiate()'`:

```toml
[deps]
HarmonicPowerModels = "56365af1-a032-4c1d-861b-3c9f15bdb54e"
LineCableModels = "dffb7669-8cd4-4628-99dd-7694e356ba37"

[sources]
LineCableModels = {url = "https://github.com/Electa-Git/LineCableModels.jl", rev = "e180b435e39495eeb7eb16dbdaa5193d374a9f60"}
```

The examples below run in order from the HPM repository root. Start Julia with
`julia -t 4 --project` for parallel network evaluations. Install `JuMP`, `Clarabel`,
`Ipopt` and `PowerModels` in the example environment as well. Each `@assert` states a value
that the example checks. The extension test kit runs these blocks as one session.

```julia
using HarmonicPowerModels
using LineCableModels
using LinearAlgebra: BLAS
using JuMP: optimizer_with_attributes
using Clarabel, Ipopt, PowerModels
BLAS.set_num_threads(1)
```

## Network models

An `AbstractNetworkModel` implements `uncertain_branches(model)` and
`branch_data(model, b, p, h, sequence)`. The latter returns per-unit
`(br_r, br_x, g_fr, b_fr, g_to, b_to)` at a named parameter point `p`, real harmonic
order `h` and `:zero`, `:positive` or `:negative` sequence. Branches outside the model
retain HPM's harmonic scaling from `replicate`.

`ScaledBranchModel(hdata, classes; lengths = Int[])` scales source resistance,
reactance and charging with `kR_<class>`, `kL_<class>` and `kC_<class>`. Selected lengths
use `len_<branch>`. Missing factors equal one. It supports positive and negative sequences.
Use `ExactPiModel` for physical frequency dependence and zero-sequence parameters.

Without LCM, `ExactPiModel(branches; f1 = 50.0)` takes records
`(zy = ..., len = ..., zbase = ...)`. The protocol is
`zy(p, f, sequence) -> (z, y)`, with frequency in Hz, `z` in Ω/m and `y` in S/m.
For length ℓ in meters, HPM computes

```math
\gamma=\sqrt{zy},\qquad Z_c=\sqrt{z/y},\qquad
Z=Z_c\sinh(\gamma\ell),\qquad Y/2=\tanh(\gamma\ell/2)/Z_c.
```

Series impedance is divided by `zbase` [Ω]. Each end's shunt admittance is multiplied
by it. This is the exact π-equivalent of a uniform line at that frequency.

## Build a cable circuit

This small cable follows `test/linecablemodels/fixtures.jl`, whose physical inputs come
from LCM's `examples/cable_30kv_630mm2_trefoil.jl`. Its dimensions are in meters and
material resistivities in Ω·m. The cable has an aluminum core, XLPE insulation and a
copper screen. The cables touch in trefoil. Each screen has connection index zero,
representing a solid bond and earth connection at both ends.

```julia
function cable_problem(p, f)
    aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.0, 20.0, 0.00403)
    copper = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    semicon = Material(:semicon, 1000.0, 1000.0)
    xlpe = Material(:insulator, 1.0e14, p.eps_xlpe)
    polyethylene = Material(:insulator, 1.0e14, 2.3)
    design = @design "manual_cable" begin
        @terminal :core begin
            core(aluminum; r = 14.16e-3, tag = :core)
        end
        screen(semicon; t = 0.8e-3, tag = :semicon_inner)
        insulation(xlpe; t = 8.0e-3, tag = :insulation)
        screen(semicon; t = 0.8e-3, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = 0.3e-3, tag = :screen)
        end
        jacket(polyethylene; t = 2.5e-3, tag = :jacket)
    end
    spacing = 2r_ex(design.outer)
    placements = trefoil(design;
        spacing, at = Placement(0.0, -1.0 + spacing / (2sqrt(3))))
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = k, sheath = 0) for k in 1:3],
        id = "manual_cable", line_length = 1.0)
    earth = homogeneous(rho = 100.0, eps_r = 10.0)
    return LineParametersProblem(system, earth; temperature = 20.0, frequencies = [f])
end

box = ParameterBox(:eps_xlpe => (2.3, 2.25, 2.5), :length => (1000.0, 990.0, 1010.0))
p = nominal_point(box)
@assert p == (eps_xlpe = 2.3, length = 1000.0)
@assert latent(box, [-1.0, 1.0]) == (eps_xlpe = 2.25, length = 1010.0)
formulation = CoaxialFormulation(transposition = true)
branches = Dict(
    1 => (problem = cable_problem, len = p -> p.length, zbase = 100.0,
          coordinates = (:eps_xlpe,)),
    2 => (problem = cable_problem, len = 1000.0, zbase = 100.0,
          coordinates = (:eps_xlpe,)))
model = ExactPiModel(branches, formulation, box; f1 = 50.0)
zy = model.branches[1].zy
@assert zy === model.branches[2].zy
sequences = zy(p, 250.0)
@assert keys(sequences) == (:zero, :positive, :negative)
@assert sequences.positive == zy(p, 250.0, :positive)
@assert isapprox(sequences.positive.z, sequences.negative.z; rtol = 1e-12)
bd = branch_data(model, 1, p, 5, :negative)
@assert bd.g_fr == bd.g_to && bd.b_fr == bd.b_to
@assert bd.br_r > 0
z50 = 1000zy(p, 50.0, :positive).z # Ω/km
@assert isapprox(z50, 0.053111 + 0.097625im; atol = 1e-6)
```

The 50 Hz impedance reference is the nominal `TOM-630` row of the committed
`examples/rhhc-variants/output/ieee30-tom-scan.txt`, printed in Ω/km.

The LCM constructor takes `(problem, len, zbase)` records and optional `coordinates`.
`problem(p, f)` must return a `LineParametersProblem` for frequency `f`. `len` is a number
or a function of `p`. Branches sharing the same problem function share one circuit,
including its computed sequences. Their lengths and base impedances may differ.

The callable stored as `zy` is private to the extension. Its two-argument call returns
all three sequences. Its three-argument call selects one. It applies the Fortescue
transformation `F M_s Fᴴ`, with `M_s = (M + transpose(M))/2`, to LCM's phase matrices.
For three phases, the diagonal order is zero, positive, negative. Ideal transposition
makes the matrices circulant: the sequence entries are `m_self + 2m_mutual`,
`m_self - m_mutual`, `m_self - m_mutual`. The extension requires three retained phases.
It warns if the transformed off-diagonal entries exceed `1e-4` times the largest diagonal
entry, and keeps the diagonal.

## Parameter boxes and stored results

`ParameterBox` stores `(nominal, lower, upper[, scale])` for each coordinate name.
The default scale is `:linear`. `:log` maps latent coordinates affinely in the logarithm
and needs positive bounds. `latent(box, ξ)` maps `[-1, 1]^d` to a named tuple.
`nominal_point(box)` gives the specified nominal, which can differ from `latent(box, zeros(d))`.

A circuit stores sequences by point and frequency. Declaring `coordinates` restricts
the key to those coordinates. The problem still receives the full point. Declare every
coordinate the problem reads. Unknown names and conflicting declarations for branches
sharing a problem are refused. Without a declaration, the key uses the full point.
Length belongs to HPM's π-model and need not enter a circuit's declaration.

Writes to the shared table use a lock. Concurrent threads can compute the same absent key,
then store the same result. Points of the box's nominal point type, at `Float64`
frequencies, use the table. Other number types, including dual and interval points,
pass through to LCM and compute directly. The input function must preserve these types.

## Robust hosting capacity

The next block uses HPM's three-bus star fixture, with the cable above on its two lines.
With one harmonic and a small dictionary, the solve runs in seconds. `hdata_at` first installs the nominal physical branch parameters,
including the fundamental. The fundamental operating point stays at nominal during
the network uncertainty calculation. `solve_rhhc` updates its input with fundamental power
flow and PCE data, which the posterior check then uses.

```julia
source = HarmonicPowerModels.parse_file("examples/shhc/data/3node_test.m")
hdata = hdata_at(HarmonicPowerModels.replicate(source; H = [1, 5]), model, p)
conic = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
nlp = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
pce = Dict{String, Any}(
    "alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
    "epsilon" => 0.05, "bound" => "gaussian",
    "harmonics" => Dict(5 => (deg2rad(4.75), deg2rad(40.05))))
result = solve_rhhc(hdata, model, box, conic;
    pce_params = deepcopy(pce), hpf_optimizer = nlp,
    branch_derivative = :dual, nbackground = 4, use_seeds = false, posterior = 8)
@assert result["signed_off"]
@assert 0 < result["lower_bound"] <= result["upper_bound"]
check = posterior_check(hdata, model, box, result["solution"]; n = 8, seed = 2)
@assert check.n == 8
@assert check.violated_fraction == 0.0
Y, bus_idx = admittance_at(hdata, model, p, 5; sequence = :negative)
@assert size(Y) == (3, 3) && length(bus_idx) == 3
resonances = modal_scan(hdata, model, p; hmin = 1.0, hmax = 6.0, step = 0.2)
followed = [locate_resonance(hdata, model, p, r) for r in resonances]
@assert length(followed) == length(resonances)
```

`solve_rhhc` alternates a conic master problem and an approximate worst-point oracle.
The master starts at `nominal_point(box)`. The dictionary starts at latent zero,
with vertices for at most ten coordinates
(or sign-predicted corners above ten), resonance seeds when `use_seeds = true`, and
`nbackground` Sobol points. `effort` controls the oracle's constraints, starts and
iterations. `signoff` controls its final search. `tol` and `maxiter` control stopping.
`branch_derivative = :dual` uses ForwardDiff for branch derivatives. `:central` is the
default. `posterior` requests an independent check in the solve result.

The result includes the master `upper_bound`, the scaled `lower_bound` and `solution`,
`master_solution`, `max_violation`, `gap`, `signed_off` and iteration `history`.
With an approximate oracle, the lower bound holds relative to the points evaluated.
A posterior check supplies additional sampled evidence. It does not exhaust the box.

The Gaussian factor is `Φ⁻¹(1-ε)` and assumes a Gaussian constrained quantity.
Cantelli uses `sqrt((1-ε)/ε)` without that distribution assumption. Both are fixed over
the box. RHHC refuses empirical factors, including factors already stored in the PCE
data, because their calibration on a nominal network does not hold at other points.

`modal_scan` finds local minima of modal admittance magnitude over real harmonic orders.
`locate_resonance` follows a mode to another point by bracketing and interpolation.
`admittance_at` returns the full harmonic admittance and bus-index map. Between harmonics,
these calls use positive sequence by default. Choose `sequence` for `admittance_at`
when a different sequence is intended.

## Other solvers at a physical network point

`networks_at` returns harmonic networks, excluding the fundamental, for moment evaluation.
`moments_at` evaluates the moments for RHHC at a point. `hdata_at` instead returns a
deep copy of the whole solver input, including the fundamental LCM parameters and the
appropriate sequence at each harmonic. Each call gives independent mutable solver data. `solve_hopf` accepts the same
physical-network input through its existing solver interface.

```julia
physical = hdata_at(hdata, model, p)
@assert physical["nw"]["1"]["branch"]["1"]["br_r"] ==
    branch_data(model, 1, p, 1, :positive).br_r
stochastic = solve_shhc(hdata_at(hdata, model, p), sHHC_SOC, conic;
    pce_params = deepcopy(pce), hpf_optimizer = nlp)
deterministic = solve_hhc(hdata_at(hdata, model, p), dHHC_SOC, conic, nlp)
powerflow = solve_hpf(hdata_at(hdata, model, p), PowerModels.IVRPowerModel, nlp)
@assert stochastic["objective_value"] > 0
@assert deterministic["objective"] > 0
@assert string(powerflow["termination_status"]) == "LOCALLY_SOLVED"
```

## Zero sequence and converter populations

`harmonic_sequence(h)` chooses positive sequence for `h % 3 == 1`, negative for
`h % 3 == 2`, and zero for triplen harmonics. The network's transformer vector groups,
neutral earthing and star points determine whether those currents have an earth return.
Declare these as described in [data](data.md) and [transformers](xfmr.md).
A floating unit drops out at that harmonic. Units with an earth return remain active.
`"recombination" => false` aligns unit contributions in the worst case, with zero
variance and no angle data required. It supports orders beyond the supplied angle tables.

The converter inputs in `examples/shhc/pce_params.jl` are the Table III diode bridge and
Table IV phase-controlled thyristor bridge. `unit_shapes` maps load IDs to their own
Beta `(alpha, beta)` shape. Missing units use `alpha_hat`, `beta_hat`. Each harmonic's
range can be a shared tuple or a dictionary from unit ID to range, in radians.
Each unit uses the same germ for its angles at every harmonic.

```julia
mixed = deepcopy(pce)
mixed["unit_shapes"] = Dict(1 => (4.0414, 4.7762), 2 => (4.139, 3.659))
mixed["harmonics"] = Dict(5 => Dict(
    1 => (deg2rad(4.75), deg2rad(40.05)),
    2 => (deg2rad(-138.23), deg2rad(-10.87))))
mixed_result = solve_shhc(hdata_at(hdata, model, p), sHHC_SOC, conic;
    pce_params = mixed, hpf_optimizer = nlp)
@assert mixed_result["objective_value"] > 0
@assert harmonic_sequence.([3, 5, 7]) == [:zero, :negative, :positive]
```

For a fixed-network SHHC calculation, empirical factors use sampled quantiles.
With `"empirical_moments" => "pce"` (the default), `sHHC_SOC` recalibrates each factor
on the expansion mean and standard deviation to equate its moment bound to the sampled
quantile. `"sample"` calibrates on sampled moments instead. The SOC result returns
`"lambdas_used"`. Pass those factors to `sHHC_QCQP` for a comparison using the same factors.
The iterative helper starts from a Gaussian solve and recalibrates the
factors at successive hosting decisions. This small sample demonstrates the call.
Choose the sample size and convergence settings for the study being performed.

```julia
calibrated = solve_shhc_calibrated(hdata_at(hdata, model, p), sHHC_SOC, conic;
    pce_params = merge(mixed, Dict("empirical_moments" => "pce")),
    hpf_optimizer = nlp, nsamp = 256, maxiter = 2)
@assert calibrated["iterations"] in 1:2
@assert haskey(calibrated["lambdas_used"], :ihd)
@assert haskey(calibrated["empirical_lambdas"], :q)
```

The [physical-input studies](linecablemodels_studies.md) select populations through case records.

## Computation and limits

At each network point and frequency, LCM constructs and solves each distinct circuit
problem. The extension reduces both phase matrices once and stores all three sequences.
HPM forms each branch's π-model with its own length and base. The table
saves repeated LCM computations for identical designs on different links.
Coordinate keys also reuse a result when unrelated inputs change.

Commit `9c625ad` measured 28.85% fewer stored LCM results with central derivatives and
16.97% fewer with dual derivatives, with identical bounds. Commit `ad7ca74` computes
dual branch derivatives per chunk of at most 32 latent partials and calls each shared
`zy` once per harmonic. Commit `d440521` measured resonance seeding falling from about
50 s to 12–15 s on 16 threads for the 23-coordinate IEEE30 study. Commit `077b2cc`
measured dictionary memory per point at 49 harmonics falling from 30.8 MB to 2.0 MB,
and dictionary scans becoming 15 times faster. The compact dictionary stores only the
harmonics used by each constraint.

The model assumes a balanced network and ideal circuit transposition. Screens are
bonded and earthed at both ends. Cross-bonding and coupling between parallel circuits
are omitted. The infeed acts as an ideal zero-sequence earth. Transformer winding and
earthing models belong to HPM. There is no three-limb core model. LCM supplies lines
and cables only.

## API reference

```@docs
replicate
AbstractNetworkModel
ScaledBranchModel
ExactPiModel
ParameterBox
latent
nominal_point
solve_rhhc
solve_shhc
shhc_moments
fast_violation_gradient
posterior_check
hdata_at
networks_at
admittance_at
moments_at
Resonance
Seed
modal_scan
locate_resonance
resonance_seeds
build_dictionary
harmonic_sequence
empirical_lambdas
solve_shhc_calibrated
HarmonicPowerModels.fairness_weights
```
