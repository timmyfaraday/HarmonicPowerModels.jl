# Physical-input studies

The studies in `examples/rhhc-variants/` connect physical constructions and their bounded
inputs to HPM's network models. Begin with [LineCableModels](linecablemodels.md) for the
extension API. Run the commands on this page from the repository root. Data paths inside
case records are relative to the study directory. Results and their committed output
inventory are in the [study README](https://github.com/amaurigmartins/HarmonicPowerModels.jl/blob/feat/LineCableModelsExt/examples/rhhc-variants/README.md).

The IEEE30 studies use the `inductive_cable33` scenario and 21 harmonic units at the
first 21 loads. Baseline cases select orders 5, 7, 11 and 13, the Gaussian factor and
risk 0.05. Their fundamental networks use nominal LCM branch parameters. Each cable
screen has a solid bond to ground at both terminals, with connection index `0`.
Each overhead earth wire has a bond at every tower and a solid connection to ground,
with index `0`. Every circuit uses `CoaxialFormulation(transposition = true)`.
The studies omit cross-bonding.

To use the other solvers, include `examples/rhhc-variants/builders.jl` and build
`s = study("cases/ieee30-physical.toml")`. Pass
`hdata_at(s.hd, s.model, nominal_point(s.box))` to `solve_hhc`, `solve_hpf`,
`solve_hopf` or `solve_shhc`. Call
`solve_hhc(hdata_at(s.hd, s.model, nominal_point(s.box)), dHHC_SOC, s.conic, s.ipopt)`.
For `solve_shhc`, also pass `pce_params = deepcopy(s.pce)` and `hpf_optimizer = s.ipopt`.
Each call to `hdata_at` returns independent solver data with the fundamental LCM network.

## Files and input ownership

```text
examples/rhhc-variants/
    Project.toml             study environment and pinned package sources
    builders.jl              cases and physical inputs to network models
    data/materials.toml      fixed materials and reference temperature
    data/specs/              primary quantities of each product
    data/families/           shared physical inputs and intervals
    data/formations/         geometry and placement conventions
    data/register/           inputs.csv, derive.py, derived.json, README.md
    cases/                   network assignments and parameter selections
    networks/                saved converted networks and conversion script
    run.jl                   capacity and posterior reports
    compare.jl               paired capacity and corner reports
    scan.jl                  sequence tables and modal scans
    sweep.jl                 physical-input admission study
    output/                  saved reports and sequence tables
    test/                    study suite
```

Specs describe individual cable or conductor products through primary quantities, SI
units and citations. Families define inputs that vary together across products.
Formations define placement and support geometry. Cases select products, formations,
routes, ratings and active box coordinates. The builder derives coupled dimensions and
material quantities from those inputs. The extension supplies LCM results to HPM.

## Physical-input register

This register lists the physical inputs of the IEEE 30-bus cables and overhead lines. Each
input has a nominal value with its source and an interval with its basis. The studies of
`examples/rhhc-variants/` cite its rows by their `id`. An audit checked the register against
licensed copies of the standards it cites. A value enters a spec only after the register holds
it.

| File | Content |
| --- | --- |
| `inputs.csv` | one row per input, with the sixteen columns below |
| `derive.py` | the arithmetic behind every derived value, in the Python standard library |
| `derived.json` | the output of `derive.py` |

Run `python3 -I examples/rhhc-variants/data/register/derive.py` from the repository root to
reproduce `derived.json`.

### Columns

| Column | Content |
| --- | --- |
| `id` | the row's name, which the specs cite |
| `construction` | the spec or the family the row describes |
| `component` | the part of the construction |
| `quantity` | the physical quantity |
| `lcm_input` | the LineCableModels input it sets, or its role as a check value |
| `unit` | the unit of `nominal`, `lower` and `upper` |
| `nominal` | the nominal value |
| `definition_class` | the class of the nominal value |
| `definition_source` | the source of the nominal value |
| `lower`, `upper` | the interval |
| `interval_class` | the class of the interval |
| `interval_basis` | the basis of the interval |
| `sensitivity` | the effect of the input on the harmonic network |
| `coordinate` | `yes` for an established coordinate, `sweep` for an interval awaiting the sensitivity study, `no` for a fixed input, an absorbed quantity or a value outside LCM. Rows outside the full IEEE30 audit retain their earlier labels |
| `notes` | remarks and the history of the value |

### Classes

| Nominal value | Definition |
| --- | --- |
| STD | a value that a standard fixes, such as a tabulated thickness or a maximum resistance |
| CAT | a value in the datasheet of a commercial product, with its page |
| DER | a value that `derive.py` computes from cited input rows, including stated assumptions |
| ASM | an engineering assumption, with its basis |

| Interval | Definition |
| --- | --- |
| N | the acceptance limit of a standard, often one-sided |
| P | the tolerance that the manufacturer prints for the product |
| E | production evidence, the values that type-test reports measured |
| M | the range of a material family, across compounds or products |
| S | a study interval inside the N, P and E limits, with its reason |
| ENV | installation, environment or operation, never a manufacturing tolerance |
| FIX | a discrete choice of construction, or an input that the result does not depend on |
| DER | an interval that follows from other rows |

An asterisk, as in `N*`, `P*` or `CAT*`, marks a value taken by analogy from another product or
standard. It also marks a clause that the audit did not check against a licensed copy.

The candidate study is the physical-candidate study of the IEEE 30-bus case that preceded this
register. It is not part of the repository. A row that names it records a choice or a value of
that study, and the row's class says how far the audit supports it.

### Citations and rows

A spec cites a row ID such as `LINK.length_m`.
`definition_source` identifies the source of the nominal value and `interval_basis`
identifies the basis of its bounds. The register's `README.md` resolves citation keys
such as IEC60228, IEC60502-2, IEC60287-1-1, EN50182, OLEX and APAR to editions,
clauses, catalogue pages and links. Its construction sections give the assignments,
analogies, coupled quantities and adopted environmental intervals. Consult those rows
when adding a product or changing its bounds.

## Case records

Numeric records use `{value = ..., register = "row.id"}` or
`{value = ..., origin = "source and location"}`. The builder decodes their values for
computation. The records retain their sources for readers. Each register row identifies
the nominal definition and interval basis separately. Keep dependent inputs coupled.
The wire diameter sets its area and resistance. Layer thickness sets outside diameter. Outside
diameter sets touching-trefoil spacing.

| Key | Meaning |
| --- | --- |
| `network` | `kind = "scenario"` with `scenario`, or `kind = "matpower"` / `"converted"` with `file`. File paths are relative to the study directory. |
| `model` | `"exact"` uses LCM and `ExactPiModel`. `"scaled"` uses source data and `ScaledBranchModel`. |
| `box_names` | Ordered active box coordinates. Other declared candidates stay nominal. Bounded branch lengths are appended in numeric branch order. |
| `[box.<name>]` | `family` and `coordinate` resolve a family input and its bounds. A direct row supplies `nominal`, `lower`, `upper` and optional `scale`. |
| `[bindings.<family>]` | Maps a physical family field to the name in the parameter point. Unbound inputs use family nominals. |
| `[branches."<index>"]` | Assigns `spec` and `formation` to a source branch. Optional `route`, `length_m` and `rating` select the length and rating. |
| `spec`, `formation` | File stems under `data/specs/` and `data/formations/`. Equal pairs share one problem function. |
| `route` | Name of a `[routes.<name>]` record containing `length_m` and its source. Choose either `route` or branch `length_m`. |
| `length_m` | Fixed value in meters or a record with `value`, `lower`, `upper` and a citation. Bounds create `length_<index>`. |
| `rating` | Name under `[ratings]`, whose `current_A` supplies the fixed branch rating. It is converted with the local voltage and network power base. |
| `classes` | In a scaled case, maps class names to branch indices for R, L and C factors. |
| `harmonics` | Harmonic orders excluding the fundamental. The builder adds order 1. |
| `units` | A count of the first loads or an explicit list of load IDs. |
| `epsilon`, `bound` | Chance-constraint risk and factor. The supplied cases use `"gaussian"`. RHHC also accepts `"cantelli"`. |
| `recombination` | Defaults to `true`. `false` uses worst-case alignment without angle statistics. |
| `population` | `"diode"`, `"thyristor"` or `"mix"`. Requires recombination and the supplied table harmonics. |
| `base` | One case stem whose top-level records are replaced by this variant's records. Nested bases are refused. |
| `f1` | Fundamental frequency in Hz. |
| `posterior` | Number of posterior points used by the runner. |
| `principle` | Optional fairness principle replacing the source value. RHHC accepts absolute or weighted equality. |
| `note` | Description of the case's physical assumptions. |
| `star_points` | Bus IDs representing star points of three-winding transformer equivalents. |
| `[transformers."<index>"]` | `f_bus`, `t_bus`, `vg`, `gnd1`, `gnd2`, `re1_ohm`, `xe1_ohm`, `re2_ohm`, `xe2_ohm`. The builder verifies branch identity and converts neutral impedances on each winding's voltage base. |
| `[generators."<index>"]` | `zero_sequence_open = true`, or `rsc0` and `xsc0` in per unit, supplies the generator zero-sequence model. |

This complete case uses the supplied three-bus star and the fixed nominal U33-630
cable construction. Its one active coordinate is XLPE permittivity. Each branch has
an example length of 1000 m and keeps its source network rating. The manual tests build this record.

```toml
# complete-case
model = "exact"
bound = "gaussian"
network = {kind = "matpower", file = "../shhc/data/3node_test.m"}
box_names = ["eps_U33"]
harmonics = {value = [5], origin = "examples/shhc/pce_params.jl, Table III"}
units = {value = [1, 2], origin = "examples/shhc/data/3node_test.m, load IDs"}
epsilon = {value = 0.05, origin = "examples/shhc/pce_params.jl"}
f1 = {value = 50.0, origin = "examples/shhc/ieee30_scenarios.jl"}
posterior = {value = 8, origin = "manual example sampling choice"}

[box.eps_U33]
family = "U33"
coordinate = "eps"

[bindings.U33]
eps = "eps_U33"

[branches."1"]
spec = "U33-630"
formation = "trefoil"
length_m = {value = 1000.0, origin = "manual example route"}
origin = "manual example assignment"

[branches."2"]
spec = "U33-630"
formation = "trefoil"
length_m = {value = 1000.0, origin = "manual example route"}
origin = "manual example assignment"
```

Save this as a TOML file and call `study(path)` after including
`examples/rhhc-variants/builders.jl`. It returns `file`, decoded `config`, materials
`mat`, `families`, `box`, `model`, branch `details`, harmonic data `hd`, converter
parameters `pce`, the `ipopt` and `conic` optimizers, and distinct LCM callables `lcm`.
`study(path; all_coordinates = true)` activates every declared candidate for sweeps.
The exact model writes the nominal fundamental branch data before replication.

The following fragment shows a length declaration. Its bounds are an example
installation choice. A study must cite its own registered interval.

```toml
[branches."1"]
spec = "U33-630"
formation = "trefoil"
length_m = {value = 1000.0, lower = 990.0, upper = 1010.0, origin = "manual example installation interval"}
```

## Converter populations and case variants

A case selects its converter population with a cited record:

```toml
population = {value = "mix", origin = "examples/shhc/paper_mix_ab.jl"}
```

`diode` and an omitted key use the existing Table III parameters, without a
`unit_shapes` entry. `thyristor` uses the Table IV Beta shape and angle ranges from
`PCE_PARAMS_THYRISTOR` in `examples/shhc/pce_params.jl`. `mix` sorts the harmonic
load ids and assigns diode bridges to odd positions and thyristor bridges to even positions. The 21-unit IEEE30
population has 11 diode bridges and 10 thyristor bridges. Each unit receives its
own angle range at each harmonic and its own Beta shape.

The builder refuses a population with `recombination = false` or a harmonic order
outside the Table III and IV data (5, 7, 11 and 13). It also refuses an unknown
population. Run studies with `SHHC_UNIT` unset or equal to `diode`. Other values are refused
when building a study. `examples/shhc/pce_params.jl` merges the thyristor inputs into
`PCE_PARAMS` when `SHHC_UNIT` is `thyristor`.

The 20 cases with recombination have `-mix` variants. Each variant contains only its
base case name and the population record:

```toml
base = "ieee30-physical"
population = {value = "mix", origin = "examples/shhc/paper_mix_ab.jl"}
```

`base` names a TOML case in `cases/`, without its path or extension. The builder
reads that file and replaces its top-level keys with the variant's keys. A record
is replaced as a whole. A base case that itself declares `base` is refused.
The reports use the variant's own filename. Its box and network model
come from the base case. The population record changes the converter inputs.

The five `ieee30-zs-*` cases use worst-case alignment without recombination and
have no population variants. `compare.jl` accepts the paired variants through the
same two paths as their bases:

```sh
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/compare.jl \
    cases/ieee9-physical-mix.toml cases/ieee9-mix.toml
```

The variants write capacity and comparison reports. Their sequences and physical
input sweeps are those of their bases. The README inventory lists those reports with their bases.

## Gaussian factor

Every case selects `bound = "gaussian"`. HPM uses
`λ = Φ⁻¹(1 - ε)` in the moment constraint `E[S] + λ sqrt(Var[S]) ≤ limit`.
This assumes a Gaussian distribution for the constrained quantity.
Cantelli instead uses `sqrt((1 - ε)/ε)`.
With `recombination = false`, the zero-sequence variants align contributions
in the worst case. Their variance is zero. The factor changes only the label.

## Dataset and construction

The LCM cases use one construction chain: `problem` resolves the family coordinates,
`arrangement` selects the formation's physical state, `design` builds the LCM conductors,
and `system` places and connects them. The bindings select the coordinates of a case.
A family input without a binding, or without a coordinate in the box, has its nominal value.
TOM, U33, ACSR and the earthwire use this same chain.

Family files own the shared varying inputs and their nominal values and intervals.
`EARTH.toml` owns earth resistivity, permittivity and the logarithmic resistivity scale.
`U33.toml` owns cable temperature, burial depth, screen resistance and material properties.
`O132.toml` owns ambient temperature, reference sags, earthwire dimensions and mechanics.
`ACSR.toml` owns wire materials, temperatures and lay ratios. `TOM.toml` defines the placeholder
radius, insulation thickness and permittivity. Materials contain fixed properties and the
common electrical reference temperature. Formation files contain fixed coordinates and
placement conventions. A formation with `structure = "O132"` selects that structure's
mechanics. O33-H has fixed sag and uses the ACSR phase temperature.

Each family field specifies a physical quantity. Case coordinates use the form
`<field>_<FAMILY>`, including `eps_xlpe_TOM`, `t_xlpe_TOM` and `r_core_TOM`.
A per-size spec contains each nominal and its registered physical bounds.
`interpolate` maps the signed normalized coordinate to this interval, with nominal at
zero and the registered endpoints at minus one and plus one. A one-sided interval uses [-1, 0] or [0, 1].
A family row with `scale = "registered"` marks this normalized convention. HPM
samples these normalized positions on its linear scale. This applies to U33 resistance, jacket thickness, conductor diameter, inner semiconductor
thickness and ACSR wire diameters. The box names `t_ins_U33`, `t_sc2_U33` and `d_wire_U33`
bind the shared physical thicknesses and screen-wire diameter in meters. All bounds and
aggregate projection limits come from the family or spec records.

The compacted U33 core has radius equal to half the resolved conductor diameter and
resistivity `R20 * π * radius^2`, with the resolved R20. The copper screen tube has area
`rho_Cu / R20_screen`. Its jacket correction preserves the physical outside diameter.
The TOM design retains its fixed metal resistivities and its family geometry.
LCM computes ACSR resistance from the realized wires. The aggregate projection and
thermal sag model use the same realized mass and metal area.

The family files define normalized jacket and diameter coordinates, held at nominal
when absent from the box.

The cases supply series impedance and shunt admittance per meter through the extension.
HPM computes the exact π-model in `branch_data` of `ExactPiModel`:

```math
Z = Z_c \sinh(\gamma\ell), \qquad Y/2 = \tanh(\gamma\ell/2)/Z_c,
\qquad \gamma = \sqrt{zy}, \qquad Z_c = \sqrt{z/y}.
```

## Coordinate declarations

The builders declare each LCM branch's `coordinates` from its spec family,
formation structure and earth family bindings, restricted to names in the box.
Branches sharing a problem function share one circuit and must declare the same tuple.
The constructor refuses unknown names. Without a declaration, the full point is the key.
The projection is set in the circuit type at construction, and the stored key type is
concrete. The problem still receives the full point. Other point types, including dual
points, compute directly without entering the table.

The study test moves each undeclared coordinate alone to latent positions -1, -0.5,
0.5 and 1 at 50 and 650 Hz. All other coordinates retain their nominal values.
It covers each base case once and skips cases that declare `base`.
It bypasses the table to check the problem's actual dependencies and checks exact
sequence equality with the projected-table result.

## Installed length

A branch's `length_m = {value = ..., register = "LINK.length_m"}` sets its fixed
length in meters. Add `lower` and `upper` to make the length uncertain. The builder
appends the linear coordinate `length_<branch index>` after `box_names`, in numeric
branch order. Without `length_m`, the branch uses its route length. Each branch
chooses either `route` or `length_m`.

`ieee30-physical-length` sets this key on its 22 cable links, with the retained route
lengths as nominal values. Register row `LINK.length_m` supplies the adopted ±1 percent
installation interval (class ENV), based on route deviations, snaking, and slack at
joints and terminations. Length coordinates enter HPM's π-model and stay outside the
LCM circuit's coordinate declaration.

## Transformers

The declarations follow Tom Van Acker's proposal. Branches 12, 15 and 36 are YNyn0,
132/33 kV.

The three-winding unit at bus 6 has star point 9. Branch 11 represents its 132 kV YN winding.
Branch 14 represents its 33 kV yn winding. Branch 13 represents its 11 kV delta winding. Branch 16 is the YNd11 step-up of the
condenser at bus 13, with the delta on the 11 kV side.

The 132 kV neutrals have solid connections to earth. Each 33 kV neutral has a resistance of 19.05 Ω: assumption to agree with Tom:
limits a 33 kV earth fault to 1.0 kA.

The case files store neutral impedances in ohms. The runner prints HPM's `vg`, `gnd1`,
`gnd2`, `re1`, `xe1`, `re2` and `xe2`. It converts the impedances to per unit on their winding
voltage bases. HPM writes YNyn0 as `Yy0` and YNd11 as `Yd11`, with the grounding flags.
The runner checks each branch's existence, transformer status and buses.

The reports print these declarations. HPM applies the vector groups to tapped branches
at every harmonic and the winding earthing in the zero sequence. The network retains
the leakage parameters of the source transformers.

## Zero-sequence variants

The `ieee30-zs-*` cases extend `ieee30-physical` to orders 2 through 50, including
zero-sequence orders 3, 6, …, 48. With `recombination = false`, the contributions of
all 21 units align in the worst case at every bus and branch. The angle ranges do not
enter this calculation. Bus 9 is the star point of the three-winding unit.

Tom Van Acker assumes a connection from each generator through a YNd step-up, with 0.12 pu
leakage on a rating of 0.2 times its short-circuit power and X/R = 40. Its network bus
sees an impedance to earth. The condensers at buses 11 and 13 are open in
the zero sequence, behind branches 13 and 16. The 19.05 Ω neutral resistance is also
Tom's assumption, corresponding to a 1 kA earth fault at 33 kV. The case records state
these origins. The variants change the windings and earthing of units with a 33 kV winding.

| Case | 132/33 kV units | 33 kV earthing |
| --- | --- | --- |
| `ieee30-zs-ynyn-solid` | YNyn0 | Solid |
| `ieee30-zs-ynyn-resistance` | YNyn0 | 19.05 Ω |
| `ieee30-zs-dyn` | Dyn11, delta on 132 kV | 19.05 Ω |
| `ieee30-zs-ynd` | YNd11, delta on 33 kV | Condenser step-up star unearthed |
| `ieee30-zs-unearthed` | YNy0 | Every 33 kV star unearthed |

The clock numbers preserve the phase shift along every 132/33 kV path. The positive-
and negative-sequence capacities stay equal across these variants. The run
report gives the nominal hosting capacity and robust hosting capacity, each summed per
sequence, and the bound used.

The zero-sequence model keeps each lumped earthing resistance constant with harmonic
order. The earthing reactance scales with frequency. It refers the from-side neutral
term to the leakage side as `3Zn1 / abs2(tap)`, while the to-side term is `3Zn2`.
Star-to-star clocks 2, 6 and 10 reverse zero-sequence polarity by 180 degrees.
The other star-to-star clocks retain positive polarity. In a floating connected part,
harmonic units drop out of the injections and allocation weights at that harmonic.
Units with an earth return retain their capacity. The reports use this corrected model.

## Run

From the repository root, with Julia 1.12 and the public LCM source available:

```sh
julia --project=examples/rhhc-variants -e 'using Pkg; Pkg.instantiate()'
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/ieee30-tom.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/ieee30-u33.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/ieee30-physical.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl cases/ieee30-tom.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl cases/ieee30-u33.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl cases/ieee30-physical.toml
```

Paths such as `cases/ieee30-tom.toml` are relative to this directory unless they already
identify a file on disk. The solver is Clarabel. Ipopt solves the fundamental power flow.
Both project files pin the public source
`https://github.com/Electa-Git/LineCableModels.jl` at commit
`e180b435e39495eeb7eb16dbdaa5193d374a9f60`. `Pkg.instantiate()` fetches this
commit. Pkg supplies the pinned LCM source.
HPM resolves from this repository through the relative source path `../..`.
For a fresh kit environment, run
`julia --project=test/linecablemodels -e 'using Pkg; Pkg.instantiate()'`.
The reports read the resolved LCM revision from the study manifest.
The scripts include the shared PCE parameters and scenario builders under `examples/shhc/`.

`run.jl` writes `output/<case>.txt`. It computes the nominal stochastic capacity with a
zero-width box, scans the nominal network, then calls `solve_rhhc` over the parameter box
with 1024 posterior points.

The nominal solve evaluates the LCM frequency dependence at
every harmonic. Each report has one header with a repository-relative case path. Solver banners and progress
logs go to the terminal. The report gives the git commit, case, threads, solver, declarations,
capacities, bounds, sign-off, posterior check and sum of the circuits' table sizes.
The commit identifies the repository base. The study and extension can contain uncommitted
changes.

The table sizes count stored point-frequency pairs. Concurrent calls can compute
the same absent pair before either stores it.

`scan.jl` writes `output/<case>-sequences.csv` and the text companion report
`output/<case>-scan.txt`. The CSV has real and imaginary parts of `z0`, `z1`, `y0` and `y1`, per meter in SI units, for every circuit.
It follows HPM's modal scan: orders `1:0.02:(1.25 * maximum(harmonics))`, plus the harmonic
orders. The text report gives the branch-to-circuit table, nominal resonances, sequences
at 50 Hz and data checks. Each branch table also gives the source reactance and scenario
capacitance on the retained route. The inductive scenario sets overhead capacitance to zero.

## Add a spec or case

Add a spec with its family, primary quantities, SI units and a source for every numeric
value.

Use the audited register to source physical inputs. Select a formation
or add sourced placement data. A new case lists its branches, formation, family bindings, ordered box names and
transformer declarations. The family records define its environment and temperatures. A family coordinate that
the case does not bind stays at its nominal value. The U33 family defines the normalized jacket and diameter
coordinates. A case without those coordinates uses their nominal values.

## Other test systems

The case selects `network.kind`: `scenario` calls the IEEE30 scenario builder,
`matpower` reads the declared MATPOWER file, and `converted` reads a saved PowerModels
JSON file. Source paths are relative to this directory. The `model` field selects
`exact` for LCM circuits or `scaled` for source branches with the case's declared factors.
`units` is a count of the first loads, or an explicit list of load IDs.
The case also declares its harmonic orders. Shared angle distributions come from
`examples/shhc/pce_params.jl` for those orders.

LCM branch bases come from the source bus voltage and network power base. Lengths come
from the IEEE30 scenario routes, from route comments in MATPOWER branch rows, or from
`length_m` in converted branch records. Physical lengths use meters.

| Case | Network source | Model | Box coordinates |
| --- | --- | --- | --- |
| `ieee33` | `../shhc/data/33bus_network_shhc.m` | Source branches with R and L factors | `kR_ohl`, `kL_ohl` |
| `ieee33-physical` | As `ieee33` | Six ACSR sizes on an unshielded O13-F pole, fitted routes and LCM charging | Shared physical inputs (11 coordinates) |
| `industrial` | `../shhc/data/industrial_network_shhc.m` | Source cables with R, L and C factors | `kR_cable`, `kL_cable`, `kC_cable` |
| `antwerp` | `../shhc/data/medium_network_shhc.m` | Source cables with R, L and C factors | `kR_cable`, `kL_cable`, `kC_cable` |
| `hello-world` | `../rhhc/data/hello_world.m` | U33-630 on branches 2 and 3, Goat on O33-H on branch 4 | `R20_U33`, `eps_U33`, `t_ins_U33`, `rho_EARTH` |
| `star3` | `../shhc/data/3node_test.m` | Source branches with R and L factors | `kR_ohl`, `kL_ohl` |
| `star5` | `../shhc/data/5node_test.m` | Source branches with R and L factors | `kR_ohl`, `kL_ohl` |
| `ieee9` | `networks/ieee9.json` | Converted source lines with R and L factors | `kR_ohl`, `kL_ohl` |
| `ieee39` | `networks/ieee39.json` | Converted source lines with R and L factors | `kR_ohl`, `kL_ohl` |

The physical variants use the same source networks with the constructions below.
Their intervals, route origins and fixed ratings are recorded in `data/register/README.md`.

| Physical case | LCM construction | Selected coordinates |
| --- | --- | --- |
| `industrial-physical` | U36-Al400 and U11-Al240 screened MV cables, plus the NYY-Cu240-120 variant as three unscreened single cores in touching trefoil | 30 |
| `antwerp-physical` | U155-Cu1600 with a round core of equal DC resistance and screens bonded at both ends | 12 |
| `star3-physical`, `star5-physical` | Falcon ACSR on O132-D, with identical physical branches in each star | 11 |
| `ieee9-physical` | Cardinal ACSR on O230-F, with the saved 230 kV routes | 10 |
| `ieee39-physical` | Cardinal ACSR on O345-F, with the saved 345 kV routes | 10 |
| `ieee33-physical` | Rabbit, Weasel, Dog, Tiger, Wolf and Bison ACSR on O13-F | 11 |

Antwerp retains the source 300 MVA branch limit. Its source construction uses a
Milliken core and cross-bonded screens. The study variant uses a round DC equivalent
and both-end bonds. The industrial LV variant omits the multicore neutral.

The factor ranges follow `examples/rhhc/ieee30_robust.jl:58-60`. They describe
network parameters. The physical cases and their input origins are documented in `data/register/README.md`. Hello-world uses the audited U33 and ACSR families,
with the U33 coordinates and earth interval of `ieee30-u33`. Its overhead family
coordinates stay at nominal values. O33-H has no earth wire.

Run a new case with the same environment and solvers:

```sh
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/hello-world.toml
for case in ieee33 industrial antwerp hello-world star3 star5 ieee9 ieee39; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl "cases/$case.toml"
done
```

The hello-world scan writes the LCM sequence table and modal resonances.
The source cases above use scaled branch data. The capacity runner also scans the
nominal network for cases that use source branches. A failed run writes its exception message to
`output/<case>.txt`.

The source network data retain their voltage limits, ratings and transformer data.
The study premise applies to the cases, with source electrical data kept while physical
specs await audit. The HV candidate ratings assume cross-bonding. These cases instead
assume solid grounding at both cable ends and ideal transposition, and omit cross-bonding.

## Conversion from PowerImpedance

`networks/convert.jl` reads the passive-network portions of the two source examples
at PowerImpedance commit `c0547d03670b2ca4432459274c54b879d345d07b`. It calls
`convert(network, ACDCPowerFlow().backend)` to obtain PowerModels data. It then adds
HPM's ideal-source flags, the study harmonic limits and the harmonic-unit load records.
The saved files include their commit, script and component mapping. A read of each saved
file checks the bus, branch, generator, shunt and load tables against the conversion.

To reproduce the files once, with the existing source environment available:

```sh
julia --project=<PowerImpedance checkout> examples/rhhc-variants/networks/convert.jl <PowerImpedance checkout>
```

The capacity and scan scripts read the saved files without PowerImpedance. The
conversion keeps the passive IEEE9 network before its MMC study and the IEEE39
network before its STATCOM study. It includes the original lines, transformers,
source impedances and capacitors. NetworkBuilder assigns bus numbers. The saved
`conversion_buses` and `conversion_elements` tables identify the original names.

Grounded impedances with active power become HPM loads at their 50 Hz equivalents.
Pure capacitors remain shunts. HPM then applies its own harmonic scaling to the saved
network. This represents a test of the HPM pipeline with the source network at 50 Hz.
It does not reproduce PowerImpedance's frequency-dependent load or transformer models.
The source conversion's branch ratings are defaults, except for `Zg39` as stated below.
Their origins are recorded in `networks/convert.jl`. They are not audited thermal ratings.

The conversion adds the absolute-equality principle and the IEC harmonic limits
used in `examples/shhc/data/3node_test.m:17,33-35`. These are study assumptions.
The JSON records use HPM's native per-unit values on the 1000 MVA source power
base. The conversion uses local voltage bases in line-to-line kV.
It converts each voltage level to its own base, with transformer taps as off-
nominal ratios. Physical route lengths use meters. Cases without a `principle`
declaration keep their source fairness principle.
Hello-world declares absolute equality for `solve_rhhc`.

The conversion reports “Transformer has no open-circuit voltage ratings” twelve times,
from `PowerImpedance src/Network/Components/transformer/transformer.jl:351`.

IEEE9 uses local voltage bases at 13.8, 16.5, 18 and 230 kV. Its transformers have
off-nominal taps near one. The conversion checks the admittance identity under the
change of voltage bases in `networks/convert.jl, normalize`, as well as the saved tables.

IEEE39 source impedance `Zg39` uses the apparent-power capability of generator 39 in
[MATPOWER `case39.m`, tag 8.1, line 136](https://github.com/MATPOWER/matpower/blob/8.1/data/case39.m#L136).
That row gives `Pmax = 1100 MW`, `Qmax = 300 MVAr` and `Qmin = -100 MVAr`.
The capability is `S = sqrt(Pmax^2 + max(Qmax^2, Qmin^2)) = 1140.175425 MVA`.
On the 1000 MVA study base, `rate_a = rate_b = rate_c = 1.140175425` pu.
The nominal-voltage current rating `c_rating_a` has the same per-unit value.
The row's `mBase` is the machine base. It does not set this rating.
The other branches keep their source conversion ratings.

## Full-uncertainty IEEE30 case

`cases/ieee30-full.toml` uses the physical case's 22 U33 cables and 12 ACSR lines.
The cable screens and the earthwire have solid connections to ground, and the circuits use ideal
transposition. Each input has one coordinate per family, shared by its branches.
The case declares 39 wire and material candidates under `box`. The sweep admits 23 into
`box_names`. The separate ACSR resistance factor has been removed.
With more than 10 coordinates, HPM omits the vertices from its oracle dictionary.

The full case selects the ACSR wire inputs. LCM computes their resistance, and the
catalogue R20 remains a check value. The derived resistance envelope is absorbed by
those wire inputs and has no additional coordinate.

EN 50182:2001 5.4 limits conductor diameter to ±1 percent at these sizes, and 5.7.3
limits total mass to ±2 percent. The mass includes the helical wire-length increments.
For a requested point, the builder first applies 5.5.6 with
`lay_outer_eff = min(lay_outer, lay_inner)`. The inner aluminum range ends at 16 and
the steel range starts at 16. Their ordering cannot bind within the registered ranges
except at equality.

A feasible wire realization maps to itself. Otherwise, the builder finds the largest
t in [0, 1] for which both aggregate limits hold along
`d_eff = d_nom + t * (d - d_nom)`. Bisection reduces both diameter deviations along
the same direction. The realized mass and metal area also enter the sag equation.
The input coordinates retain their rectangular bounds. This mapping couples the
physical wires inside the builder, so HPM can keep its parameter box without an
additional constrained-domain protocol.

The U33 coordinates describe core resistance, conductor diameter, insulation thickness,
permittivity, jacket thickness, conductor temperature, screen resistance, both semiconductor
thicknesses and semiconductor permittivity. ACSR coordinates describe
phase temperature, aluminum and steel wire dimensions, their resistivities, steel
permeability and the two aluminum lay ratios. The O132 coordinates describe ambient
temperature and earthwire diameter, resistivity and permeability. Earth resistivity has
one coordinate. The normalized per-size coordinates interpolate the individual products' registered
limits. Shared physical inputs take their values directly from the family records.

LCM's public `LineParametersProblem` provides one temperature for every component.
The full-case builder evaluates each metal's resistivity using its own temperature and
registered electrical temperature coefficient. It passes the reference temperature of
20 °C to LCM, which applies no further correction. The cable core and screen use the
U33 temperature. ACSR phases use the maximum of their temperature and ambient, while
the earthwire uses ambient. The other cases use the same temperature evaluation with their own nominal family inputs.

The phase and earthwire sag follow the register's elastic thermal state equation at each
evaluation, with a fixed 300 m span. The positive sag root uses bisection in the input
number type. Phase mass follows the catalogue component masses, wire areas and helical
lengths. The earthwire uses its steel density and helical wire lengths. The lay ratio of the outer aluminum layer
is capped by the inner ratio, as EN 50182 requires. The sag calculation preserves `BigFloat`. The suite checks endpoints against `derived.json`.

Each sweep moves one candidate to each endpoint with all others nominal. It evaluates
every affected circuit at 50 and 650 Hz through the extension's `zy(p, f)` protocol.
The maximum of |Δ|/|nominal| over z, y and all three sequences admits the coordinate
when it reaches 0.001. This rule applies equally to the previously recommended inputs.
The sweep table records the endpoint and result that give that maximum. The complete report
is `output/ieee30-full-sweep.txt`, and the register records each decision.

The scan binds the thermal coordinate and evaluates a point at the electrical reference
temperature or the family's registered upper temperature. It leaves family nominals unchanged.
Strand spans, steel-layer counts and earthwire identity come from the specs.

## Study tests

From the repository root, run the physical model, corner adversary and reporting tests:

```sh
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/test/runtests.jl
```

`compare.jl` compares capacity runs and evaluates their solutions with random posterior
points and corners predicted from nominal constraint sensitivities. It also accepts a
single source case, such as `cases/ieee30-full.toml`, for its own-case corner check.
`scan.jl` and `sweep.jl` evaluate the selected case's sequences and physical inputs.

## Reproduce every output

The commands below run from the repository root after instantiation. They use the
same case settings as the saved reports. Each comparison also writes the capacity
reports for its two cases. A scan writes both its text report and sequence CSV.
A sweep writes the complete endpoint admission table.

```sh
for case in ieee30-tom ieee30-u33 ieee30-physical ieee30-full hello-world ieee30-physical-length \
    ieee30-zs-dyn ieee30-zs-unearthed ieee30-zs-ynd ieee30-zs-ynyn-resistance ieee30-zs-ynyn-solid; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl "cases/$case.toml"
done
for case in antwerp ieee33 ieee39 ieee9 industrial star3 star5; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/compare.jl \
        "cases/$case-physical.toml" "cases/$case.toml"
done
for case in ieee30-tom ieee30-u33 ieee30-physical ieee30-full hello-world ieee30-physical-length; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl "cases/$case-mix.toml"
done
for case in antwerp ieee33 ieee39 ieee9 industrial star3 star5; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/compare.jl \
        "cases/$case-physical-mix.toml" "cases/$case-mix.toml"
done
for case in ieee30-tom ieee30-u33 ieee30-physical ieee30-full hello-world ieee30-physical-length \
    antwerp-physical ieee33-physical ieee39-physical ieee9-physical industrial-physical star3-physical star5-physical; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl "cases/$case.toml"
done
for case in ieee30-full antwerp-physical ieee33-physical ieee39-physical ieee9-physical \
    industrial-physical star3-physical star5-physical; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/sweep.jl "cases/$case.toml"
done
```

A zero-sequence IEEE30 case needs about 11 GB of memory.

