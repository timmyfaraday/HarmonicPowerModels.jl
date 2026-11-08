# Network studies with LineCableModels

The IEEE30 studies use the `inductive_cable33` scenario and 21 harmonic units at the first 21 loads.
The baseline cases use harmonic orders 5, 7, 11 and 13, Cantelli bounds and a risk level of 0.05. The fundamental
network uses the nominal LCM branch parameters. HPM computes each branch's exact π-model
from the extension's series impedance and shunt admittance per meter at each frequency.

Each cable screen has a solid bond to ground at both terminals, with connection index `0`.
Each overhead earth wire has a bond at every tower and a solid connection to ground, with index `0`.
Every circuit uses `CoaxialFormulation(transposition = true)`. The studies omit cross-bonding.

| Case | 33 kV cables | 132 kV lines | Box |
| --- | --- | --- | --- |
| `ieee30-tom` | Tom's placeholder `TOM-630`, trefoil in contact, lower cables at 1 m depth | Scenario data, zero shunt capacitance | Tom's permittivity, insulation thickness, core radius and earth resistivity |
| `ieee30-u33` | Four Olex U33 sizes, trefoil in contact, group center at 0.8 m depth | Scenario data, zero shunt capacitance | Normalized family resistance, permittivity, insulation thickness and earth resistivity |
| `ieee30-physical` | As `ieee30-u33` | Lion, Goat, Bear and Panther ACSR on `O132-D`, with one steel earth wire | As `ieee30-u33` |

The cables on branches 20, 31, 33 and 34 use U33-400. Branch 26 uses U33-630 and branch 21
uses U33-800. The other 33 kV cables use U33-500. The case files give every assignment.
Cable lengths follow the scenario's retained routes. The overhead lengths follow
`length_km = x_pu * zbase / OHL_X[132.0]`, with `OHL_X[132.0] = 0.40 Ω/km`.

The equivalent and surrogate levels live on the branch `experimental/surrogates`.

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

The cases supply series impedance and shunt admittance per meter through the extension.
HPM computes the exact π-model in `branch_data` of `ExactPiModel`:

```math
Z = Z_c \sinh(\gamma\ell), \qquad Y/2 = \tanh(\gamma\ell/2)/Z_c,
\qquad \gamma = \sqrt{zy}, \qquad Z_c = \sqrt{z/y}.
```

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

## Run

From the repository root, with Julia and the LCM worktree available:

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
The studies and kit use LCM at `e180b435e`. Resolve their environments against a checkout
of that commit.
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

`scan.jl` writes `output/<case>-sequences.csv` and `output/<case>-scan.txt`. The CSV has real
and imaginary parts of `z0`, `z1`, `y0` and `y1`, per meter in SI units, for every circuit.
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

## Checks and results

Check each U33 nominal overall diameter against the catalogue within 0.05 mm. Compare its
low-frequency resistance with R20, its resistance at 90 °C with Olex's AC value, and its
50 Hz reactance and capacitance with Olex. The hot calculation heats both core and screen.
The low-frequency diagnostic uses 0.01 Hz at 20 °C. It approaches DC and includes the
circuit's bonding and earth-return model.

For ACSR, compare the low-frequency resistance with catalogue R20. Report the 50 Hz
reactance and capacitance beside the case's source values, without a catalogue comparison.
For the earth wire, compare with 3.55 Ω/km and the register interval 3.30–4.01 Ω/km.
The earth-wire diagnostic uses three copies as phases at the same formation heights.

The data checks below come from `output/ieee30-tom-scan.txt`. All resistance and reactance
values are in Ω/km. Capacitance is in µF/km. The listed specs use nominal family coordinates.

| Spec | DC LCM | R20 catalogue | R50 at 20 °C | R50 at 90 °C | Olex R90 |
| --- | --- | --- | --- | --- | --- |
| TOM-630 | 0.044870 | n/a | 0.053111 | n/a | n/a |
| U33-400 | 0.077800 | 0.077800 | 0.089553 | 0.109096 | 0.101000 |
| U33-500 | 0.060500 | 0.060500 | 0.072368 | 0.087005 | 0.079700 |
| U33-630 | 0.046900 | 0.046900 | 0.058964 | 0.069724 | 0.062900 |
| U33-800 | 0.036700 | 0.036700 | 0.049131 | 0.056947 | 0.050700 |

| Spec | X50 LCM | Olex X50 | C50 LCM | Olex C | Diameter LCM/catalogue mm |
| --- | --- | --- | --- | --- | --- |
| TOM-630 | 0.097625 | n/a | 0.298556 | n/a | 53.120 / n/a |
| U33-400 | 0.106870 | 0.116000 | 0.267196 | 0.267000 | 50.900 / 50.900 |
| U33-500 | 0.102364 | 0.111000 | 0.296939 | 0.297000 | 54.600 / 54.600 |
| U33-630 | 0.099162 | 0.108000 | 0.329092 | 0.329000 | 58.600 / 58.600 |
| U33-800 | 0.095515 | 0.104000 | 0.365399 | 0.366000 | 63.400 / 63.400 |

U33 DC resistance agrees with the catalogue input at the printed precision. All nominal
overall diameters agree within 0.05 mm.

Capacitance differs from Olex by at most 0.16%. The catalogue's AC resistance excludes the
screen loss that both-end bonding causes. This loss accounts for LCM's higher R90.
The checks add the circulating-current loss to Olex R90. IEC 60287-1-1:2006, Section 2.3.1 gives:

```math
\Delta R = \frac{R_s X_m^2}{R_s^2 + X_m^2}, \qquad
X_m = \omega \, 2\times10^{-7} \ln(2S/d_s).
```

Here `S` is the overall cable diameter, the spacing in trefoil with cables in contact.
The screen's mean diameter `d_s` equals its inner diameter plus one wire diameter.
The resistance `R_s` uses the screen's R20 and copper temperature coefficient at 90 °C.
That section gives `X_m` and the loss. The reactance below is the classical inductance of a
solid round conductor in trefoil, with an internal part of 0.05 mH/km. The reduction follows
from the same loop of core and screen:

```math
X_{\mathrm{classical}} = \omega [0.05 + 0.2\ln(2S/d_c)]\times10^{-3}
\quad [\Omega/\mathrm{km}], \qquad
\Delta X = \frac{X_m^3}{R_s^2 + X_m^2} \quad [\Omega/\mathrm{m}].
```

LCM's X agrees with the classical formula after the screen-current reduction.
The catalogue does not state the basis of Olex's reactance. Its 6–8% excess over the
classical values remains unexplained.

The corrected comparisons below use `R_s` at 90 °C. The LCM resistance uses 90 °C.
The LCM reactance uses the cases' nominal temperature of 20 °C. R and X are in Ω/km.

| Spec | Corrected R90 | LCM R90 | Difference % | Corrected X50 | LCM X50 | Difference % |
| --- | --- | --- | --- | --- | --- | --- |
| U33-400 | 0.109148 | 0.109096 | -0.047 | 0.107886 | 0.106870 | -0.941 |
| U33-500 | 0.087723 | 0.087005 | -0.819 | 0.103424 | 0.102364 | -1.025 |
| U33-630 | 0.070802 | 0.069724 | -1.523 | 0.100310 | 0.099162 | -1.144 |
| U33-800 | 0.058565 | 0.056947 | -2.763 | 0.096823 | 0.095515 | -1.351 |

| Spec | DC LCM | Catalogue R20 | R50 LCM at 20 °C | X50 LCM | Source X50 | C50 LCM | Scenario C50 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Panther | 0.132840 | 0.136300 | 0.133342 | 0.400219 | 0.400000 | 0.009062 | 0 |
| Lion | 0.118227 | 0.121200 | 0.118780 | 0.396553 | 0.400000 | 0.009149 | 0 |
| Bear | 0.106532 | 0.109300 | 0.107136 | 0.393275 | 0.400000 | 0.009228 | 0 |
| Goat | 0.086861 | 0.089100 | 0.087580 | 0.386846 | 0.400000 | 0.009386 | 0 |

ACSR DC resistance is 2.5% below catalogue R20. BS EN 50182:2001, Section 5.10 bases nominal
DC resistance on aluminum wires and aluminum-clad steel. It excludes galvanized steel.
LCM's terminal also conducts in its galvanized steel core. The checks print the steel's
share of DC conductance from wire areas and resistivities:

```math
q_{\mathrm{st}} = \frac{A_{\mathrm{st}}/\rho_{\mathrm{st}}}
{A_{\mathrm{st}}/\rho_{\mathrm{st}} + A_{\mathrm{al}}/\rho_{\mathrm{al}}}.
```

The steel's conductance accounts for LCM's lower DC resistance. The diagnostic uses
wire areas without lay corrections. The steel share is
3.328002% for Panther, Lion, Bear and Goat.

The earth wire gives 3.545151 Ω/km against 3.55 Ω/km, within
3.30–4.01 Ω/km. The source X50 column uses the scenario's 0.40 Ω/km route convention.
The scenario's overhead capacitance is zero. These are case comparisons for X and C.

## Nominal scans

Each circuit has 763 frequency rows from 50 to 812 Hz, with the harmonic orders included.
The following resonances come from HPM's modal scan.

| Case | Circuits | Sequence rows | First resonance h | Other resonances h |
| --- | --- | --- | --- | --- |
| ieee30-tom | 1 | 763 | 6.8574 | none in the scan |
| ieee30-u33 | 4 | 3052 | 6.9504 | none in the scan |
| ieee30-physical | 8 | 6104 | 6.2720 | 12.3169 |

At 50 Hz, the sequence impedances below are in Ω/km. The cable circuits appear in both
U33 cases. The overhead circuits appear in `ieee30-physical`.

| Circuit | z0 | z1 |
| --- | --- | --- |
| TOM-630/trefoil-tom | 0.409113 + j0.119907 | 0.053111 + j0.097625 |
| Lion/O132-D | 0.318170 + j1.348817 | 0.118780 + j0.396553 |
| Goat/O132-D | 0.286970 + j1.339110 | 0.087580 + j0.386846 |
| Bear/O132-D | 0.306526 + j1.345539 | 0.107136 + j0.393275 |
| Panther/O132-D | 0.332732 + j1.352484 | 0.133342 + j0.400219 |
| U33-500/trefoil | 0.317285 + j0.085548 | 0.072368 + j0.102364 |
| U33-400/trefoil | 0.334408 + j0.089387 | 0.089553 + j0.106870 |
| U33-800/trefoil | 0.294047 + j0.079804 | 0.049131 + j0.095515 |
| U33-630/trefoil | 0.303939 + j0.083011 | 0.058964 + j0.099162 |

## Capacity runs

The reports below use 16 Julia threads and Clarabel. Each posterior check uses 1024 points.

| Case | Nominal pu | `lower_bound` pu | Upper bound pu | Gap | Iterations | Time s | Sign-off | Posterior violated | Stored pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ieee30-tom | 0.050179 | 0.043700 | 0.043720 | 4.4e-04 | 2 | 25 | signed off | 0.0% | 12228 |
| ieee30-u33 | 0.049971 | 0.042919 | 0.042927 | 2.0e-04 | 3 | 30 | signed off | 0.0% | 51044 |
| ieee30-physical | 0.051871 | 0.041474 | 0.041474 | 0.0e+00 | 2 | 48 | signed off | 0.0% | 85584 |

| Case | h5 ratio | h7 ratio | h11 ratio | h13 ratio | Posterior worst | Constraint |
| --- | --- | --- | --- | --- | --- | --- |
| ieee30-tom | 0.834 | 0.895 | 0.867 | 0.905 | -1.12e-3 | IHD, h7, bus 30 |
| ieee30-u33 | 0.906 | 0.966 | 0.785 | 0.859 | -7.51e-3 | IHD, h11, bus 30 |
| ieee30-physical | 0.884 | 0.483 | 0.914 | 0.623 | -1.91e-3 | IHD, h11, bus 7 |

Tom's saved report uses Mosek, an earlier LCM and 10 threads. Its nominal capacity, lower
and upper bounds, gap, resonance, harmonic ratios and posterior metric agree with this run
at the printed precision.

Tom's report records 7 s and 12613 stored pairs. The current report gives its elapsed time
and stored-pair count. Solver tolerances and the LCM version can affect finite-difference
search points and their reuse. The thread count and compilation can affect elapsed time.

## Interpretation

All three cases completed with solver sign-off and no violated posterior points. This is
the numerical result under the stated physical premise and input ranges.

The ACSR checks use catalogue R20 for DC and case values for X and C. The earth-wire diagnostic uses the same low-frequency method as the phase specs.
The scan writes a text companion report for the branch table, data checks and resonances.
The family files also define normalized jacket and diameter coordinates, held at nominal when absent from the box. The physical inputs cite their sources.

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

The factor ranges follow `examples/rhhc/ieee30_robust.jl:58-60`. They describe
network parameters. The pending-audit section of `data/register/README.md` lists the
candidate physical realizations. Hello-world uses the audited U33 and ACSR families,
with the U33 coordinates and earth interval of `ieee30-u33`. Its overhead family
coordinates stay at nominal values. O33-H has no earth wire.

Run a new case with the same environment and solvers:

```sh
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/hello-world.toml
for case in ieee33 industrial antwerp hello-world star3 star5 ieee9 ieee39; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl "cases/$case.toml"
done
```

Only hello-world has LCM circuits among these new cases. Its scan writes the sequence
table and modal resonances. The capacity runner also scans the nominal network for
cases that use source branches. A failed run writes its exception message to
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
is capped by the inner ratio, as EN 50182 requires. The endpoint checks agree
with `derived.json` within 1e-6 m, and the sag calculation preserves `BigFloat`.

Each sweep moves one candidate to each endpoint with all others nominal. It evaluates
every affected circuit at 50 and 650 Hz through the extension's `zy(p, f)` protocol.
The maximum of |Δ|/|nominal| over z, y and all three sequences admits the coordinate
when it reaches 0.001. This rule applies equally to the previously recommended inputs.
The table records the endpoint and result that give that maximum. The complete report
is `output/ieee30-full-sweep.txt`, and the register records each decision.

| Coordinate | Maximum change % | Circuit | Sequence | Quantity | Hz | End | Admitted |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Rscreen_U33` | 3.78283 | U33-800/trefoil | zero | z | 50 | lower | yes |
| `theta_ACSR` | 6.93177 | Panther/O132-D | positive | z | 50 | upper | yes |
| `eps_semicon_U33` | 0.554279 | U33-630/trefoil | positive | y | 50 | lower | yes |
| `t_ins_U33` | 4.40394 | U33-800/trefoil | zero | y | 650 | lower | yes |
| `rho_inner_U33` | 0 | U33-500/trefoil | zero | z | 50 | lower | no |
| `rho_EARTH` | 14.5575 | Goat/O132-D | zero | z | 50 | upper | yes |
| `rho_copper_U33` | 0.0355754 | U33-400/trefoil | negative | z | 650 | upper | no |
| `t_sc1_U33` | 0.538328 | U33-400/trefoil | zero | y | 50 | lower | yes |
| `eps_U33` | 12.494 | U33-500/trefoil | negative | y | 650 | lower | yes |
| `eps_jacket_U33` | 4.61172e-14 | U33-800/trefoil | zero | y | 50 | lower | no |
| `lay_steel_ACSR` | 0.00724018 | Panther/O132-D | negative | z | 50 | lower | no |
| `E_O132` | 0 | Lion/O132-D | zero | z | 50 | lower | no |
| `ambient_O132` | 3.15003 | Panther/O132-D | positive | z | 50 | upper | yes |
| `R20_U33` | 4.38329 | U33-400/trefoil | negative | z | 50 | lower | yes |
| `tand_xlpe_U33` | 0 | U33-500/trefoil | zero | z | 50 | lower | no |
| `rho_aluminium_ACSR` | 0.390942 | Panther/O132-D | positive | z | 50 | lower | yes |
| `lay_outer_ACSR` | 0.254882 | Panther/O132-D | positive | z | 50 | lower | yes |
| `rho_outer_U33` | 0 | U33-500/trefoil | zero | z | 50 | lower | no |
| `d_wire_U33` | 0.0241173 | U33-630/trefoil | positive | z | 50 | lower | no |
| `lay_inner_ACSR` | 0.483406 | Panther/O132-D | positive | z | 50 | lower | yes |
| `expansion_ACSR` | 0 | Lion/O132-D | zero | z | 50 | lower | no |
| `depth_U33` | 0.000618794 | U33-800/trefoil | zero | z | 50 | upper | no |
| `theta_U33` | 25.9813 | U33-400/trefoil | zero | z | 50 | upper | yes |
| `d_steel_ACSR` | 0.128726 | Goat/O132-D | positive | z | 50 | lower | yes |
| `diameter_O132` | 0.150977 | Goat/O132-D | zero | z | 50 | upper | yes |
| `E_ACSR` | 0 | Lion/O132-D | zero | z | 50 | lower | no |
| `rho_steel_ACSR` | 0.120517 | Panther/O132-D | negative | z | 50 | upper | yes |
| `rho_xlpe_U33` | 0 | U33-500/trefoil | zero | z | 50 | lower | no |
| `mu_steel_ACSR` | 0.129885 | Goat/O132-D | positive | z | 50 | upper | yes |
| `rho_steel_O132` | 0.530452 | Goat/O132-D | zero | z | 50 | upper | yes |
| `shield_sag_O132` | 0.0298361 | Goat/O132-D | zero | z | 650 | upper | no |
| `rho_jacket_U33` | 0 | U33-500/trefoil | zero | z | 50 | lower | no |
| `d_c_U33` | 1.27447 | U33-400/trefoil | zero | y | 50 | upper | yes |
| `t_sc2_U33` | 0.458577 | U33-500/trefoil | zero | z | 650 | lower | yes |
| `mu_steel_O132` | 3.11603 | Goat/O132-D | zero | z | 650 | upper | yes |
| `t_jkt_U33` | 2.0171 | U33-800/trefoil | positive | z | 50 | upper | yes |
| `lay_O132` | 0.0316726 | Goat/O132-D | zero | z | 50 | lower | no |
| `phase_sag_O132` | 0.0296214 | Goat/O132-D | zero | z | 650 | upper | no |
| `d_aluminium_ACSR` | 0.619491 | Panther/O132-D | negative | z | 50 | lower | yes |

Loss tangent and dielectric DC resistivities have zero effect under the selected lossless
laws. Jacket permittivity has no measurable effect with the screen grounded. The copper
resistivity, screen wire diameter, burial depth, steel lay ratio, earthwire lay ratio and
reference sags have changes smaller than the threshold. They remain at nominal. Modulus and expansion
also remain nominal because the sweep changes them at the reference temperature, where
the thermal strain is zero. Their effects at other temperatures are interactions that
this prescribed sweep does not test. Fixed support dimensions, earth permittivity and
other fixed register inputs remain nominal for their documented reasons.

The old cases' nominal sequence results are bit-identical to HEAD at 50 and 650 Hz.
The sag endpoint checks agree with the register within 1e-6 m.


The Bear corner with both wire diameters at their upper tolerance maps with
`t = 0.825704225352`. Its conductor diameter is 23.6845 mm, 1 percent above nominal,
and its mass is 1.23708333959 kg/m, 1.985436 percent above nominal. The sweep checks
these values through the realized LCM design. All 128 endpoint combinations of wire
diameters and lay ratios satisfy the aggregate limits and lay order.

Removing the separate ACSR resistance coordinate reduces the selected box from 24 to 23 coordinates. The remaining
sweep decisions are unchanged. The inner-lay maximum changes from 0.115048 percent
to 0.483406 percent with the outer-lay cap, and remains above the threshold.

The capacity run uses 16 Julia threads and Clarabel. The nominal stochastic capacity
is 0.051871 pu with OPTIMAL termination. The capacity calculation signs off after five
iterations in 316 s, with lower bound 0.039870 pu, upper bound 0.039882 pu and gap
2.9e-4. The final master status is OPTIMAL. The 1024-point posterior check reports
zero violations and a worst metric of -3.48e-2 for IHD at harmonic 11, bus 7.
The first nominal resonances are h = 6.2720 and h = 12.3169. The coordinate count
is 23, so HPM omits the vertices from its dictionary. The reports are
`output/ieee30-full.txt`, `output/ieee30-full-scan.txt` and
`output/ieee30-full-sweep.txt`.

The normalized coordinates preserve the full case's physical intervals. Matching points
under the previous and current coordinate representations agree within 1.5e-14 relative
in cable sequence impedance and admittance. Floating-point evaluation order changes the
oracle search trajectory: the reported lower and upper bounds remain 0.039870 and
0.039882 pu, while the final gap and stored-pair count change. The earlier U33 boxes now
reach each size's own registered resistance lower bound.

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
