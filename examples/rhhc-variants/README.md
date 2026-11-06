# Network studies with LineCableModels

The IEEE30 studies use the `inductive_cable33` scenario and 21 harmonic units at the first 21 loads.
They use harmonic orders 5, 7, 11 and 13, Cantelli bounds and a risk level of 0.05. The fundamental
network uses the nominal LCM branch parameters. HPM computes each branch's exact π-model
from the extension's series impedance and shunt admittance per meter at each frequency.

Each cable screen has a solid bond to ground at both terminals, with connection index `0`.
Each overhead earth wire has a bond at every tower and a solid connection to ground, with index `0`.
Every circuit uses `CoaxialFormulation(transposition = true)`. The studies omit cross-bonding.

| Case | 33 kV cables | 132 kV lines | Box |
| --- | --- | --- | --- |
| `ieee30-tom` | Tom's placeholder `TOM-630`, trefoil in contact, lower cables at 1 m depth | Scenario data, zero shunt capacitance | Tom's permittivity, insulation thickness, core radius and earth resistivity |
| `ieee30-u33` | Four Olex U33 sizes, trefoil in contact, group center at 0.8 m depth | Scenario data, zero shunt capacitance | Family resistance factor, permittivity, insulation thickness factor and earth resistivity |
| `ieee30-physical` | As `ieee30-u33` | Lion, Goat, Bear and Panther ACSR on `O132-D`, with one steel earth wire | As `ieee30-u33` |

The cables on branches 20, 31, 33 and 34 use U33-400. Branch 26 uses U33-630 and branch 21
uses U33-800. The other 33 kV cables use U33-500. The case files give every assignment.
Cable lengths follow the scenario's retained routes. The overhead lengths follow
`length_km = x_pu * zbase / OHL_X[132.0]`, with `OHL_X[132.0] = 0.40 Ω/km`.

## Dataset and construction

`data/materials.toml` holds the material properties. Each file under `data/specs/` holds the
primary quantities of one cable or conductor. The family files declare the coordinates,
their bounds and the fields they act on. Formation files give the placement references,
support coordinates and sags. Numeric quantities use SI units and cite a row of the
physical-input register in `data/register/`, or their source. Transformer assumptions cite Tom Van Acker's proposal.

A compacted core has radius `(d_c + dd)/2` and equivalent resistivity
`kR * R20 * π * radius^2`. Its diameter and resistance jointly set this resistivity
(Gustavsen 2005, Section III-A). The copper screen is a tube with area `rho_Cu / R20_screen`
(TB 531, Section 4.2.2.5 and Gustavsen 2016, conclusion 4). The jacket thickness is
`kj * t_jkt + d_wire - t_tube`. These expressions preserve the coupling of the inputs.
Tom's placeholder retains his core resistivity and varies the radius directly.

HPM evaluates the exact π-model's hyperbolic functions in
`branch_data(model::ExactPiModel, ...)`, in `src/robust/network_model.jl`:

```math
Z = Z_c \sinh(\gamma\ell), \qquad Y/2 = \tanh(\gamma\ell/2)/Z_c,
\qquad \gamma = \sqrt{zy}, \qquad Z_c = \sqrt{z/y}.
```

The studies and the extension supply only `z` in Ω/m and `y` in S/m.

An ACSR terminal contains the steel core from `stranded` and two aluminum wire courses.
LCM's `flatten` reduces them together. The family factor multiplies both metals' resistivities.
Their relative conductance is unchanged. These cases fix the factor at its nominal value. The overhead mean height is the support height minus two thirds of the
sag, the mean of the parabolic span profile.

`builders.jl` lifts TOML records into typed named tuples. It creates one problem function per
spec, formation and earth combination. Every branch of that circuit receives the same
function. The extension computes the sequences and stores the shared results.

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

The reports print these declarations. The network retains the source transformer parameters.
HPM's harmonic model reads vector groups from its `xfmr` table, which has no tap, whereas
the IEEE30 transformers are tapped branches. This proposal needs an upstream
change to the harmonic model. The zero-sequence scans describe the circuits per unit length, with the stated
bonding.

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
or add sourced placement data. A new case lists its branches, formation, earth, temperature,
family bindings, ordered box names and transformer declarations. A family coordinate that
the case does not bind stays at its nominal value. The U33 family defines the jacket factor and diameter
offset. These cases fix both at their nominal values.

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
| ieee30-u33 | 0.049971 | 0.042918 | 0.042927 | 2.0e-04 | 3 | 40 | signed off | 0.0% | 50440 |
| ieee30-physical | 0.051871 | 0.041474 | 0.041474 | 0.0e+00 | 2 | 48 | signed off | 0.0% | 85584 |

| Case | h5 ratio | h7 ratio | h11 ratio | h13 ratio | Posterior worst | Constraint |
| --- | --- | --- | --- | --- | --- | --- |
| ieee30-tom | 0.834 | 0.895 | 0.867 | 0.905 | -1.12e-3 | IHD, h7, bus 30 |
| ieee30-u33 | 0.906 | 0.966 | 0.785 | 0.859 | -7.52e-3 | IHD, h11, bus 30 |
| ieee30-physical | 0.884 | 0.483 | 0.914 | 0.623 | -1.91e-3 | IHD, h11, bus 7 |

Tom's saved report uses Mosek, an earlier LCM and 10 threads. Its nominal capacity, lower
and upper bounds, gap, resonance, harmonic ratios and posterior metric agree with this run
at the printed precision.

Tom's report records 7 s and 12613 stored pairs. The current report gives its elapsed time
and stored-pair count. Solver tolerances and the LCM version can affect finite-difference
search points and their reuse. The thread count and compilation can affect elapsed time.

## Interpretation

All three cases completed with solver sign-off and no violated posterior points. This is
the numerical result under the stated physical premise and input ranges. The transformer
declarations remain outside the harmonic network model.

The ACSR checks use catalogue R20 for DC and case values for X and C. The earth-wire diagnostic uses the same low-frequency method as the phase specs.
The scan writes a text companion report for the branch table, data checks and resonances.
The family files also give the fixed jacket factor and diameter offset. The physical inputs cite their sources.

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
| `industrial` | `../shhc/data/industrial_network_shhc.m` | Source cables with R, L and C factors | `kR_cable`, `kL_cable`, `kC_cable` |
| `antwerp` | `../shhc/data/medium_network_shhc.m` | Source cables with R, L and C factors | `kR_cable`, `kL_cable`, `kC_cable` |
| `hello-world` | `../rhhc/data/hello_world.m` | U33-630 on branches 2 and 3, Goat on O33-H on branch 4 | `kR_U33`, `eps_U33`, `kt_U33`, `rho_earth` |
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
