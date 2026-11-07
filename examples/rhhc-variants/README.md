# IEEE30 physical-input study

These studies compute harmonic hosting capacity with frequency-dependent cable and
overhead-line parameters from LineCableModels (LCM). The IEEE30 cases use HPM's
`inductive_cable33` scenario, 21 harmonic units, harmonic orders 5, 7, 11 and 13,
Cantelli bounds and a risk level of 0.05. The fundamental network uses the nominal
LCM parameters. HPM forms the exact π-model at each harmonic frequency.

The nominal capacity comes from `solve_rhhc` with a zero-width parameter box.
The robust HHC is its lower bound over the physical-input box. The cutting-set
method reports its upper bound, relative gap, iteration count and sign-off.
Sign-off is the oracle's numerical stopping condition. A separate posterior check
evaluates 1024 random points of the box, as measured evidence at those points.

| Case | Physical inputs added |
| --- | --- |
| `ieee30-tom` | Tom's placeholder 630 mm² aluminum cable in touching trefoil at 1 m depth, with his four uncertain inputs. The 132 kV lines keep the scenario data without charging. |
| `ieee30-u33` | Four Olex 33 kV cable sizes in touching trefoil, group center at 0.8 m depth, with registered intervals on resistance, insulation thickness, permittivity and earth resistivity. |
| `ieee30-physical` | The U33 cables plus Lion, Goat, Bear and Panther ACSR conductors on O132-D towers, with a steel earth wire and line charging. |
| `ieee30-full` | The physical case with 23 coordinates shared by physical family, including wire dimensions, lay ratios, materials and temperatures. |
| `hello-world` | HPM's hello-world network with three harmonic units, U33-630 cables and Goat lines on O33-H. |

The physical replacements are study assumptions. They change the scenario's network,
including its line charging and resonances. They describe credible constructions, which
need not match the historical equipment of the IEEE30 system.

Cable screens are solidly grounded at both terminals. Connection index `0` lets Kron
reduction eliminate the screen, and overhead earth wires also use index `0`. Every
circuit uses `CoaxialFormulation(transposition = true)`, which applies ideal
transposition. Cross-bonding is not modeled.

The cases declare the transformer vector groups per tapped branch. Branches 12, 15 and
36 are YNyn0. The three-winding unit uses branches 11 and 14 as Yy0 and branch 13 as
Yd11, and branch 16 is YNd11. HPM receives `Yy0` or `Yd11` and keeps the tap magnitude.
The case records declare solid 132 kV neutral earthing and resistance earthing at
33 kV, sized for a 1 kA earth fault. HPM applies vector groups and delta blocking in the
harmonic admittance, while winding earthing impedances enter only its IV formulation.
The run report lists these declarations and their per-unit impedances.

## Inputs and their sources

`data/register/inputs.csv` records nominal values, intervals, units, source classes,
citations, coupling rules and coordinate decisions. Its `README.md` defines the classes
and the reference documents. Every value in `data/materials.toml`, `data/specs/`,
`data/families/` and `data/formations/` cites a register row or states its origin.
The case files give branch assignments, network origins, solver inputs and transformer
assumptions. `data/register/derive.py` produces `derived.json` from these registered
inputs.

`builders.jl` reads these records and constructs the LCM problems. U33 core resistivity
follows the conductor diameter and DC resistance jointly, and the screen tube keeps the
registered screen resistance. Family coordinates interpolate each product's own bounds.
ACSR wire dimensions and lay ratios determine resistance, mass and metal area. EN 50182
diameter and mass limits constrain the realized wires. The builder projects wire
diameters toward their nominal values and caps the outer lay ratio by the inner ratio,
so the parameter domain remains a box.

Each metal's resistivity uses its own temperature and electrical temperature
coefficient, and LCM receives the 20 °C reference temperature. Phase temperature is at
least ambient, and the earth wire takes the ambient temperature. Sag follows the
registered elastic thermal state equation. `sweep.jl` moves each candidate coordinate
to both endpoints, with the other coordinates nominal, at 50 and 650 Hz. A relative
change of at least 0.001 in any sequence impedance or admittance selects a coordinate.
This endpoint sweep does not test interactions between coordinates.

## Reproduce the computations

Use Julia 1.12. From the repository root, instantiate the environments:

```sh
julia --project=examples/rhhc-variants -e 'using Pkg; Pkg.instantiate()'
julia --project=test/linecablemodels -e 'using Pkg; Pkg.instantiate()'
```

Both environments load this checkout of HPM through a relative path, and LCM from
`https://github.com/Electa-Git/LineCableModels.jl` at
`59cc3f2329e5df9b33929f52009e996cee1d4784`, on `release/v0.1.0-datamodel`. Clarabel
solves the capacity problems and Ipopt the fundamental power flow. The scripts set BLAS
to one thread. Run with 16 Julia threads:

```sh
for case in ieee30-tom ieee30-u33 ieee30-physical ieee30-full hello-world; do
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl "examples/rhhc-variants/cases/$case.toml"
    julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl "examples/rhhc-variants/cases/$case.toml"
done
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/sweep.jl examples/rhhc-variants/cases/ieee30-full.toml
```

`output/<case>.txt` gives capacities, bounds, gap, iterations, elapsed solve time,
nominal resonances, harmonic ratios, posterior results and the sum of the LCM circuit
table sizes. That sum counts stored point-frequency results, and concurrent misses can
compute the same key more than once. The oracle's label counts cover nominal, vertex,
corner, seed, background and refined points. Above ten coordinates, HPM uses
sign-predicted corners in place of all vertices.

`output/<case>-scan.txt` reports branch assignments, sequence parameters at 50 Hz,
resonances and catalogue comparisons. `output/<case>-sequences.csv` contains the
frequency scans, and `output/ieee30-full-sweep.txt` the endpoint sensitivities and
wire-construction checks. Each report gives the HPM base commit and the LCM revision.
The scripts write every output file directly.

Run the extension kit and the study tests from the repository root:

```sh
julia --project=test/linecablemodels test/linecablemodels/runtests.jl
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/test/runtests.jl
```

For every IEEE30 case, the study tests call `runstudy` and check four results. The
zero-width box converges in one iteration at HPM's 0.001 stopping tolerance. The 50 Hz
per-unit values have the magnitudes of cables and overhead lines. Each harmonic capacity
compares with its nominal reference. The run signs off with zero posterior violations.
The 50 Hz checks distinguish the handoff's placeholder cable from the registered
physical replacements. The scan's catalogue comparisons explain the correction for
screen circulating currents under the declared bonding.
