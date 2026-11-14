# Network studies with LineCableModels

These studies join physical cable and overhead-line inputs to HPM's harmonic hosting
capacity models. The [Physical-input studies manual](../../docs/src/linecablemodels_studies.md)
describes the folder layout, input register, case format, construction method and scripts.
The [LineCableModels manual](../../docs/src/linecablemodels.md) gives a runnable extension session.

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

## Checks and results

Check each U33 nominal overall diameter against the catalogue within 0.05 mm. Compare its
low-frequency resistance with R20, its resistance at 90 °C with Olex's AC value, and its
50 Hz reactance and capacitance with Olex. The hot calculation heats both core and screen.
The low-frequency diagnostic uses 0.01 Hz at 20 °C. It approaches DC and includes the
circuit's bonding and earth-return model.

The ACSR checks use catalogue R20 for DC and case values for 50 Hz reactance
and capacitance.
For the earth wire, compare with 3.55 Ω/km and the register interval 3.30–4.01 Ω/km.
The earth-wire diagnostic uses the same low-frequency method as the phase specs,
with three copies as phases at the same formation heights.

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

At ε = 0.05, the Gaussian factor is about 1.644854. The former Cantelli factor
was about 4.358899. With worst-case alignment the variance is zero, so the factor
changes only the label.

The reports use 16 Julia threads, one BLAS thread, Clarabel and 1024 posterior points.
The nominal stochastic capacity uses the nominal point. The robust
hosting capacity varies the declared box.

| Case | Nominal pu | Robust hosting capacity pu | Capacity ratio | Posterior violated | Stored pairs |
| --- | --- | --- | --- | --- | --- |
| ieee30-full | 0.055640 | 0.042112 | 0.757 | 0.0% | 173820 |
| ieee30-physical | 0.055640 | 0.043813 | 0.787 | 0.0% | 73460 |
| ieee30-tom | 0.052095 | 0.045243 | 0.868 | 0.0% | 11778 |
| ieee30-u33 | 0.051815 | 0.044278 | 0.855 | 0.0% | 47532 |
| hello-world | 0.021128 | 0.017518 | 0.829 | 0.0% | 26074 |
| ieee30-physical-length | 0.055640 | 0.042947 | 0.772 | 0.0% | 75840 |

The installed-length case has the same nominal capacity as `ieee30-physical`.
Its robust hosting capacity is 0.000866 pu lower under the adopted length interval.
Both runs signed off with no violated posterior points.

Tom Van Acker's review of PR #35 gives 0.055640 pu nominal and 0.042112 pu robust
hosting capacity for the Gaussian run, a ratio of 0.757. `ieee30-full` reproduces
both capacities.

## Other test systems

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

## Full-uncertainty sweep results

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

The capacity results appear above and in `output/ieee30-full.txt`.
The scan and sweep reports are `output/ieee30-full-scan.txt` and
`output/ieee30-full-sweep.txt`. Their sequence parameters and admission decisions
do not depend on the chance-constraint factor.

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

## Cases and output inventory

The table lists the case files and the outputs of the reproduction commands.
A comparison file is listed with its physical case. It also evaluates the source case.

| Case file | Generated outputs |
| --- | --- |
| `cases/antwerp-physical.toml` | `output/antwerp-physical-comparison.txt`, `output/antwerp-physical-scan.txt`, `output/antwerp-physical-sequences.csv`, `output/antwerp-physical-sweep.txt`, `output/antwerp-physical.txt` |
| `cases/antwerp.toml` | `output/antwerp.txt` |
| `cases/hello-world.toml` | `output/hello-world-scan.txt`, `output/hello-world-sequences.csv`, `output/hello-world.txt` |
| `cases/ieee30-full.toml` | `output/ieee30-full-scan.txt`, `output/ieee30-full-sequences.csv`, `output/ieee30-full-sweep.txt`, `output/ieee30-full.txt` |
| `cases/ieee30-physical-length.toml` | `output/ieee30-physical-length-scan.txt`, `output/ieee30-physical-length-sequences.csv`, `output/ieee30-physical-length.txt` |
| `cases/ieee30-physical.toml` | `output/ieee30-physical-scan.txt`, `output/ieee30-physical-sequences.csv`, `output/ieee30-physical.txt` |
| `cases/ieee30-tom.toml` | `output/ieee30-tom-scan.txt`, `output/ieee30-tom-sequences.csv`, `output/ieee30-tom.txt` |
| `cases/ieee30-u33.toml` | `output/ieee30-u33-scan.txt`, `output/ieee30-u33-sequences.csv`, `output/ieee30-u33.txt` |
| `cases/ieee30-zs-dyn.toml` | `output/ieee30-zs-dyn.txt` |
| `cases/ieee30-zs-unearthed.toml` | `output/ieee30-zs-unearthed.txt` |
| `cases/ieee30-zs-ynd.toml` | `output/ieee30-zs-ynd.txt` |
| `cases/ieee30-zs-ynyn-resistance.toml` | `output/ieee30-zs-ynyn-resistance.txt` |
| `cases/ieee30-zs-ynyn-solid.toml` | `output/ieee30-zs-ynyn-solid.txt` |
| `cases/ieee33-physical.toml` | `output/ieee33-physical-comparison.txt`, `output/ieee33-physical-scan.txt`, `output/ieee33-physical-sequences.csv`, `output/ieee33-physical-sweep.txt`, `output/ieee33-physical.txt` |
| `cases/ieee33.toml` | `output/ieee33.txt` |
| `cases/ieee39-physical.toml` | `output/ieee39-physical-comparison.txt`, `output/ieee39-physical-scan.txt`, `output/ieee39-physical-sequences.csv`, `output/ieee39-physical-sweep.txt`, `output/ieee39-physical.txt` |
| `cases/ieee39.toml` | `output/ieee39.txt` |
| `cases/ieee9-physical.toml` | `output/ieee9-physical-comparison.txt`, `output/ieee9-physical-scan.txt`, `output/ieee9-physical-sequences.csv`, `output/ieee9-physical-sweep.txt`, `output/ieee9-physical.txt` |
| `cases/ieee9.toml` | `output/ieee9.txt` |
| `cases/industrial-physical.toml` | `output/industrial-physical-comparison.txt`, `output/industrial-physical-scan.txt`, `output/industrial-physical-sequences.csv`, `output/industrial-physical-sweep.txt`, `output/industrial-physical.txt` |
| `cases/industrial.toml` | `output/industrial.txt` |
| `cases/star3-physical.toml` | `output/star3-physical-comparison.txt`, `output/star3-physical-scan.txt`, `output/star3-physical-sequences.csv`, `output/star3-physical-sweep.txt`, `output/star3-physical.txt` |
| `cases/star3.toml` | `output/star3.txt` |
| `cases/star5-physical.toml` | `output/star5-physical-comparison.txt`, `output/star5-physical-scan.txt`, `output/star5-physical-sequences.csv`, `output/star5-physical-sweep.txt`, `output/star5-physical.txt` |
| `cases/star5.toml` | `output/star5.txt` |
| `cases/antwerp-mix.toml` | `output/antwerp-mix.txt` |
| `cases/antwerp-physical-mix.toml` | `output/antwerp-physical-mix-comparison.txt`, `output/antwerp-physical-mix.txt` |
| `cases/hello-world-mix.toml` | `output/hello-world-mix.txt` |
| `cases/ieee30-full-mix.toml` | `output/ieee30-full-mix.txt` |
| `cases/ieee30-physical-length-mix.toml` | `output/ieee30-physical-length-mix.txt` |
| `cases/ieee30-physical-mix.toml` | `output/ieee30-physical-mix.txt` |
| `cases/ieee30-tom-mix.toml` | `output/ieee30-tom-mix.txt` |
| `cases/ieee30-u33-mix.toml` | `output/ieee30-u33-mix.txt` |
| `cases/ieee33-mix.toml` | `output/ieee33-mix.txt` |
| `cases/ieee33-physical-mix.toml` | `output/ieee33-physical-mix-comparison.txt`, `output/ieee33-physical-mix.txt` |
| `cases/ieee39-mix.toml` | `output/ieee39-mix.txt` |
| `cases/ieee39-physical-mix.toml` | `output/ieee39-physical-mix-comparison.txt`, `output/ieee39-physical-mix.txt` |
| `cases/ieee9-mix.toml` | `output/ieee9-mix.txt` |
| `cases/ieee9-physical-mix.toml` | `output/ieee9-physical-mix-comparison.txt`, `output/ieee9-physical-mix.txt` |
| `cases/industrial-mix.toml` | `output/industrial-mix.txt` |
| `cases/industrial-physical-mix.toml` | `output/industrial-physical-mix-comparison.txt`, `output/industrial-physical-mix.txt` |
| `cases/star3-mix.toml` | `output/star3-mix.txt` |
| `cases/star3-physical-mix.toml` | `output/star3-physical-mix-comparison.txt`, `output/star3-physical-mix.txt` |
| `cases/star5-mix.toml` | `output/star5-mix.txt` |
| `cases/star5-physical-mix.toml` | `output/star5-physical-mix-comparison.txt`, `output/star5-physical-mix.txt` |

## Script inventory

| Script | Role |
| --- | --- |
| `builders.jl` | Decode cases, construct physical inputs and declare circuit coordinates. |
| `run.jl` | Write nominal capacity, robust hosting capacity and posterior checks. |
| `compare.jl` | Write paired capacity reports, transfer checks and corner checks. |
| `scan.jl` | Write sequence tables, nominal resonances and physical checks. |
| `sweep.jl` | Write endpoint sensitivities and coordinate admission decisions. |
| `networks/convert.jl` | Reproduce the saved IEEE9 and IEEE39 network inputs from PowerImpedance. |
| `data/register/derive.py` | Reproduce the register arithmetic in `data/register/derived.json`. |
| `test/runtests.jl` | Run the study suite. |
| `test/populations.jl` | Check converter parameters, variant networks and refused combinations. |
| `test/coordinates.jl` | Check declarations, undeclared-coordinate probes and installed length. |
| `test/stars.jl` | Check Falcon constructions, exact π-models and sequence references. |
| `test/transmission.jl` | Check the physical transmission cases and Cardinal construction. |
| `test/distribution.jl` | Check the physical feeder and industrial cables. |
| `test/infeasible.jl` | Check reports for infeasible nominal operating points. |
| `test/antwerp.jl` | Check Antwerp cables and IEEE30 cable sequence references. |
| `test/hdata_at.jl` | Check the physical network inputs of the other solvers. |
| `test/manual.jl` | Parse the manual TOML examples and build complete cases. |
