# Uncertainty levels

The generated independent equivalent box has **35 coordinates**, compared with
23 source coordinates. It does not reach HPM's vertex threshold of ten.
Exact lumping and verified sharing produce 40 candidates. The 0.1 percent rule
omits five radial widths. The independent equivalent intervals lose source
correlations and increase the dimension by twelve.

The corner adversary finds positive own-case violations for both reported
solutions: 0.0029374463 on the source box and 0.0092172411 on the equivalent box.
The random checks missed these points. The violations invalidate the feasibility
certificates needed for the reported lower bounds.

The case is [ieee30-full-equivalent.toml](cases/ieee30-full-equivalent.toml),
generated from `ieee30-full.toml`. Its network and branch assignments are unchanged.
The [coordinate data](data/equivalent/ieee30-full.toml) contain the definitions,
physical members, nominal values, bounds, scales and the source records. The
[reduction report](output/ieee30-full-equivalent-reduction.txt) retains the coordinate
table, sharing groups, subdivision cells, containment and five largest added
extremes. The separate comparison script computes the capacity checks. Raw numerical tables and logs stay in
`tmp/`.

The generated data identify source-case revision `85cdb95`. The capacity
comparison uses HPM `2262246` and LCM pinned at
`59cc3f2329e5df9b33929f52009e996cee1d4784` in the detached
`LineCableModels-hpm` worktree. The study and kit environments declare that
source path. Generated data and reports record the pin. Each case has eight
circuit designs, solid grounding and ideal transposition.

The plain-tube reconstruction preserves the coaxial computation before screening.
The admission rule then fixes five varying widths at nominal. Containment
applies to the admitted equivalent members. It does not establish exact containment
of the full source physics after screening. The cross-check below measures the
resulting solutions on both models.

## Effective domain and subdivision

The source maps produce the effective inputs consumed by the designs. `realize` applies
the lay-order map and the EN 50182 radial projection. `effective` resolves phase
temperature with the ambient constraint. The design consumes the resulting
values. The effective set keeps the family ranges, outer lay at most inner lay,
phase temperature at least ambient, wire tolerances, conductor diameter limits
of one percent and mass limits of two percent.

The observation vector contains the 636 coaxial inventory entries and 452
entries for the equivalent parameterization. The analysis evaluates nominal,
46 axis endpoints and 256 Sobol points, with HPM seed 1. Sensitivities use feasible
effective perturbations with a half step of 0.0001 of the coordinate interval.
The normalized sensitivity zero band remains 0.000001.

For a cell with conflicting signs, `split_cell` tests splits on its driving
effective coordinates. Its score counts disagreeing derivative pairs left in
the children. It chooses the smallest score, then bisects a bracketed derivative
zero when the witnesses permit that refinement. The recursion stops at depth
six and sends a still inconsistent cell to the interior search described below.
The algorithm uses measured values and
sensitivities to select its coordinate and threshold.

The phase heights split once at **20.000000178892172 °C**, within
0.00000018 °C of the 20 °C reference in the input records. This result applies to
Bear, Goat, Lion and Panther. Each of their three phase heights has the cells:

| Cell | Effective phase temperature °C | Samples | Depth |
| --- | --- | ---: | ---: |
| Cold | [-5, 20.000000178892172] | 89 | 1 |
| Hot | [20.000000178892172, 75] | 214 | 1 |

The remaining effective constraints restrict each cell. The cells apply to
the twelve phase heights, represented in both the raw and equivalent inventory.
Of the 1088 observations, 24 have two cells, 1040 have one cell, and 24 fixed
infinite earth properties remain constant. The maximum observed depth is one.
The reduction report lists these subdivision cells for every overhead design.

## Hull search and containment

The search closes each quantity's driving set over the effective constraints.
Independent monotone coordinates take their oriented endpoints. It enumerates
the constrained corners and intersects infeasible corner rays with order and
aggregate limits by bisection. It removes coordinatewise dominated candidates.
The search refines the retained candidates on coordinate lines and active faces.
On a face it solves a companion coordinate by bisection and searches the free
coordinate with golden-section steps. Order ties remain active faces, including
the edge where inner and outer lays agree. The face search probes nine values
before its scalar search of 32 golden-section steps.
A final refinement repeats up to four passes.

Affine aliases within a circuit share a numerical search only after their driving
sets and normalized profiles agree within 1e-11. Their bounds include the measured
alias error. Finite search bounds use 64 machine spacings for outward rounding.
Fixed values use a roundoff envelope. The temporary hull and cell tables record the method, bounds and cell
intervals per quantity. The overall interval endpoints are the extreme cell endpoints.

The containment run uses **4096 additional Sobol points per circuit**, with seed
3, mapped through the source maps. A point must lie inside the value interval of
a cell containing its effective inputs. The comparison has no tolerance.
The result is **zero escapes on all eight circuits**. A separate check of the
303 analysis points per circuit also has zero escapes. The largest normalized
outward margin is **-7.389644451904987e-15**. A negative margin means the point is
inside the cell interval. The margin is the larger bound violation divided by
the larger magnitude of the overall interval endpoints.

The bounds are numerical search results with the stated sample verification.
The numerical settings above also apply to the fixed and unsplit quantities.

## Coordinate definitions and exact lumping

Each coordinate is a signed normalized position. Zero gives the nominal value,
minus one the lower hull bound and plus one the upper bound. Positive and negative sides use
their own widths. The material and geometry values use SI units.

The equivalent design has plain tubes and annular layers. Positive radial widths
reconstruct the radii and their common interfaces. A tube's operating DC
resistance replaces its operating resistivity through
`rho = resistance * pi * (outer^2 - inner^2)`, with the reference-temperature
conversion retained. Its relative permeability preserves LCM's equivalent
geometric mean radius. The tube resistances, widths and permeabilities together
retain the source coaxial computation. The source-value round trip passes the
1e-12 sequence tolerance at the 303 points of all eight circuits. The direct blueprint
reconstruction has zero measured difference in both blueprint and sequence checks.

The U33 core DC resistance combines R20 and temperature in one value. The screen
resistance combines screen R20 and its temperature. The stranded conductors'
DC resistance combines the wire conductances and lay effects that LCM flattens.
Their radii and equivalent permeabilities remain separate quantities. The height
combines the complete sag calculation into one placement value per circuit.

Contact trefoils keep LCM's spacing law with their reconstructed outside radius.
The formation center comes from the inverse public trefoil placement. The box
keeps radial adjacency and the formation law. Independent widths can produce
cumulative radii outside the source image, as their source correlations are lost.
Phase heights within an overhead circuit share one signed position and
keep their fixed height offsets.

Sharing across designs requires agreement of normalized positions at **4399
source points**: the 303 analysis points and 4096 verification points. Values within 64 machine spacings of nominal use position zero. The
complete-link tolerance is **1e-10**. The largest error within a shared group is
**1.1357581541915351e-12**. Distinct groups of the same family and quantity retain
their measured differences in `equivalent-sharing.csv`.

The table lists the coordinates and their physical bounds. Height bounds use
common shifts from nominal, in meters. Other bounds use the physical values.
Repeated phases have the same tube and layer values. The generated coordinate data give the nominal and bounds for each physical field.

| Coordinate | Definition and circuits | Bounds | Unit | Maximum change % | Admit |
| --- | --- | --- | --- | ---: | --- |
| `width_ACSR_1` | `cable.1.tube.1.width` on Bear | [0.01160775, 0.01184225] | m | 0.1664268 | yes |
| `resistance_ACSR_1` | `cable.1.tube.1.resistance` on Bear | [9.176931151e-05, 0.0001352692586] | Ω/m | 7.0234559 | yes |
| `mu_r_ACSR_1` | `cable.1.tube.1.mu_r` on Bear | [0.7938111248, 0.8745211544] | 1 | 0.21569331 | yes |
| `width_O132_1` | `cable.4.tube.1.width` on Bear, Goat, Lion, Panther | [0.00465, 0.0048] | m | 0.014512891 | no |
| `resistance_O132_1` | `cable.4.tube.1.resistance` on Bear, Goat, Lion, Panther | [0.002920941552, 0.004488268038] | Ω/m | 1.0117059 | yes |
| `mu_r_O132_1` | `cable.4.tube.1.mu_r` on Bear, Goat, Lion, Panther | [8.286102801, 44.00092835] | 1 | 3.1160342 | yes |
| `y_ACSR_1` | `cable.1.y` on Bear | [-1.260762515, 0.6017071499] | m | 1.1467578 | yes |
| `y_O132_1` | `cable.4.y` on Bear, Goat, Lion, Panther | [-0.4042071722, 0.4054566001] | m | 0.34389857 | yes |
| `rho_EARTH_1` | `earth.rho` on Bear, Goat, Lion, Panther, U33-400, U33-500, U33-630, U33-800 | [10, 1000] | Ω·m | 14.557496 | yes |
| `width_ACSR_2` | `cable.1.tube.1.width` on Goat | [0.01285515, 0.01311485] | m | 0.16928823 | yes |
| `resistance_ACSR_2` | `cable.1.tube.1.resistance` on Goat | [7.482371521e-05, 0.0001102912108] | Ω/m | 5.8738207 | yes |
| `mu_r_ACSR_2` | `cable.1.tube.1.mu_r` on Goat | [0.7934555602, 0.8746559128] | 1 | 0.22170948 | yes |
| `y_ACSR_2` | `cable.1.y` on Goat | [-1.260679911, 0.6016396079] | m | 1.1578442 | yes |
| `width_ACSR_3` | `cable.1.tube.1.width` on Lion | [0.0110187, 0.0112413] | m | 0.16500378 | yes |
| `resistance_ACSR_3` | `cable.1.tube.1.resistance` on Lion | [0.0001018433901, 0.0001501185925] | Ω/m | 7.6802799 | yes |
| `mu_r_ACSR_3` | `cable.1.tube.1.mu_r` on Lion | [0.793583551, 0.8746071825] | 1 | 0.2128802 | yes |
| `y_ACSR_3` | `cable.1.y` on Lion | [-1.260400888, 0.6014115466] | m | 1.1407411 | yes |
| `width_ACSR_4` | `cable.1.tube.1.width` on Panther | [0.010395, 0.010605] | m | 0.1634402 | yes |
| `resistance_ACSR_4` | `cable.1.tube.1.resistance` on Panther | [0.0001144312332, 0.0001686732505] | Ω/m | 8.4736109 | yes |
| `mu_r_ACSR_4` | `cable.1.tube.1.mu_r` on Panther | [0.7933148859, 0.874709761] | 1 | 0.20943265 | yes |
| `y_ACSR_4` | `cable.1.y` on Panther | [-1.260636399, 0.6016040346] | m | 1.1347783 | yes |
| `width_U33_1` | `cable.1.tube.1.width` on U33-400, U33-500, U33-630, U33-800 | U33-400: [0.0115, 0.0117], U33-500: [0.01325, 0.01345], U33-630: [0.01495, 0.01515], U33-800: [0.0171, 0.0173] | m | 1.2744714 | yes |
| `resistance_U33_1` | `cable.1.tube.1.resistance` on U33-400 | [6.871452e-05, 9.974738e-05] | Ω/m | 15.578143 | yes |
| `width_U33_2` | `cable.1.tube.2.width` on U33-400 | [0.0004687207227, 0.0005221717484] | m | 0.065109031 | no |
| `resistance_U33_2` | `cable.1.tube.2.resistance` on U33-400, U33-500, U33-630, U33-800 | [0.0002411357, 0.0003353513] | Ω/m | 22.587672 | yes |
| `width_U33_3` | `cable.1.layer.1.width` on U33-400, U33-500, U33-630, U33-800 | U33-400/U33-500: [0.000765, 0.000935], U33-630/U33-800: [0.000945, 0.001155] | m | 0.53832816 | yes |
| `eps_r_U33_1` | `cable.1.layer.1.eps_r, cable.1.layer.3.eps_r` on U33-400, U33-500, U33-630, U33-800 | [100, 10000] | 1 | 0.55427886 | yes |
| `width_U33_4` | `cable.1.layer.2.width` on U33-400, U33-500, U33-630, U33-800 | [0.0076, 0.0084] | m | 4.4039383 | yes |
| `eps_r_U33_2` | `cable.1.layer.2.eps_r` on U33-400, U33-500, U33-630, U33-800 | [2.1, 2.5] | 1 | 12.494039 | yes |
| `width_U33_5` | `cable.1.layer.3.width` on U33-400, U33-500, U33-630, U33-800 | [0.00081, 0.00099] | m | 0.46156282 | yes |
| `width_U33_6` | `cable.1.layer.4.width` on U33-400 | [0.003677828252, 0.004731279277] | m | 1.6908036 | yes |
| `resistance_U33_3` | `cable.1.tube.1.resistance` on U33-500 | [5.345529e-05, 7.756705e-05] | Ω/m | 13.385707 | yes |
| `width_U33_7` | `cable.1.tube.2.width` on U33-500 | [0.000434863113, 0.00048231521] | m | 0.059207296 | no |
| `width_U33_8` | `cable.1.layer.4.width` on U33-500 | [0.00381768479, 0.004895136887] | m | 1.8036029 | yes |
| `resistance_U33_4` | `cable.1.tube.1.resistance` on U33-630 | [4.136307e-05, 6.013049e-05] | Ω/m | 11.156839 | yes |
| `width_U33_9` | `cable.1.tube.2.width` on U33-630 | [0.000402890182, 0.0004457068612] | m | 0.054479586 | no |
| `width_U33_10` | `cable.1.layer.4.width` on U33-630 | [0.003954293139, 0.005067109818] | m | 1.8955202 | yes |
| `resistance_U33_5` | `cable.1.tube.1.resistance` on U33-800 | [3.243786e-05, 4.705307e-05] | Ω/m | 9.2310187 | yes |
| `width_U33_11` | `cable.1.tube.2.width` on U33-800 | [0.0003722087909, 0.0004100607692] | m | 0.04853768 | no |
| `width_U33_12` | `cable.1.layer.4.width` on U33-800 | [0.004239939231, 0.005437791209] | m | 2.0354937 | yes |

The maximum change uses the `sweep.jl` measure. It compares each hull endpoint
with nominal, with the other equivalent coordinates nominal, over all applicable
circuits, zero, positive and negative sequences, z and y, at 50 and 650 Hz.
HPM's `ExactPiModel` obtains those sequence results through the extension.
The admission threshold is 0.1 percent. The omitted widths are the earthwire
radius and the four U33 screen tube widths. Their reported changes range from
0.0145 to 0.0652 percent. The omitted widths keep their nominal values in the generated parameterization.
The admission test measures each coordinate independently.
The remaining 35 coordinates comprise 16 ACSR,
15 U33, three earthwire and structure and one earth coordinate.

## Source correlations lost by boxing

The dependence analysis lists the admitted equivalent coordinates driven by
each source coordinate. Examples with several outputs include:

- `theta_ACSR`: four tube DC resistances, four equivalent permeabilities and
  four phase-height coordinates.
- `ambient_O132`: those twelve phase coordinates, earthwire DC resistance and
  earthwire height, through the temperature-order map and the earthwire state.
- `theta_U33`: four core DC resistances and the shared screen DC resistance.
- `R20_U33`: four core DC resistances, which fail the normalized sharing test.
- `Rscreen_U33`: screen DC resistance and four jacket widths.
- The ACSR wire diameters and lays: radius, resistance, permeability and height
  coordinates of several designs, including the source projection's coupling.

The Bear box admits a tube DC resistance of 0.00013526925858673306 Ω/m
together with the highest, cold phase placement of 14.601707149929336 m. The resistance maximum needs the hot phase
state, whereas the highest placement needs the cold phase state. The independent
box drops this relation. The tests construct that equivalent combination and
check both quantities through LCM's public blueprint path.

The analysis finds 181 admitted coordinate pairs with common source inputs. The smallest reported joint residual is
0.0033113770650823326. The Bear resistance and height example has a residual of
0.9153060659381409. The reduction report gives the five largest added extremes, their shared inputs,
boxed endpoints and physical values. The joint
residual is the maximum discrepancy of the two signed normalized positions.
The source search starts at the closest of the 4399 sampled points and refines
on the driving source coordinates. If it reaches the first selected corner,
it tries another corner. The residual records the smallest value found by that
search. The raw pair table retains the closest source points in `tmp/`.

## Capacity comparison

All four cases use 16 Julia threads and Clarabel with the unchanged runner
settings. HPM returns OPTIMAL and signs off all runs. The independent corner
checks refute the own-case feasibility reported for the source and equivalent
solutions. The surrogate solutions pass the sampled own-case checks.

| Quantity | Source | Equivalent | Factor | Scaled |
| --- | ---: | ---: | ---: | ---: |
| Coordinates | 23 | 35 | 6 | 6 |
| Nominal stochastic capacity, pu | 0.0518710788968 | 0.0518710788967 | 0.0518710788968 | 0.0395783138504 |
| Reported lower bound, pu | 0.0398702070185 | 0.0398798701882 | 0.0383378103716 | 0.0259375640772 |
| Reported upper bound, pu | 0.0398819424728 | 0.0398896113566 | 0.0383378103716 | 0.0259375640772 |
| Relative gap | 0.000294254834209 | 0.000244203141187 | 0 | 0 |
| Iterations | 5 | 2 | 2 | 2 |
| Stored LCM point-frequency pairs | 659896 | 311792 | 7720 | 8 |
| Robust solve time, s | 282.779629618 | 155.068624618 | 9.24357309 | 9.503212659 |

The runner labels the stored-pair count as LCM solves. It sums the circuit table
sizes. Solve time covers the capacity calculation with uncertainty. These wall
times include concurrent studies checks. The reported source and equivalent
intervals overlap, and their measured difference does not establish an increase
in capacity from boxing.

The ratios below compare each level's harmonic current with its own nominal
current. The factor nominal capacity equals the source nominal capacity. The
scaled nominal capacity changes because its frequency law changes the networks.

| Harmonic | Source robust HHC / nominal | Equivalent robust HHC / nominal | Factor robust HHC / nominal | Scaled robust HHC / nominal |
| --- | ---: | ---: | ---: | ---: |
| 5 | 0.827083765 | 0.825754342 | 0.800646012 | 0.769870546 |
| 7 | 0.448395565 | 0.445278044 | 0.421842240 | 0.487712312 |
| 11 | 0.903713307 | 0.906295442 | 0.875112569 | 0.592787197 |
| 13 | 0.579774842 | 0.578463642 | 0.541944938 | 0.699804195 |

Each `posterior_check` uses 1024 random points, seed 2 and tolerance 0.0001.
Negative maxima satisfy the checked constraints. The scaled solution violates
source limits at 124 sampled points, or 12.109375 percent of the sample.

| Solution | Model and box checked | Maximum normalized violation | Constraint | Violated points |
| --- | --- | ---: | --- | ---: |
| Source | Source | -0.0347936668576 | IHD, harmonic 11, bus 7 | 0 |
| Equivalent | Equivalent | -0.0210248161584 | IHD, harmonic 11, bus 7 | 0 |
| Equivalent | Source | -0.0292701170284 | IHD, harmonic 11, bus 7 | 0 |
| Factor | Factor | -0.0368241455083 | IHD, harmonic 5, bus 30 | 0 |
| Factor | Source | -0.0949206614722 | IHD, harmonic 11, bus 7 | 0 |
| Scaled | Scaled | -0.0475348757066 | IHD, harmonic 5, bus 30 | 0 |
| Scaled | Source | 0.406182843924 | IHD, harmonic 13, bus 30 | 124 |

The factor solution passes both source transfer checks. This measured result
establishes feasibility only at the checked source points. The family-factor
set remains neither an inner nor an outer approximation of the source set.
The scaled solution fails transfer despite its lower total capacity. Its
harmonic-13 current is 0.0003569555 pu, compared with 0.0002682928 pu for the
source solution. The different harmonic allocation violates source IHD limits.

The [scan report](output/ieee30-full-equivalent-scan.txt) and
[sequence table](output/ieee30-full-equivalent-sequences.csv) retain the electrical
checks. The source, equivalent and factor nominal first resonances are harmonic
orders 6.27 and 12.32. Each capacity report contains only the output of `run.jl`.
The [comparison report](output/ieee30-full-comparison.txt) contains all levels,
decision vectors, posterior checks and the corner adversary.
The script `compare.jl` writes that report directly.

## Corner adversary

HPM omits dictionary vertices for the source and equivalent boxes because their
dimensions exceed ten. It keeps the 64 vertices of each six-coordinate surrogate
box. The nominal-sign adversary evaluates the distinct predicted points below.
Positive maxima violate the indicated constraints.

| Solution | Case | Corners | Maximum normalized violation | Constraint |
| --- | --- | ---: | ---: | --- |
| Source | Source | 88 | 0.0029374462668 | IHD, harmonic 5, bus 30 |
| Equivalent | Equivalent | 134 | 0.00921724112709 | IHD, harmonic 11, bus 7 |
| Equivalent | Source | 90 | 0.00753486516911 | IHD, harmonic 11, bus 7 |
| Factor | Factor | 28 | -8.79063488668e-12 | IHD, harmonic 11, bus 30 |
| Factor | Source | 87 | -0.0601556299595 | IHD, harmonic 5, bus 30 |
| Scaled | Scaled | 25 | -4.46839232282e-11 | IHD, harmonic 11, bus 30 |
| Scaled | Source | 82 | 0.769040553209 | IHD, harmonic 13, bus 30 |

The source and equivalent own-case certificates fail this check. The equivalent
and scaled transfers to the source also fail. Both surrogate own-case checks
have negative maxima. The factor solution passes the source corner check.
A clean sampled check does not establish feasibility throughout a box. The
comparison report gives every coordinate of each worst point, including signed
positions for registered source intervals.

For each solution and target case, the adversary evaluates every positive-limit
constraint at the case's nominal point. It uses `moments_at`, `shhc_constraints`
and `normalized_violation` with squared currents `g_h = I_h²`. Finite differences
respect each coordinate's linear or logarithmic scale, use one-sided differences
at an endpoint and classify normalized derivatives as negative, zero or positive with `METHOD.zero_band`.
The derivative scale corresponds to a signed latent coordinate on [-1, 1].

Each constraint predicts an adverse point: upper bound for a positive sign,
lower bound for a negative sign and nominal for zero sensitivity. The adversary
removes duplicate points and evaluates all remaining points with `posterior_check`
on a zero-width `ParameterBox`, with one sample per point. The source builders
apply the same physical maps as at every other source point. The report gives
the count, maximum normalized violation, constraint and physical point for each
of the own-case and source-transfer combinations.

A positive own-case violation refutes the reported feasibility of that solution
on its own box. A positive transfer violation refutes feasibility of the level's
solution on the source box, without by itself refuting its own-case certificate.
The adversary is not a proof: resonance constraints can change sensitivity sign
and have larger violations away from the predicted corners.

## Construction and reproduction

The generated case uses the existing `problem` → `arrangement` → `design` →
`system` chain. Its loader reads normalized coordinate declarations into the
families. The existing bindings, `coordinates` and `interpolate` resolve each
physical member. The equivalent `design` method builds plain tubes and annular
layers through the public LCM API. The overhead arrangement supplies the saved
equivalent mean heights with zero remaining sag. The existing `system` places
and connects the conductors. All runners consume the same construction chain.

From the repository root:

```sh
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/reduce.jl cases/ieee30-full.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/test/runtests.jl
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/ieee30-full-equivalent.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl cases/ieee30-full-equivalent.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/compare.jl cases/ieee30-full.toml cases/ieee30-full-equivalent.toml cases/ieee30-full-factor.toml cases/ieee30-full-scaled.toml
```

Reduction writes the case, coordinate data and compact report. Its temporary
evidence directory defaults to `tmp/equivalent`. A second argument
selects another directory. It recomputes subdivision, hulls, sharing, admission and
decorrelation. Sorted TOML output and a fixed table order make regeneration
deterministic. The source commit and LCM source fingerprint identify the input code.

The reduction report depends only on the generated reduction evidence.
`compare.jl` calls the runner's solve routine for each case with its nominal,
resonance and capacity settings. It captures the prepared network data through
that routine, then calls public HPM functions for the posterior and corner
checks. The runner's standalone report contains no comparison section.

`METHOD` in `builders.jl` defines the method parameters once. Reduction,
admission, subdivision, sharing, verification, the corner adversary and their tests
read those parameters. The finite-difference half step is 0.0001, the sensitivity
zero band is 0.000001, the maximum subdivision depth is six and the sharing
tolerance is 1e-10. Hull verification uses 4096 points and admission requires a
relative endpoint effect of 0.001. These method parameters leave the runner's
solver and posterior defaults unchanged.

The studies tests call the reconstruction, effective-map, subdivision, hull and
generation code. They check equality at 303 source points on eight circuits at
50 and 650 Hz, containment at 4096 mapped verification points per circuit, sharing,
and a boxed resistance and height combination. They regenerate all three files
into `tmp/` and compare their bytes. They also load the generated case, compare
its nominal sequences with the source, check containment against its saved members
and exercise its coordinate interpolation away from nominal.

## Surrogate factors on resistance, inductance and capacitance

The depth cap remains six. A cell that still has conflicting sensitivity signs
uses an interior search. Every driving-coordinate corner and the center starts
both a minimum and a maximum search. Infeasible starts project toward a feasible
sample within the cell. Bounded coordinate sweeps use golden-section steps and
retain the best values found. Monotone cells retain the earlier search unchanged.

`METHOD.interior_starts` selects the corners and center. The relative search
tolerance is 1e-10 and the iteration limit is 32 for coordinate sweeps and scalar
golden-section searches. The cells retain the label `interior search` in the
reduction report, with their depths and starts per extreme. Each affected
quantity also receives 4096 additional Sobol samples distributed among its
inconsistent cells. Sampling projects onto the cell's effective feasible set
when necessary. The report gives each cell's sample count and largest normalized
outward margin. These numerical searches and containment checks do not prove a
global extremum.

Both-end-bonded screens give a physical reason for an interior resistance
maximum. The circulating-current contribution is about
`R_s * X_m^2 / (R_s^2 + X_m^2)`. Its derivative changes sign at `R_s = X_m`.
The operating screen resistance depends on reference resistance and temperature,
while mutual reactance depends on geometry and frequency. The peak condition
defines a curved surface in the effective coordinates.

At any finite subdivision depth, a cell can retain points on both sides of
`R_s = X_m`. The numerical screen-loss evidence below compares an LCM resistance
peak with this estimate.

The generated cases [ieee30-full-factor](cases/ieee30-full-factor.toml) and
[ieee30-full-scaled](cases/ieee30-full-scaled.toml) use six
coordinates: `kR`, `kL` and `kC` for each of U33 and ACSR. HPM keeps the vertices
for these six-dimensional boxes. Both cases use the same generated
data file [data/surrogate/ieee30-full.toml](data/surrogate/ieee30-full.toml),
family classes and factor intervals.
The source case defines the physical designs and the fundamental branch data.

| Coordinate | Nominal | Lower | Upper |
| --- | ---: | ---: | ---: |
| `kR_ACSR` | 1 | 0.8623096825708542 | 1.2682189163872188 |
| `kL_ACSR` | 1 | 0.9962618928782407 | 1.005315726528468 |
| `kC_ACSR` | 1 | 0.9981387764988817 | 1.0018618770161116 |
| `kR_U33` | 1 | 0.9031342364062734 | 1.2323618743221811 |
| `kL_U33` | 1 | 0.9205434341153378 | 1.1543234302904137 |
| `kC_U33` | 1 | 0.8324072347466105 | 1.10520498539662 |

Each interval combines twenty design-frequency ratios.
Panther at 50 Hz sets both endpoints of `kR_ACSR`, and U33-400 at 50 Hz sets
both endpoints of `kR_U33`. Interior maxima at other frequencies still enter the
individual design-frequency hulls.

The source verification
has zero escapes at 4096 mapped points on each of eight designs. Its largest
normalized outward margin is -2.549795008242146e-6. The interior search handles
Bear R at 350 Hz, Goat R at 250 Hz and U33-400 R at 250 Hz. Other quantities use
the earlier cell hull search or its verified affine alias. The inconsistent
cells have depth six.

| Design | Quantity | Cell | Starts per extreme | Extra samples | Escapes | Largest margin |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Bear | R at 350 Hz | 6 | 8193 | 2048 | 0 | -0.016144699608372805 |
| Bear | R at 350 Hz | 7 | 8193 | 2048 | 0 | -0.016157598092505998 |
| Goat | R at 250 Hz | 7 | 8193 | 4096 | 0 | -0.008530009161791089 |
| U33-400 | R at 250 Hz | 5 | 513 | 2048 | 0 | -0.0033373622757674084 |
| U33-400 | R at 250 Hz | 9 | 513 | 2048 | 0 | -0.003105966575570989 |

The additional checks total 4096 feasible cell samples per affected quantity,
with zero escapes. Their largest normalized outward margin is
-0.003105966575570989. The unchanged equivalent case, coordinate data and
reduction report regenerate byte-identically under the LCM pin.

At 250 Hz, LCM gives an interior R-ratio maximum for U33-400 as screen resistance
varies with the other effective inputs fixed. The operating screen-resistance
range is [0.0002411357, 0.0003353513] Ω/m. The peak occurs at
0.0002735974623 Ω/m, compared with the geometric mutual-reactance estimate
0.0002687993415 Ω/m, a 1.785 percent difference. The estimate uses
`μ₀ * f * log(spacing / screen_outer_radius)`. The R ratios at the two reference
screen-resistance bounds are 1.0125457752 and 1.0126164850. The interior peak is
1.0127378077. These observations support the circulating-current explanation.
The approximate loss term does not replace LCM in the hull computation.

For each source design and frequency, the positive-sequence per-meter quantities
are `R = real(z)`, `L = imag(z)/(2πf)` and `C = imag(y)/(2πf)`. A ratio divides
each quantity by its nominal value at the same frequency. Each family interval
is the hull over all its designs and the case frequencies: 50, 250, 350, 550 and
650 Hz. The extension's callable evaluates every physical point through LCM.

The ratio observations use the same effective axes, feasible-set constraints,
sensitivity parameters, subdivision and scalar hull search as the equivalent
reduction. The face search remains active on order ties and aggregate limits.
Containment checks 4096 additional Sobol points per design, mapped through the
physical source maps. The generated reduction report records the cells, search
methods and containment results. The raw hull table stays in `tmp/`.

The factor case uses the existing `exact` model kind. HPM's public
`ExactPiModel(branches; f1)` receives callables that scale the extension's nominal
per-meter data at each requested frequency:

```text
z = kR * real(z_nominal) + im * kL * imag(z_nominal)
y = real(y_nominal) + im * kC * imag(y_nominal)
```

The real shunt conductance stays nominal. At factors of one, the frequency
behavior is exact. This variant approximates only the uncertainty by the family
factors. Its LCM calls always use the source nominal point. The extension owns
the table of those nominal results.

The scaled case uses the existing `scaled` kind and HPM's `ScaledBranchModel`.
Its base is the LCM nominal exact π data at 50 Hz, as the source builder writes
them into the fundamental network. HPM then scales resistance by `sqrt(h)` and
reactance and shunt susceptance by `h`, with the same family factors. This variant
adds frequency approximation to the uncertainty approximation. The builders need
no additional model kind.

Harmonics 5 and 11 have negative sequence, and 7 and 13 have positive sequence.
With ideal transposition, those two sequence impedances and admittances agree.
The fitted hulls use positive sequence because the harmonic set omits zero sequence. If requested, the factor callable applies the same factors to LCM's nominal
zero sequence, without a zero-sequence containment claim. `ScaledBranchModel`
rejects zero sequence because its base contains only positive-sequence data.

One factor moves every design of a family together, at every frequency. The
surrogate sets are neither inner nor outer approximations of the source set.
A family ratio hull does not preserve the source dependence between designs,
frequencies or R, L and C. Transfers to the source box need their own checks.

`reduce.jl` generates both cases, their shared data and the
surrogate [reduction report](output/ieee30-full-surrogate-reduction.txt).
The source commit in generated headers identifies the last revision of the
source case file. A studies commit that leaves the source case unchanged does not change these headers.
`compare.jl <source> <level>...` runs each level with the runner's settings and
writes each standalone report through `runstudy`. The comparison contains each
level's own random and corner checks and its transfer to the source case.

The factor variant gives a reported capacity of 0.03833781037 pu and passes the
sampled source transfer checks. The scaled variant gives 0.02593756408 pu but
fails source transfer at harmonic 13, bus 30, in both checks. Its lower total
capacity does not make its harmonic allocation feasible on the source box.

The surrogate tests compare factor sequences exactly with source values on
every circuit at 50 and 650 Hz when all factors equal one. Scaled fundamental
branch data must also agree exactly. Tests exercise factors away from nominal
and check source and interior-cell containment. Regeneration tests compare the
bytes of both cases, the shared coordinate data and the reduction report.

## Quantities and consumers

Paths in this table are relative to LCM's `src/engine/`, except for `flatten` in
`src/datamodel/blueprint.jl`. The analysis calls public names and reads the public
records. It does not call the internal input constructor or formula modules.

| Quantity | Unit | Code that reads it |
| --- | --- | --- |
| Tube inner and outer radius | m | `input.jl`, `assembly.jl:impedance!`, `internalimpedance/schelkunoff1934.jl` |
| Tube resistivity at the problem temperature | Ω·m | `compute.jl:solve!`, `resistivity/linear.jl`, `internalimpedance/schelkunoff1934.jl` |
| Tube relative permeability | 1 | `assembly.jl:impedance!`, `internalimpedance/schelkunoff1934.jl` |
| Interval relative permeability with lay correction | 1 | `input.jl`, `assembly.jl:impedance!`, `insulationimpedance/ametani1980.jl` |
| Layer inner and outer radius | m | `input.jl`, `shunt/series.jl:potential!` |
| Layer relative permittivity | 1 | `compute.jl:admittivity`, `insulationadmittance/lossless.jl`, `semiconadmittance/lossless.jl` |
| Layer resistivity and loss tangent | Ω·m and 1 | Kept by `input.jl`. The selected lossless laws do not use them |
| Layer relative permeability | 1 | `flatten` combines it into the interval permeability |
| Assembly centers and cable placements | m | `input.jl` composes them. `assembly.jl` and `earthreturn/martins2026.jl` read the centers |
| Cable outside radius | m | `input.jl`, `earthreturn/martins2026.jl` |
| Earth resistivity, relative permittivity and permeability | Ω·m, 1 and 1 | `earthproperties/constant.jl`, `earthproperties/bottommost.jl`, `compute.jl:to_buffer!`, `earthreturn/martins2026.jl` |
| Earth layer thickness | m | `input.jl` locates the cable's earth layer. The case uses homogeneous earth |

For this formulation, source temperature and reference resistivity enter a tube's
surface impedance through its operating resistivity. There is no separate
post-blueprint lay input: `flatten` puts its effect into the tube material and
interval permeability. Geometry remains an independent input to the surface
impedance. A resistance coordinate alone cannot replace the tube radii.

