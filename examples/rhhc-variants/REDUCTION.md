# Equivalent uncertainty

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
`/tmp`.

The study uses HPM base `85cdb95` and LCM pinned at
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
six and reports an unresolved quantity. The algorithm uses measured values and
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
search. The raw pair table retains the closest source points in `/tmp`.

## Capacity comparison

Both cases use 16 Julia threads and Clarabel with the unchanged runner settings.
HPM returns OPTIMAL and signs off both runs. The own-case corner violations
refute the reported feasibility of both solutions.
The source column is a fresh run with the same loaded LCM code as the equivalent
run. Its bounds and stored-pair count agree with the earlier source report.

| Quantity | Source | Equivalent |
| --- | ---: | ---: |
| Coordinates | 23 | 35 |
| Nominal stochastic capacity, pu | 0.051871078897 | 0.051871078897 |
| Reported lower bound, pu | 0.039870207018 | 0.039879870188 |
| Reported upper bound, pu | 0.039881942473 | 0.039889611357 |
| Relative gap | 0.000294255 | 0.000244203 |
| Iterations | 5 | 2 |
| Stored LCM point-frequency pairs | 659896 | 311792 |
| Robust solve time, s | 283.896 | 141.786 |

The runner labels the stored-pair count as LCM solves. It sums the circuit table
sizes. Solve time covers the capacity calculation. The reported intervals overlap.
The measured difference does not establish an increase in capacity from boxing.

| Harmonic | Source robust HHC / nominal | Equivalent robust HHC / nominal |
| --- | ---: | ---: |
| 5 | 0.827084 | 0.825754 |
| 7 | 0.448396 | 0.445278 |
| 11 | 0.903713 | 0.906295 |
| 13 | 0.579775 | 0.578464 |

Each `posterior_check` uses 1024 random points, seed 2 and tolerance 0.0001.
Negative maxima satisfy the checked constraints. All four checks report zero
violated points. Their worst constraint is IHD at harmonic 11, bus 7.

| Solution | Model and box checked | Maximum normalized violation |
| --- | --- | ---: |
| Source | Source | -0.034793667 |
| Equivalent | Equivalent | -0.021024816 |
| Equivalent | Source | -0.029270117 |
| Source | Equivalent | -0.026595283 |

The random cross-checks do not detect a sampled violation due to decorrelation. They test the
reported solutions at the sampled points. The independent intervals still admit
the 181 added pairwise extremes found by the source searches.

The [scan report](output/ieee30-full-equivalent-scan.txt) and
[sequence table](output/ieee30-full-equivalent-sequences.csv) retain the electrical
checks. The nominal first resonances are harmonic orders 6.27 and 12.32.
The [capacity report](output/ieee30-full-equivalent.txt) contains only the output
of `run.jl`. The [comparison report](output/ieee30-full-comparison.txt) contains
both solves, decision vectors, the four posterior checks and the corner adversary.
The script `compare.jl` writes that report directly.

## Corner adversary

HPM omits dictionary vertices for these boxes because their dimensions exceed ten.
The nominal-sign adversary evaluates the following distinct predicted points.
Positive maxima violate the indicated constraints.

| Solution | Case | Corners | Maximum normalized violation | Constraint |
| --- | --- | ---: | ---: | --- |
| Source | Source | 88 | 0.002937446267 | IHD, harmonic 5, bus 30 |
| Equivalent | Equivalent | 134 | 0.009217241127 | IHD, harmonic 11, bus 7 |
| Equivalent | Source | 90 | 0.007534865169 | IHD, harmonic 11, bus 7 |
| Source | Equivalent | 134 | 0.012961966592 | IHD, harmonic 7, bus 30 |

Both own-case certificates fail this check. Both transfers also fail. A clean
1024-point random check did not establish feasibility on either box. The
comparison report gives every coordinate of each worst point. These are case
coordinate values, including signed positions for registered intervals.

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
of the four solution and case combinations.

A positive own-case violation refutes the reported feasibility of that solution
on its own box. A positive equivalent-on-source violation refutes transfer of the
equivalent solution to the source case. A positive source-on-equivalent violation
refutes transfer of the source solution to the independent equivalent box.
Cross-case violations alone do not refute either own-case certificate.
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
julia --project=examples/rhhc-variants examples/rhhc-variants/reduce.jl cases/ieee30-full.toml
julia --project=examples/rhhc-variants examples/rhhc-variants/test/runtests.jl
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/scan.jl cases/ieee30-full-equivalent.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/run.jl cases/ieee30-full-equivalent.toml
julia -t 16 --project=examples/rhhc-variants examples/rhhc-variants/compare.jl cases/ieee30-full.toml cases/ieee30-full-equivalent.toml
```

Reduction writes the case, coordinate data and compact report. Its temporary
evidence directory defaults to `/tmp/hpm-equivalent-step3`. A second argument
selects another directory. It recomputes subdivision, hulls, sharing, admission and
decorrelation. Sorted TOML output and a fixed table order make regeneration
deterministic. The source commit and LCM source fingerprint identify the input code.

The reduction report depends only on the generated reduction evidence.
`compare.jl` calls the runner's solve routine for each case with its nominal,
resonance and capacity settings. It captures the prepared network data through
that routine, then calls public HPM functions for all four posterior and corner
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
into `/tmp` and compare their bytes. They also load the generated case, compare
its nominal sequences with the source, check containment against its saved members
and exercise its coordinate interpolation away from nominal.

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

