# Uncertainty levels at the paper's harmonic orders

The equivalent box has **35 coordinates**, compared with **23 physical inputs**.
Exact lumping and measured sharing give 40 candidates, of which five radial widths
have endpoint changes below the admission threshold. The factor and scaled boxes have six coordinates.

The source is `ieee30-full-exact`: 23 physical coordinates, characteristic orders
5, 7, 11, 13, 17, 19, 23 and 25, exact squared-magnitude moments, and the Gaussian
factor. The equivalent, factor and scaled cases inherit it through `base`.
Their network and branch assignments come from the base. Installed lengths and converter population also
follow the base, together with its risk and angle tables. Each circuit has solid grounding and ideal transposition.

The generated files are:

| Level | Case | Data | Reduction evidence |
| --- | --- | --- | --- |
| Equivalent | `cases/ieee30-full-exact-equivalent.toml` | `data/equivalent/ieee30-full-exact.toml` | `output/ieee30-full-exact-equivalent-reduction.txt` |
| Factor | `cases/ieee30-full-exact-factor.toml` | `data/surrogate/ieee30-full-exact.toml` | `output/ieee30-full-exact-surrogate-reduction.txt` |
| Scaled | `cases/ieee30-full-exact-scaled.toml` | The same factor data | The same reduction evidence |

`reduce.jl` generates these files from the physical inputs through the current
builders and the pinned LCM extension. Raw numerical tables go to the supplied
evidence directory. The data identify their source case, revision, LCM pin and
regeneration command. The source commit identifies the last change to the case
file. Capacity reports record the current HPM revision. Uncommitted study changes
can affect a run.

## Exact lumping and the equivalent box

The public LCM blueprint replaces stranded conductors with plain tubes and
annular layers. Terminal names, materials, connections and placements are retained.
The direct reconstruction preserves the coaxial calculation. A second reconstruction
uses positive radial widths, operating DC resistance per meter and relative material
properties. For each conducting tube,

```math
\rho = R_{\mathrm{dc}}\,\pi(r_{\mathrm{out}}^2-r_{\mathrm{in}}^2).
```

Here resistivity is in Ω m, resistance per unit length is in Ω/m and radii are
in meters. The material reference-temperature conversion is retained. The core
and screen resistances combine their reference resistance and temperature effects.
Overhead heights combine construction and sag. Contact trefoils retain LCM's
spacing law as their outside radii change.

The physical source maps impose lay ordering, radial projection under diameter
and mass limits, and phase temperature at least ambient temperature. Sensitivities
are taken on these effective inputs, with feasible perturbations. The half step
is `1e-4` of the coordinate interval and the normalized zero band is `1e-6`.
The initial observations comprise nominal, both endpoints of each source axis,
and 256 Sobol points with seed 1: 303 points for this source box.

A coordinate represents a signed position: zero at nominal, minus one at the lower
hull bound and plus one at the upper hull bound. Interpolation uses the appropriate
width on each side of nominal. Sharing requires agreement of normalized profiles
at all 303 analysis points and 4096 independent verification points, within `1e-10`.
Each member of a group is compared with every other member. Values within 64 machine
spacings of nominal use position zero. Sharing retains the measured dependence
between those members. Separate coordinates discard their source correlations.

Admission tests each candidate's two endpoints, with other coordinates nominal,
through the extension at 50 Hz and every case harmonic frequency. The maximum
relative change across circuits, sequences, series impedance and shunt admittance
must reach `0.001`. The reduction report records every admitted coordinate and its
physical members. Because omitted members remain nominal, containment of admitted members does not
establish containment of the complete source physics after admission. Independent
widths also admit combinations of cumulative radii outside the physical source image.

Generated member records select the equivalent coordinates that each circuit reads.
The extension keys its stored results on those coordinates. Installed length remains
an input to HPM's exact π-model. The reducer also declares the effective inputs read
by its per-unit-length models.

## Recursive subdivision, hulls and containment

For conflicting sensitivity signs, `split_cell` tests splits on driving effective
coordinates. It minimizes the number of disagreeing derivative pairs in the two
children, then refines a bracketed derivative zero when the observations permit it.
Recursion stops at depth six. A still inconsistent cell uses an interior search.
The overhead phase-height split occurs near the reference temperature, where the
sensitivity to construction changes sign.

For each cell, the search closes the driving set over its effective constraints.
Independent monotone coordinates take oriented endpoints. Constrained corners are
projected to feasible order and aggregate limits by bisection. Dominated candidates
are removed. Coordinate lines and active faces are refined with golden-section
searches, including the face where inner and outer lays agree. A face search probes
nine values before 32 scalar iterations, followed by at most four refinement passes.
Interior searches start at feasible projected corners and the center, with tolerance
`1e-10` and 32 iterations. Finite endpoints receive 64 machine spacings of outward
rounding. Affine aliases share a search only when their driving sets and normalized
profiles agree within `1e-11`. Their bounds include the measured alias error.

Containment evaluates 4096 additional source-mapped Sobol points per circuit with
seed 3. Each quantity must lie in the interval of a cell containing its effective
inputs, without a comparison tolerance. Inconsistent cells receive another 4096
feasible projected Sobol samples per quantity. Generation refuses escapes. These
are numerical hulls with the recorded search and sample verification, rather than
interval-arithmetic proofs of global extrema.

The equivalent run has zero escapes at the 4096 verification points of each of
the eight circuits. Its largest normalized outward margin is
`-7.389644451904987e-15`. The largest sharing difference is
`1.1357581541915351e-12`.

The reduction also searches pairs of equivalent endpoints that share physical
drivers. Its joint residual measures how far the physical source can approach the
added combination. The report lists the five largest residuals.

## Family factors and square-root scaling

For each design and frequency, the positive-sequence quantities are

```math
R = \Re z,\qquad L = \frac{\Im z}{2\pi f},\qquad
C = \frac{\Im y}{2\pi f}.
```

Their units are Ω/m, H/m and F/m. Each factor interval is the hull of the ratios to
the same design's nominal value at the same frequency. It includes the four designs of
its family at 50, 250, 350, 550, 650, 850, 950, 1150 and 1250 Hz. The fundamental
supports the power-flow network. The other frequencies are the paper's harmonics.
ACSR and U33 each have one R, L and C factor, giving six coordinates.

The factor model multiplies the real and imaginary parts of nominal LCM series
impedance by the family R and L factors, and shunt susceptance by the C factor.
Shunt conductance retains its nominal value. `ExactPiModel` then computes the
branch parameters at the installed length. The positive, negative and zero sequences use the same factors
when requested. Positive sequence defines the fitted intervals.

The scaled model uses the source's nominal 50 Hz exact π branch parameters.
HPM's existing `ScaledBranchModel` applies

```math
r_h=k_R r_1\sqrt h,\qquad x_h=k_L x_1h,\qquad b_h=k_C b_1h.
```

This model supplies positive and negative sequences. A request for zero sequence
is refused. The characteristic orders used here need only positive and negative
sequences. Both surrogates read the same generated family intervals through the
builders. `hdata_at` constructs their nominal harmonic networks.

The interior search matters for screen loss: the approximate term
`Rs * Xm^2 / (Rs^2 + Xm^2)` has an interior maximum at `Rs = Xm`.
The report compares an LCM maximum with the geometric estimate
`Xm = mu0 * f * log(spacing / screen outer radius)`. LCM supplies every value used
in the hull search. The approximate expression only explains the observed maximum.

One family factor moves every design and frequency together, while its separate
R, L and C coordinates discard their physical correlations. Neither surrogate is
an inner or outer approximation of the source set. Their solutions need
checks on the physical model.

## Finding 3 and capacity checks

`findings.jl` reads `ieee30-full-exact-scaled` through `study`. It uses the same
scaled model and factor box as the level comparison. The former factor intervals
were extrema of 50 Hz branch-parameter ratios at predicted corners and 1024 random
points. They were a sampled inner envelope at one frequency. The generated intervals
instead use per-unit-length R, L and C ratios, with subdivision, interior searches
and containment checks across every design and case frequency.

The factor run has zero escapes in 32,768 source samples and 12,288 additional
interior samples. The largest normalized outward margins are
`-2.549795008242146e-6` for source samples and `-0.003105966575570989`
for interior samples. Each family interval contains 36 contributing ratios:
four designs at nine frequencies. The screen-loss check retains its interior
maximum for U33-400 at 250 Hz.

The previous finding printed its sampled intervals to four decimal places.
The generated endpoints below retain the precision of the data file.

| Family | Quantity | Previous finding, sampled at 50 Hz | Verified hull through harmonic 25 |
| --- | --- | --- | --- |
| ACSR | R | [0.8688, 1.2658] | [0.8623096825708542, 1.2682189163872188] |
| ACSR | L | [0.9976, 1.0039] | [0.9962618928782407, 1.005315726528468] |
| ACSR | C | [0.9982, 1.0018] | [0.9981387764625804, 1.0018618771950403] |
| U33 | R | [0.9045, 1.2309] | [0.9031342364062734, 1.2323618743221811] |
| U33 | L | [0.9810, 1.0438] | [0.9205434341153378, 1.1543234302904137] |
| U33 | C | [0.8324, 1.1052] | [0.8324072347466095, 1.1052049853966204] |

The R and L endpoints equal the archived h13 hulls from
`experimental/surrogates` at `6059a62`. Only the C endpoints change, by less
than `1.8e-10`. Their previous values were:

| Family | Previous verified C interval through harmonic 13 |
| --- | --- |
| ACSR | [0.9981387764988817, 1.0018618770161116] |
| U33 | [0.8324072347466105, 1.10520498539662] |

`compare.jl` runs the source and all three levels with the same solver settings.
It reports nominal capacity, robust hosting capacity, upper bounds, gaps,
iterations, termination and the solver's sign-off status. Each solution is checked
on each of the four models at 1024 random points. The corner adversary differentiates
every positive-limit constraint at nominal, chooses the adverse endpoint of each
coordinate from its derivative sign, and evaluates unique predicted corners.
Zero-sensitivity coordinates remain nominal. Positive own-case violations refute
feasibility of that reported solution. Positive transfer violations refute its
feasibility on the target model. Passing these checks is sampled evidence.

The comparison report is `output/ieee30-full-exact-comparison.txt`. Individual
reports for the three levels use their case stems. The source run is retained in
the comparison and its existing individual report is preserved. Finding 3 is in
`output/ieee30-full-exact-findings.txt`.

All four runs use eight threads. The nominal and final master problems report
`OPTIMAL`, and all four solutions have the solver status `signed_off = true`.
The table gives total harmonic hosting capacity in per unit. The comparison report
retains the unrounded values.

| Level | Coordinates | Nominal | Reported lower bound | Upper bound | Relative gap | Iterations |
| --- | --- | --- | --- | --- | --- | --- |
| Source | 23 | 0.099250 | 0.076917 | 0.076943 | 0.000331 | 3 |
| Equivalent | 35 | 0.099250 | 0.076728 | 0.076762 | 0.000447 | 6 |
| Factor | 6 | 0.099250 | 0.072568 | 0.072570 | 2.77e-05 | 3 |
| Scaled | 6 | 0.078624 | 0.048238 | 0.048251 | 0.000286 | 6 |

The following tables give the maximum normalized constraint violation. Rows select
the solution and columns select the evaluation model. A positive value is a violation.
The posterior check uses 1024 points per pair, seed 2, and tolerance `1e-4`.

| Posterior check | Source | Equivalent | Factor | Scaled |
| --- | --- | --- | --- | --- |
| Source | -0.0311009 (0%) | -0.0320371 (0%) | +0.329727 (29.4922%) | +7.78084 (100%) |
| Equivalent | -0.0349415 (0%) | -0.0334533 (0%) | +0.324456 (26.6602%) | +7.772 (100%) |
| Factor | -0.0621701 (0%) | -0.0548741 (0%) | -0.0197281 (0%) | +7.42783 (100%) |
| Scaled | +2.4929 (100%) | +2.39921 (100%) | +3.425 (100%) | -0.0482512 (0%) |

Percentages give the fraction of sampled points with a violation above tolerance.
The next table gives the corner maximum, with the number of unique corners in parentheses.

| Corner check | Source | Equivalent | Factor | Scaled |
| --- | --- | --- | --- | --- |
| Source | -0.000654548 (139) | +0.0305246 (222) | +0.385599 (42) | +8.17649 (56) |
| Equivalent | -0.00152759 (138) | +0.000492392 (224) | +0.351036 (42) | +8.16854 (56) |
| Factor | -0.023656 (140) | -0.023056 (221) | -5.54876e-05 (41) | +7.96534 (59) |
| Scaled | +2.72212 (137) | +2.85381 (228) | +3.95246 (40) | -0.000572351 (47) |

The equivalent solution fails its own corner check at IHD(19, 5), with violation
`0.0004923919623234241`, despite solver sign-off. Its reported lower bound
has a detected feasibility failure. The source and factor solutions pass their own
sample and corner checks. Both the equivalent and factor solutions pass these checks
on the physical source model.

The scaled solution passes its own checks but violates
the physical constraints at every sampled point. Its largest physical corner violation
is `2.7221245917169345` at IHD(23, 3). A smaller capacity alone does not establish
physical feasibility, because the harmonic current allocation also changes.

The comparison corner adversary uses the candidate currents and normalized constraint
derivatives. Finding 1 retains its original HPM network-dictionary corner construction.
Their corner counts differ.

The findings rerun uses HPM revision `4aec30b05762160dac375a747807c2cf079b6a19`.
The previous report used `ed4e0c95a3d43c4d0fb5bebf87fa7bba927afa61`.
The LCM pin is unchanged. The interval table above gives every factor endpoint
before and after the change. The report prints those endpoints to four decimal
places and names the families `ACSR` and `U33`, previously `line` and `cable`.
The remaining changed numbers in `ieee30-full-exact-findings.txt` are:

| Report field | Previous | Rerun |
| --- | --- | --- |
| Threads | 10 | 8 |
| With-corners elapsed time, s | 180 | 527 |
| Without-corners elapsed time, s | 253 | 784 |
| Binding IHD(23, 3) violation at the refined point | +4.44e-16 | -3.33e-16 |
| Scaled capacity over the factor box | 0.055999 | 0.048238 |
| Scaled solution, worst sampled physical violation | +6.50e+00 | +2.49e+00 |

The text diff shows no further numerical changes. The physical capacities
remain `0.076917` with corners and `0.076980` without corners. The nominal scaled
and physical capacities remain `0.078624` and `0.099250`. The nominal scaled
solution's physical violation remains `+8.68e+00`. The scaled solution over the
factor box still violates at all 1024 sampled physical points. Elapsed times come
from runs with different thread counts on a shared machine.

## Reproduction and tests

Set `TMPDIR` to the chosen scratch directory before launching Julia. Set
`SHHC_HMAX=25` before loading the builders: the angle tables are fixed at load time.
From the repository root, with an existing study environment:

```sh
SHHC_HMAX=25 julia -t 4 --project=examples/rhhc-variants \
    examples/rhhc-variants/reduce.jl cases/ieee30-full-exact.toml "$TMPDIR/reduction"
SHHC_HMAX=25 julia -t 4 --project=examples/rhhc-variants - <<'JL'
include("examples/rhhc-variants/compare.jl")
paths = ["cases/ieee30-full-exact.toml",
    "cases/ieee30-full-exact-equivalent.toml",
    "cases/ieee30-full-exact-factor.toml", "cases/ieee30-full-exact-scaled.toml"]
mktempdir() do directory
    reports = [joinpath(directory, "source.txt");
        [joinpath(STUDIES, "output", splitext(basename(path))[1] * ".txt")
            for path in paths[2:end]]]
    open(joinpath(STUDIES, "output", "ieee30-full-exact-comparison.txt"), "w") do io
        comparison(io, reports, paths...)
    end
end
JL
SHHC_HMAX=25 julia -t 4 --project=examples/rhhc-variants examples/rhhc-variants/findings.jl
SHHC_HMAX=13 JULIA_NUM_THREADS=4 julia --project=examples/rhhc-variants \
    examples/rhhc-variants/test/runtests.jl
```

The study suite launches the exact-table checks in a child process. By default,
the generators read the saved search evidence in
[`data/reduction/ieee30-full-exact/`](data/reduction/ieee30-full-exact/README.md)
and reproduce the cases, data and reduction reports byte for byte. The generators
read the search results on this path. Fresh containment checks use 256 mapped Sobol
points per design, with seed 3 and the original zero tolerance. The equivalent
check uses the saved cell bounds. The factor check uses the saved scalar hulls
and the family hulls. The factor checks also retain the nominal point, source-axis
endpoints and 256 analysis points.

Equivalent sharing uses fresh observations at the analysis and verification
points. Its source dependencies come from the saved sensitivity results. The
full derivation recomputes those dependencies.

Blueprint and sequence reconstruction, effective maps, subdivision metadata,
sharing, independent combinations, nominal identity, factor application and
screen-loss assertions remain in the default suite. It also checks the recorded
counts and margins of the full interior containment run. These stored results
record the original derivation. The 256-point checks are fresh model evaluations.

Full rederivation is explicit. The `reduce.jl` command above always performs the
searches, admission and joint-extreme searches, with 4096 source verification points
and the additional samples inside inconsistent cells. It writes the numerical
tables and `factor-verification.toml` to its evidence directory. To refresh the
stored evidence, supply `examples/rhhc-variants/data/reduction/ieee30-full-exact`
as that directory.

The full test path retains every previous search and containment assertion with
its original tolerance. Run it with:

```sh
RHHC_REDERIVE=1 SHHC_HMAX=13 JULIA_NUM_THREADS=4 \
    julia --project=examples/rhhc-variants examples/rhhc-variants/test/runtests.jl
```

The child inherits `RHHC_REDERIVE` and sets `SHHC_HMAX=25` before loading the
builders. This test path writes fresh evidence and regenerated files under
`TMPDIR`. It does not overwrite the stored study outputs.
