# Robust harmonic hosting capacity: plugging in LineCableModels

This guide is for plugging LineCableModels.jl into the IEEE 30-bus test case of the robust harmonic hosting capacity paper. It assumes no background in robust optimization. Everything lives on the branch `feature/rhhc-cutting-plane` of HarmonicPowerModels.jl.

## The problem in five lines

1. Many converters at many buses inject harmonic currents. Their phase angles are random, so the harmonic voltages they cause are random too.
2. The *hosting capacity* is the largest harmonic current per converter such that every voltage and current limit holds, in probability, at every harmonic. That is the stochastic problem (SHHC) of the previous paper.
3. The cable parameters are not known exactly; your tolerances define a box of possible cables. Near a resonance the voltages are very sensitive to them: on this case, a few percent of cable capacitance moves the resonance across a harmonic.
4. The *robust* hosting capacity has to hold for every cable in the box.
5. The method guesses a capacity, searches the box for the cable that breaks a limit most, adds that cable to the problem, and repeats until no cable in the box breaks a limit.

## Where LineCableModels fits

Step 5 needs, again and again, the network at a point `p` of the box and at a harmonic frequency `f`. For every cable branch, HarmonicPowerModels asks your code for

```julia
z1, y1 = zy(p, f, sequence)    # positive-sequence series impedance [Ω/m] and shunt admittance [S/m]
```

It then builds the exact π-model of the cable at its length, Z = Z_c·sinh(γℓ) and Y/2 = tanh(γℓ/2)/Z_c, and puts it into the harmonic network.

That callback is the whole interface. You never see the optimization; HarmonicPowerModels never sees the cable.

## Quick start

```
git clone -b feature/rhhc-cutting-plane https://github.com/timmyfaraday/HarmonicPowerModels.jl
julia> ] activate <a fresh env>
julia> ] develop <path to HarmonicPowerModels.jl>
julia> ] develop <path to LineCableModels.jl>     # main, v0.1.0
julia> ] add PowerModels JuMP Ipopt MosekTools     # or Clarabel instead of MosekTools
julia -t 10 --project=<env> HarmonicPowerModels.jl/examples/rhhc/ieee30_lcm.jl
```

`ieee30_lcm.jl` runs as is. It uses an 18/30 kV XLPE 1x630 Al cable in trefoil (my placeholder design) for all 22 cables, and a box over the XLPE permittivity, the insulation thickness, the core radius and the earth resistivity. It prints the nominal and robust capacities and writes them to `output/ieee30_lcm.txt`. Use `-t` with your number of cores: the search over the box is multithreaded.

With the placeholder design (10 threads, about 20,000 LineCableModels solves, under a minute including compilation):

| result | value |
|---|---|
| resonance at the nominal point | h = 6.86 |
| nominal stochastic capacity | 0.0502 pu |
| robust capacity | 0.0437 pu (−13 %; 2 iterations, 13 s; gap 4.6e-4) |
| robust / nominal per harmonic | 5th 0.83, 7th 0.90, 11th 0.87, 13th 0.91 |
| a-posteriori check, 1024 random points | no violation |

## What to edit

Only two blocks in `ieee30_lcm.jl`, both marked:

1. **`cable_zy(p, f)`.** Your cable design, formation, bonding and earth. It must return `(z1, y1)` per meter at frequency `f` for the parameter values in the named tuple `p`. The starter takes z1 = zs − zm and y1 = ys − ym from the 3×3 matrices, with the sheaths grounded at both ends and ideal transposition.
2. **`box`.** One entry per uncertain input: `name => (nominal, lower, upper)`, with `:log` as a fourth element for inputs that span orders of magnitude, like earth resistivity. Every name in the box must be a field that `cable_zy` reads from `p`. Inputs not in the box keep the nominal values in `CABLE_NOMINAL`.

Two common extensions:

- **Different cable types per feeder.** Give each type its own builder and its own box names (e.g., `eps_xlpe_A`, `eps_xlpe_B`), and point each branch at its builder in the `ExactPiModel` dictionary.
- **Uncertain installed length** (your suggestion). Replace `len = 1e3 * km` by a function of `p`, e.g., `len = p -> 1e3 * km * p.len_17`, and add `:len_17 => (1.0, 0.95, 1.05)` to the box for branch 17. The length enters the π-model, not your per-meter data.

## The test case

IEEE 30-bus system, scenario `inductive_cable33` (`examples/shhc/ieee30_scenarios.jl`), with:

- base 100 MVA; voltage levels 132 kV, 33 kV and 11 kV;
- all 33 kV lines as cables, with the route length recovered from the original line reactance;
- the 132 kV lines as overhead lines without charging;
- 21 harmonic units (converters) at the first 21 loads.

The cables (33 kV, from the case data):

| branch | buses | km | branch | buses | km | branch | buses | km |
|---|---|---|---|---|---|---|---|---|
| 17 | 12–14 | 8.0 | 25 | 10–20 | 6.5 | 32 | 23–24 | 8.4 |
| 18 | 12–15 | 4.1 | 26 | 10–17 | 2.6 | 33 | 24–25 | 10.2 |
| 19 | 12–16 | 6.2 | 27 | 10–21 | 2.3 | 34 | 25–26 | 11.8 |
| 20 | 14–15 | 6.2 | 28 | 10–22 | 4.7 | 35 | 25–27 | 6.5 |
| 21 | 16–17 | 6.0 | 29 | 21–22 | 0.7 | 37 | 27–29 | 12.9 |
| 22 | 15–18 | 6.8 | 30 | 15–23 | 6.3 | 38 | 27–30 | 18.8 |
| 23 | 18–19 | 4.0 | 31 | 22–24 | 5.6 | 39 | 29–30 | 14.1 |
| 24 | 19–20 | 2.1 | | | | | | |

The other network elements:

- **Overhead lines (132 kV, case data):** branches 1–10, 40 and 41, from 16.5 km to 87.1 km. You can model them as well, with a bare-conductor design in air. Note that the case removes their charging on purpose; a physical model brings it back, which changes the resonances.
- **Transformers:** branches 11–16 and 36. They are not touched.
- **Harmonic units at buses:** 2, 3, 4, 5, 7, 8, 10, 12, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 26, 29, 30.

With the placeholder data of the case (r = 0.10, x = 0.11 Ω/km, c = 0.30 µF/km), the network has parallel resonances near h = 6.6 and h = 13.7. The script prints the resonances of your cables, which will differ.

## Performance

- **One evaluation of the network at a point `p`** needs your cables at the 4 harmonics (250, 350, 550 and 650 Hz). A robust solve evaluates a few thousand points.
- **The resonance scan** at the start evaluates the nominal cable at about 800 frequencies between 50 Hz and 812.5 Hz.
- **The cache in the script** keys on `(p, f)`. All 22 cables share one design, so one LineCableModels solve serves all of them. This matters: without it, each evaluation costs 22 × 4 solves.
- **The calls come from several threads.** The cache is locked. If `compute` is not thread-safe, run with `-t 1` (slower) and tell me; I can serialize the calls.
- **Batching would be faster.** If a frequency vector costs about the same as one frequency in LineCableModels, computing all harmonic frequencies of a point `p` in one call would help. The callback is per frequency today, but it can change.

## Checks before trusting a run

1. **Per-unit sanity.** At 50 Hz, `branch_data(model, b, nominal_point(box), 1.0, :positive)` returns `br_r`, `br_x` and `b_fr` in per unit on Z_base = 33²/100 = 10.89 Ω. Compare them with the case's assumed values times the length.
2. **Zero-width box.** Set lower = nominal = upper for every input. The robust capacity must equal the nominal stochastic capacity, in one iteration. This is also how the script computes the nominal capacity. Calling `solve_shhc` on the replicated data instead would scale the cables' 50 Hz data by r·√h, x·h and b·h rather than use your values at each harmonic. On this case that changes the 13th-harmonic capacity by a factor of two, so the frequency dependence of your data matters.
3. **Resonances.** Physical cable data should still place the first resonance between the 5th and the 13th harmonic. If it does not, check the units of `y1` (S/m, imaginary part ωC).
4. **Signed off.** The output line says `signed off`, and the a-posteriori check reports no violated points.

## What is not covered yet

- **Zero sequence (triplen harmonics).** It is not used yet. When you have z0 = zs + 2zm and y0 = ys + 2ym, return them for `sequence == :zero`; the transformer vector groups and the earthing still need agreeing.
- **The fundamental power flow** uses the nominal cables only. The cable uncertainty acts on the harmonics.
- **Derivatives.** LineCableModels does not need to be generic for dual numbers or intervals: the search uses finite differences.

## Files

| file | what |
|---|---|
| `ieee30_lcm.jl` | starter: LineCableModels cables in the IEEE30 robust case |
| `ieee30_robust.jl` | the paper study with an analytic cable model (scale factors on R, L, C and lengths) |
| `lcm_network_model.jl`, `data/hello_world.m` | the small 33 kV hello-world network with LineCableModels |
| `data/lc_minimal.m` | the minimal L–C example (tests) |
| `output/` | results of the study runs |
| `src/robust/` (package) | the method: `network_model.jl` (`ParameterBox`, `ExactPiModel`), `sampling.jl`, `oracle.jl`, `cutting_set.jl` (`solve_rhhc`, `posterior_check`) |
