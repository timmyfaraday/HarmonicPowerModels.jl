# HarmonicPowerModels.jl — Code Audit

**Branch:** `shhc-pce` (created from `main`, HEAD `c332226`)  
**Audit date:** 2026-05-16

---

## A. Project.toml

```toml
name = "HarmonicPowerModels"
uuid = "56365af1-a032-4c1d-861b-3c9f15bdb54e"
authors = ["Hakan Ergun", "Frederik Geth", "Tom Van Acker"]
version = "0.2.1"

[deps]
InfrastructureModels = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
PowerModels = "c36e90e8-916a-50a6-bd94-075b64ef4655"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
SignalDecomposition = "11a47235-7b84-4c7c-b885-fc3e2a9cf955"

[compat]
InfrastructureModels = "~0.6, ~0.7"
Interpolations = "0.15"
JuMP = "1.23"
MathOptInterface = "1.31"
PowerModels = "0.21"
ProgressMeter = "1.10.0"
SignalDecomposition = "1.1"
julia = "1.9"

[extras]
Clarabel = "61c947e1-3e6d-4ee4-985a-eec8c727bd6e"
Dierckx = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test", "Ipopt", "Clarabel", "Dierckx"]
```

---

## B. Package entry point (`src/HarmonicPowerModels.jl`)

```julia
# using
using ProgressMeter

# import
import JuMP
import MathOptInterface
import PowerModels
import InfrastructureModels
import SignalDecomposition
import Interpolations
import InfrastructureModels: replicate   # overwrite target

# include (in order)
include("core/base.jl")
include("core/types.jl")
include("core/constraint_template.jl")
include("core/data.jl")
include("core/variable.jl")
include("form/iv.jl")
include("prob/hopf.jl")
include("prob/hpf.jl")
include("prob/hhc.jl")
include("util/imp.jl")
include("util/init.jl")
include("util/ref.jl")
include("util/sol.jl")
include("util/xfmr_magn.jl")

# export
export BASE_DIR
export dHHC_NLP, dHHC_SOC
export replicate
export solve_hpf, solve_hopf, solve_hhc
export calculate_pos_seq_harmonic_impedance
```

---

## C. Type hierarchy (`src/core/types.jl`)

| name | kind | supertype |
|------|------|-----------|
| `dHHC_NLP` | mutable struct | `_PMs.AbstractIVRModel` (via `<: _PMs.AbstractIVRModel`) |
| `dHHC_SOC` | mutable struct | `_PMs.AbstractIVRModel` |

Both use the `_PMs.@pm_fields` macro for their fields. There are no abstract types defined in `types.jl` itself. The inheritance chain is:

```
dHHC_NLP  <:  PowerModels.AbstractIVRModel  <:  PowerModels.AbstractACPowerModel  <:  ...
dHHC_SOC  <:  PowerModels.AbstractIVRModel  <:  ...
```

---

## D. Variable functions (`src/core/variable.jl` + `src/form/iv.jl`)

### Leaf-level scalar variable declarations (`variable.jl`)

| function name | signature | JuMP variable symbol(s) |
|---|---|---|
| `variable_fairness_principle` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cmh` (maximin), `:fh` (Kalai-Smorodinsky) |
| `variable_bus_voltage_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:vr` |
| `variable_bus_voltage_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:vi` |
| `variable_branch_current_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cr` |
| `variable_branch_current_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:ci` |
| `variable_branch_series_current_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:csr` |
| `variable_branch_series_current_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:csi` |
| `variable_filter_current_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:crf` |
| `variable_filter_current_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cif` |
| `variable_xfmr_voltage_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:vrx` |
| `variable_xfmr_voltage_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:vix` |
| `variable_xfmr_voltage_excitation_real` | `(pm::AbstractPowerModel; nw, bounded, report, epsilon)` | `:erx` |
| `variable_xfmr_voltage_excitation_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report, epsilon)` | `:eix` |
| `variable_xfmr_current_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:crx` |
| `variable_xfmr_current_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cix` |
| `variable_xfmr_current_series_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:csrx` |
| `variable_xfmr_current_series_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:csix` |
| `variable_xfmr_current_magnetizing_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cmrx` |
| `variable_xfmr_current_magnetizing_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cmix` |
| `variable_gen_current_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:crg` |
| `variable_gen_current_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cig` |
| `variable_load_current_real` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:crd` |
| `variable_load_current_imaginary` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cid` |
| `variable_load_current_magnitude` | `(pm::AbstractPowerModel; nw, bounded, report)` | `:cmd` |

### Compound dispatchers (`form/iv.jl`)

| function name | dispatches to |
|---|---|
| `variable_bus_voltage` | `variable_bus_voltage_real`, `variable_bus_voltage_imaginary` |
| `variable_branch_current` | `:cr`, `:ci`, `:csr`, `:csi` |
| `variable_filter_current` | `:crf`, `:cif` |
| `variable_xfmr_voltage` | `:vrx`, `:vix`, `:erx`, `:eix` |
| `variable_xfmr_current` | `:crx`, `:cix`, `:csrx`, `:csix`, `:cmrx`, `:cmix` |
| `variable_gen_current` | `:crg`, `:cig` |
| `variable_load_current` | `:crd`, `:cid`, `:cmd` (`:cmd` only for `nw ≠ fundamental`) |

---

## E. Constraint functions

### Template level (`src/core/constraint_template.jl`)

| function name | signature | description |
|---|---|---|
| `constraint_voltage_ref_bus` | `(pm::AbstractPowerModel, i::Int; nw)` | Sets vref=1 (fundamental) or 0 (harmonics), calls IV dispatch |
| `constraint_voltage_rms_limit` | `(pm::AbstractPowerModel, i::Int)` | Reads `vminrms`, `vmaxrms` from ref; dispatches to IV form |
| `constraint_voltage_rms_limit` | `(pm::dHHC_SOC, i::Int)` | SOC version: reads `vmaxrms`, `vm` from fundamental; dispatches to SOC form |
| `constraint_voltage_thd_limit` | `(pm::AbstractPowerModel, i::Int)` | Reads `thdmax` from ref; dispatches |
| `constraint_voltage_thd_limit` | `(pm::dHHC_SOC, i::Int)` | SOC version: reads `thdmax` + fundamental `vm` |
| `constraint_voltage_ihd_limit` | `(pm::AbstractPowerModel, i::Int; nw)` | Reads `ihdmax`; skips if fundamental |
| `constraint_voltage_ihd_limit` | `(pm::dHHC_SOC, i::Int; nw)` | SOC version: reads `ihdmax` + fundamental `vm` |
| `constraint_current_balance` | `(pm::AbstractPowerModel, i::Int; nw)` | Gathers bus arcs, xfmr arcs, filters, gens, loads, shunts; dispatches |
| `constraint_current_rms_limit` | `(pm::AbstractPowerModel, b::Int)` | Branch RMS: reads `c_rating`; dispatches to NLP form |
| `constraint_current_rms_limit` | `(pm::dHHC_SOC, b::Int)` | SOC version: reads `c_rating`, `cm_fr`, `cm_to` |
| `constraint_fairness_principle` | `(pm::AbstractPowerModel; nw)` | Gathers sorted load_ids; dispatches |
| `constraint_active_filter_current` | `(pm::AbstractPowerModel, f::Int)` | Only for active (`a/p == "a"`) filters; dispatches |
| `constraint_gen_current` | `(pm::AbstractPowerModel, g::Int; nw)` | Reads `inf`, `gsc`, `bsc`; dispatches if not infinite bus |
| `constraint_gen_current_rms_limit` | `(pm::AbstractPowerModel, g::Int)` | Reads `c_rating`; dispatches NLP form |
| `constraint_gen_current_rms_limit` | `(pm::dHHC_SOC, g::Int)` | Reads `c_rating`, fundamental `cm`; dispatches SOC form |
| `constraint_load_current` | `(pm::AbstractPowerModel, l::Int; nw)` | Fundamental: constant power; harmonic: current angle |
| `constraint_load_power` | `(pm::AbstractPowerModel, l::Int; nw)` | Fundamental: constant power; harmonic: constant current (uses `multiplier`) |
| `constraint_xfmr_core_magnetization` | `(pm::AbstractPowerModel, x::Int; nw)` | Non-linear interpolation if `Hᴵ` set, else zero |
| `constraint_xfmr_core_voltage_drop` | `(pm::AbstractPowerModel, x::Int; nw)` | Reads `xsc`; dispatches |
| `constraint_xfmr_core_voltage_phase_shift` | `(pm::AbstractPowerModel, x::Int; nw)` | Reads `tr`, `ti`; dispatches |
| `constraint_xfmr_core_current_balance` | `(pm::AbstractPowerModel, x::Int; nw)` | Reads `tr`, `ti`, `gsh`; dispatches |
| `constraint_xfmr_winding_config` | `(pm::AbstractPowerModel, x::Int; nw)` | Reads `r1`,`r2`,`re1`,`re2`,`xe1`,`xe2`,`gnd1`,`gnd2`; dispatches per winding |
| `constraint_xfmr_winding_current_balance` | `(pm::AbstractPowerModel, x::Int; nw)` | Reads `r1`,`r2`,`cnf1`,`cnf2`; dispatches per winding |
| `constraint_xfmr_current_rms_limit` | `(pm::AbstractPowerModel, x::Int)` | Reads `c_rating`; dispatches for each winding index |
| `constraint_xfmr_current_rms_limit` | `(pm::dHHC_SOC, x::Int)` | SOC version: reads `ctm_fr`, `ctm_to` |

### Implementation level (`src/form/iv.jl`) — selected key forms

| function name | description |
|---|---|
| `constraint_voltage_ref_bus(pm::AbstractIVRModel, n, i, vref)` | `vr == vref`, `vi == 0` |
| `constraint_voltage_rms_limit(pm::AbstractIVRModel, i, vminrms, vmaxrms)` | Quadratic RMS bounds over all harmonics |
| `constraint_voltage_rms_limit(pm::dHHC_SOC, i, vmaxrms, vmfund)` | `[sqrt(vmaxrms²-vmfund²); vcat(vr,vi)] ∈ SOC` |
| `constraint_voltage_thd_limit(pm::AbstractIVRModel, i, thdmax)` | Quadratic THD |
| `constraint_voltage_thd_limit(pm::dHHC_SOC, i, thdmax, vmfund)` | `[thdmax*vmfund; vcat(vr,vi)] ∈ SOC` |
| `constraint_voltage_ihd_limit(pm::AbstractIVRModel, n, i, ihdmax)` | Quadratic IHD |
| `constraint_voltage_ihd_limit(pm::dHHC_SOC, n, i, ihdmax, vmfund)` | `[ihdmax*vmfund; vcat(vr,vi)] ∈ SOC` |
| `constraint_current_balance(pm::AbstractIVRModel, n, i, ...)` | KCL: branches + xfmrs + filters = gens - loads - shunts |
| `constraint_current_rms_limit(pm::AbstractIVRModel, f_idx, t_idx, c_rating)` | Quadratic RMS for branch |
| `constraint_current_rms_limit(pm::dHHC_SOC, f_idx, t_idx, c_rating, cm_fr, cm_to)` | SOC form for harmonic-only currents |
| `constraint_fairness_principle(pm::AbstractIVRModel, n, load_ids)` | Dispatches by `pm.data["principle"]` |
| `constraint_load_constant_power` | `pd = vr*crd + vi*cid`, `qd = vi*crd - vr*cid` |
| `constraint_load_current_angle` | `cmd*sind(aref) == cid`, `cmd*cosd(aref) == crd` |
| `constraint_load_constant_current` | `crd = mult*fund_crd`, `cid = mult*fund_cid` |

---

## F. Data functions (`src/core/data.jl`)

| function name | signature | description |
|---|---|---|
| `ihd_current_limit_ieee_519` | `(voltage, i_ratio, harmonic)` | Returns per-unit IHD current limit per IEEE 519-2022 table look-up |
| **`_HPM.replicate`** | `(data::Dict{String,Any}; bus_id, H, xfmr_magn)` | **Builds the multi-network harmonic dict for the deterministic HHC problem.** Computes branch/xfmr current ratings, applies THD/IHD standards, calls `_PMs.replicate`, rescales branch/gen/xfmr impedances per harmonic, sets transformer sequence shifts |

Constants also defined here:
- `ihd_limits::Dict` — IHD voltage limit tables (IEC/IEEE standards)
- `thd_limits::Dict` — THD voltage limit scalars

---

## G. HHC problem structure

### Top-level solve function for `dHHC_SOC`

```julia
solve_hhc(hdata, model_type::Type, hhc_optimizer, hpf_optimizer; kwargs...)
```

(Two-argument form: takes both `hhc_optimizer` and `hpf_optimizer`.)

### Build function passed to `solve_model`

```julia
build_hhc  # same name for both dHHC_NLP and dHHC_SOC, dispatched by model type
```

Called as: `_PMs.solve_model(hdata, model_type, hhc_optimizer, build_hhc; ...)`

### Inside `build_hhc(pm::dHHC_SOC)` — in order

**Pre-step:** `JuMP.add_bridge(pm.model, _MOI.Bridges.Constraint.SOCtoNonConvexQuadBridge)`

**Variable calls** (loop over `nw_ids(pm)`, only for `n ≠ fundamental`):
1. `variable_fairness_principle(pm, nw=n, bounded=true)`
2. `variable_bus_voltage(pm, nw=n, bounded=true)`
3. `variable_xfmr_voltage(pm, nw=n, bounded=true)`
4. `variable_branch_current(pm, nw=n, bounded=true)`
5. `variable_xfmr_current(pm, nw=n, bounded=true)`
6. `variable_filter_current(pm, nw=n, bounded=true)`
7. `variable_load_current(pm, nw=n, bounded=true)`
8. `variable_gen_current(pm, nw=n, bounded=true)`

**Objective:**
- `objective_maximum_hosting_capacity(pm)`

**Overall constraint loops** (commented out in `dHHC_SOC`, i.e., NOT active):
- `constraint_voltage_rms_limit` — commented out
- `constraint_voltage_thd_limit` — commented out
- `constraint_current_rms_limit` (branch) — commented out
- `constraint_xfmr_current_rms_limit` — commented out

**Harmonic constraint calls** (loop over `nw_ids(pm)`, only for `n ≠ fundamental`):
1. `constraint_fairness_principle(pm, nw=n)`
2. `constraint_voltage_ref_bus(pm, i, nw=n)` (for each ref bus, conditional on `fix_refbus_angle`)
3. `constraint_current_balance(pm, i, nw=n)` (for each bus)
4. `constraint_voltage_ihd_limit(pm, i, nw=n)` (for each bus)
5. `_PMs.constraint_current_from(pm, b, nw=n)` (for each branch)
6. `_PMs.constraint_current_to(pm, b, nw=n)` (for each branch)
7. `_PMs.constraint_voltage_drop(pm, b, nw=n)` (for each branch)
8. `constraint_xfmr_core_magnetization(pm, x, nw=n)` (for each xfmr)
9. `constraint_xfmr_core_voltage_drop(pm, x, nw=n)`
10. `constraint_xfmr_core_voltage_phase_shift(pm, x, nw=n)`
11. `constraint_xfmr_core_current_balance(pm, x, nw=n)`
12. `constraint_xfmr_winding_config(pm, x, nw=n)`
13. `constraint_xfmr_winding_current_balance(pm, x, nw=n)`
14. `constraint_load_current(pm, l, nw=n)` (for each load)

> **Note:** In `dHHC_SOC`, the overall (cross-harmonic) RMS and RMS current constraints are **completely commented out**. Only the per-harmonic IHD limit is enforced. The `dHHC_NLP` version does include THD, RMS, and branch/gen/xfmr RMS limits.

---

## H. Data field names

| quantity | exact `ref()` / dict path |
|---|---|
| **IHD voltage limit** | `_PMs.ref(pm, nw, :bus, i, "ihdmax")` — set by `replicate` from `ihd_limits[std][nh]` |
| **THD voltage limit** | `_PMs.ref(pm, fundamental(pm), :bus, i, "thdmax")` — set by `replicate` from `thd_limits[std]` |
| **Current / ampacity limit (branch)** | `branch["c_rating"]` — set in `replicate` as `rate_a / sqrt(3) / vmmin` |
| **Current / ampacity limit (xfmr)** | `xfmr["c_rating"]` — set in `replicate` as `rateA / baseMVA / vmmin` |
| **Current / ampacity limit (gen/load)** | `gen["c_rating"]` / `load["c_rating"]` |
| **Bus reference voltage (fundamental)** | `_PMs.ref(pm, fundamental(pm), :bus, i, "vm")` (from HPF results written into hdata) |
| **Bus vmin/vmax** | `bus["vmin"]`, `bus["vmax"]` — from data; harmonic `vmin` set to `0.0` in `replicate` |
| **Harmonic unit entries** | Component name: **`:load`** — `_PMs.ids(pm, :load, nw=n)`, accessed as `_PMs.ref(pm, nw, :load, l)` |

---

## I. Variable symbols

| quantity | exact `:symbol` |
|---|---|
| Real part of bus voltage | `:vr` |
| Imaginary part of bus voltage | `:vi` |
| Lifted voltage-squared auxiliary | **none** — no `w` or `wr`/`wi` variable exists; SOC constraints use `vr`/`vi` directly |
| Load current magnitude (hosting capacity budget) | `:cmd` (one per load per harmonic network) |
| Fairness variable — maximin | `:cmh` (one per harmonic network) |
| Fairness variable — Kalai-Smorodinsky | `:fh` (one per harmonic network; bounded [0,1]) |

---

## J. Test data

### `test/hhc.jl` — NLP Industrial Network

- **Data file:** `test/data/matpower/industrial_network_hhc.m`
- **Solve function:** `HPM.solve_hhc(hdata, dHHC_NLP, solver_nlp)` (single-optimizer form)
- **Fields accessed in `result["solution"]`:**
  - `["nw"]["$nh"]["bus"][nb]["vm"]`
  - `["nw"]["$nh"]["bus"][nb]["vr"]`, `["vi"]`
  - `["nw"]["$nh"]["branch"]["$nb"]["csr_fr"]`, `["csi_fr"]`, `["cr_fr"]`, `["ci_fr"]`, `["cr_to"]`, `["ci_to"]`
  - `["nw"]["$nh"]["xfmr"]["$nx"]["erx"]`, `["eix"]`, `["vrx_fr"]`, `["vix_fr"]`, `["vrx_to"]`, `["vix_to"]`, `["cmrx"]`, `["cmix"]`, `["crx_fr"]`, `["cix_fr"]`, `["crx_to"]`, `["cix_to"]`, `["csrx_fr"]`, `["csix_fr"]`, `["csrx_to"]`, `["csix_to"]`

### `test/hhc.jl` — SOC Industrial Network

- **Data file:** `test/data/matpower/industrial_network_hhc.m`
- **Solve function:** `HPM.solve_hhc(hdata, dHHC_SOC, solver_soc, solver_nlp)` (two-optimizer form)
- **Fields accessed:** same bus/branch/xfmr fields as NLP; for SOC the fundamental `vm` comes from `hdata["nw"]["1"]["bus"][nb]["vm"]` (not from `result`)

### `test/hhc.jl` — Equivalence test

- **Data file:** `test/data/matpower/industrial_network_hhc.m`
- **Solve functions:** both `dHHC_NLP` and `dHHC_SOC` variants

---

## K. Large-case content

**No file** under `src/` or `examples/` references a "large", "transmission", or 1880-bus case. The largest test network referenced is the **IEEE 30-bus system** (`test/data/matpower/30bus_network_hhc.m`), used only in the `examples/dhhc_pes_gm/ieee_30_bus_system_hhc.jl` example script.

---

## L. Stochastic content

**None found.** A comprehensive grep of all `.jl` files under `src/` and `examples/` for:

```
PolyChaos, Beta01OrthoPoly, Tensor, computeSP2, shhc, SHHC, pce, PCE
```

returned **zero matches**. There is no stochastic HHC (sHHC) or polynomial chaos expansion (PCE) code anywhere in the repository.

---

## M. Julia version

| field | value |
|---|---|
| Julia compat string | `julia = "1.9"` |
| Package version | `version = "0.2.1"` |

---

## N. Dependency versions

| dependency | compat bound |
|---|---|
| InfrastructureModels | `~0.6, ~0.7` |
| Interpolations | `0.15` |
| JuMP | `1.23` |
| MathOptInterface | `1.31` |
| PowerModels | `0.21` |
| ProgressMeter | `1.10.0` |
| SignalDecomposition | `1.1` |
| julia | `1.9` |

Test-only (not in `[compat]`): `Clarabel`, `Dierckx`, `Ipopt`, `Test`

---

## O. Implications

The following facts, discovered from reading the actual source, are non-obvious and must be known before writing new code for the `shhc-pce` branch:

1. **No stochastic code exists yet.** The entire `PolyChaos` / PCE / sHHC layer must be written from scratch. There are no stubs, no partial implementations, no relevant imports.

2. **`PolyChaos` is not a dependency.** It appears in neither `[deps]` nor `[compat]`. It must be added to `Project.toml` before any PCE code can run. The test block will also need a `[extras]`/`[targets]` entry if PCE tests are separate.

3. **The multi-network key is `nw`, indexed by harmonic number as a string.** The fundamental is always network `"1"` (Int `1`). Harmonics are `"3"`, `"5"`, etc. The `replicate` function explicitly deletes non-selected networks: `delete!(hdata["nw"], "$nh")`. A PCE extension must preserve this structure and add scenario networks or augment existing ones.

4. **`_HPM.replicate` is the sole data-preparation function.** It is exported as `replicate` and overrides `InfrastructureModels.replicate`. Any stochastic variant must either extend this function (add keyword arguments) or provide a new wrapper that calls it and then augments `hdata`.

5. **The harmonic unit (load) is the only decision variable holder for hosting capacity.** The symbol `:cmd` (load current magnitude) is the quantity being maximized. A stochastic extension must decide whether PCE coefficients are stored per `:cmd`, alongside it, or as a separate variable collection.

6. **No lifted voltage-squared (`w`/`wr`/`wi`) auxiliary variable exists.** The SOC formulation uses `vr`/`vi` directly inside `SecondOrderCone()` constraints. Any new SOC constraint on harmonic voltages must follow the same pattern.

7. **`dHHC_SOC` currently omits all cross-harmonic (RMS, THD) constraints.** The overall `constraint_voltage_rms_limit`, `constraint_voltage_thd_limit`, `constraint_current_rms_limit` (branch), and `constraint_xfmr_current_rms_limit` calls are **commented out** in `build_hhc(pm::dHHC_SOC)`. Only per-harmonic IHD limits are active. A sHHC formulation must decide whether to include these or keep the same omission.

8. **Two distinct `solve_hhc` signatures exist.** The single-optimizer form runs a pure NLP. The two-optimizer form first runs a fundamental HPF (`update_hdata_with_fundamental_hpf_results!`) to fix fundamental voltages and currents, then solves the SOC HHC problem with those values pre-loaded into `hdata`. Any stochastic wrapper must choose which pipeline to extend.

9. **The bus `"vm"` field used in SOC constraint templates is written by `update_hdata_with_fundamental_hpf_results!`**, not from the original data file. It is `floor(..., digits=10)` of the HPF result. The SOC constraint `constraint_voltage_rms_limit(dHHC_SOC)`, etc., reads `_PMs.ref(pm, fundamental(pm), :bus, i, "vm")`. This means the SOC model treats the fundamental as fixed/deterministic even in a putative stochastic extension.

10. **The fairness principle is stored in `pm.data["principle"]` as a plain string** (`"maximum efficiency"`, `"absolute equality"`, `"maximin"`, `"Kalai-Smorodinsky bargaining"`). The variable and constraint dispatch is via `if/elseif` string comparisons, not Julia multiple dispatch. A sHHC extension that needs its own principle must add another `if` branch in `variable_fairness_principle`, `constraint_fairness_principle`, and `objective_maximum_hosting_capacity` in `form/iv.jl`.

11. **The `constraint_load_current` / `constraint_load_current_angle` formulation fixes the angle of harmonic current injection** via `load["ref_angle"]`. In a stochastic context, whether this angle should be a random variable or remain deterministic needs a design decision.

12. **`Clarabel` is the test SOC solver; `Ipopt` is the test NLP solver.** The `dHHC_SOC` path in tests uses `solver_soc = Clarabel.Optimizer`. PCE-based stochastic problems may require a different solver class (e.g., an NLP solver if PCE constraints are polynomial).

13. **`JuMP = "1.23"` is pinned.** `add_nonlinear_constraint` (used in transformer magnetization) is the JuMP 1.x API. Ensure any new PCE constraints use the same JuMP API version.

14. **`MathOptInterface = "1.31"` is pinned.** The `SOCtoNonConvexQuadBridge` bridge is explicitly added in `build_hhc(dHHC_SOC)`. Any new SOC constraints will benefit from this bridge automatically.

15. **The `"load"` component is the harmonic unit key**, not `"generator"` or a custom key. All HHC objective and fairness logic iterates over `_PMs.ids(pm, :load, nw=n)`. Any stochastic extension that models uncertainty in device parameters must operate on `load` data fields.

16. **The branch current limit is derived, not from data directly.** `c_rating` is computed in `replicate` as `rate_a / sqrt(3) / vmmin`. The original MATPOWER `rate_a` field is the source. This derivation happens once at data-preparation time and is stored in the multi-network dict.

17. **`InfrastructureModels` compat allows both `~0.6` and `~0.7`.** This is an unusual double-tilde compat; verify that the PCE library (if it uses `IM` internals) is compatible with both.

18. **The `dhhc_pes_gm` example script is the most complete usage reference.** It covers harmonics 1–50, two network variants, IEEE 519, IEC 61000, and the HHC optimization approach side-by-side. It is not a test and is not run by `runtests.jl`.

---

## REQUIRED FIXES (identified during sHHC-PCE development)

### RF-1: `is_zero_sequence` / `is_pos_sequence` / `is_neg_sequence` break in PCE-expanded multi-network

**Affected functions:** `constraint_xfmr_winding_config`, `constraint_xfmr_winding_current_balance` (both in `src/form/iv.jl`).

**Root cause:** These functions test `is_zero_sequence(n)` where `n` is the JuMP network id passed by the build loop. In the deterministic dHHC multi-network, `n` equals the harmonic number (1, 3, 5, …), so `is_zero_sequence(3) == (3 % 3 == 0) == true`. In the PCE-expanded multi-network produced by `build_mn_pce_data`, network ids are encoded as `nw_id(h_idx, k, P_size)`. For h=3, k=0, P_size=3: `nw_id = 7`. Then `is_zero_sequence(7) == (7 % 3 == 0) == false` — **wrong**.

**Required fix:** Add a helper `harmonic_from_nw(pm, nw)` that reads `ref(pm, nw)["harmonic_idx"]` (stored by `build_mn_pce_data`) and returns the true harmonic number. Replace `is_zero_sequence(n)` calls in these two functions with `is_zero_sequence(harmonic_from_nw(pm, n))`, or add `AbstractSHHCModel`-dispatched overloads.

**Status:** Not yet fixed. The sHHC_SOC build function must not call `constraint_xfmr_winding_config` or `constraint_xfmr_winding_current_balance` until this fix is applied.

### RF-2: `constraint_gen_current` compares `nw::Int` to the function `fundamental` (not its return value)

**Location:** `src/core/constraint_template.jl:143`: `if iszero(inf) && nw ≠ fundamental`

**Root cause:** Missing `(pm)` call — `fundamental` is a function object, so `nw ≠ fundamental` is always `true`. The intended test is `nw ≠ fundamental(pm)`.

**Required fix:** Change to `nw ≠ fundamental(pm)`. One-line fix — safe to apply.
