<a href="https://github.com/timmyfaraday/HarmonicPowerModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/HarmonicPowerModels.jl/workflows/CI/badge.svg"></img></a>
<a href="https://timmyfaraday.github.io/HarmonicPowerModels.jl/"><img src="https://github.com/timmyfaraday/HarmonicPowerModels.jl/workflows/Documentation/badge.svg"></img></a>

# HarmonicPowerModels.jl

HarmonicPowerModels.jl is an extension package of PowerModels.jl for Steady-State 
Power System Optimization with Power Harmonics. 

## Core Problem Specification
- Balanced Harmonic Power Flow (hpf)
  - IVR (`IVRPowerModel`)
- Balanced Harmonic Optimal Power Flow (hopf)
  - IVR (`IVRPowerModel`)
- Balanced Harmonic Hosting Capacity (hhc)
  - Deterministic NLP (`dHHC_NLP <: IVRPowerModel`)
  - Deterministic SOC (`dHHC_SOC <: IVRPowerModel`)

## Installation

The package requires `Julia 1.9` or newer. The latest stable release of `HarmonicPowerModels` can be installed using the Julia package manager with

```julia
] add HarmonicPowerModels
```

Test that the package works by running

```julia
] test HarmonicPowerModels
```

## Examples

The numerical illustrations and case studies of the papers written as part of the package development can be found in the `examples`-folder.

## Acknowledgements
The primary developer is Tom Van Acker, BASF Antwerp, ([@timmyfaraday](https://github.com/timmyfaraday)), with support from the following contributors: 
  - Hakan Ergun, KU Leuven, ([@hakanergun](https://github.com/hakanergun)), and
  - Frederik Geth, GridQube, ([@frederikgeth](https://github.com/frederikgeth)).

## License
This code is provided under a BSD 3-Clause License.