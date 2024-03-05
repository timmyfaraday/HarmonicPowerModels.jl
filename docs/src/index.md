# HarmonicPowerModels.jl Documentation

```@meta
CurrentModule = HarmonicPowerModels
```

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

At least one solver is required for running HarmonicPowerModels.  The open-source solver Ipopt is recommended, as it is fast, scaleable and can be used to solve a wide variety of the problems provided in HarmonicPowerModels. The Ipopt solver can be installed via the package manager with

```julia
] add Ipopt
```

Test that the package works by running

```julia
] test HarmonicPowerModels
```

## Acknowledgements
The primary developer is Tom Van Acker, BASF Antwerp, ([@timmyfaraday](https://github.com/timmyfaraday)), with support from the following contributors: 
  - Hakan Ergun, KU Leuven, ([@hakanergun](https://github.com/hakanergun)), and
  - Frederik Geth, GridQube, ([@frederikgeth](https://github.com/frederikgeth)).

## License
This code is provided under a BSD 3-Clause License.