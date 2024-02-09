<a href="https://github.com/timmyfaraday/HarmonicPowerModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/HarmonicPowerModels.jl/workflows/CI/badge.svg"></img></a>
<a href="https://timmyfaraday.github.io/HarmonicPowerModels.jl/"><img src="https://github.com/timmyfaraday/HarmonicPowerModels.jl/workflows/Documentation/badge.svg"></img></a>

# HarmonicPowerModels.jl

HarmonicPowerModels.jl is an extension package of PowerModels.jl for Steady-State 
Power Network Optimization with Power Harmonics. 

## Core Problem Specification
- Balanced Harmonic Power Flow (hpf)
  - IVR (`IVRPowerModel`)
- Balanced Harmonic Optimal Power Flow (hopf)
  - IVR (`IVRPowerModel`)
- Balanced Harmonic Hosting Capacity (hhc)
  - Deterministic NLP (`dHHC_NLP <: IVRPowerModel`)
  - Deterministic SOC (`dHHC_SOC <: IVRPowerModel`)

## Installation

For now, HarmonicPowerModels is unregistered. Nevertheless, you can install it through

```
] add https://github.com/timmyfaraday/HarmonicPowerModels.jl.git
```

At least one solver is required for running PowerModels.  The open-source solver 
Ipopt is recommended, as it is fast, scaleable and can be used to solve a wide 
variety of the problems and network formulations provided in HarmonicPowerModels. The Ipopt solver can be installed via the package manager with

```julia
] add Ipopt
```

Test that the package works by running

```julia
] test HarmonicPowerModels
```

## Acknowledgements
This code has been developed by the following developers:
  - Tom Van Acker, BASF Antwerp, ([@timmyfaraday](https://github.com/timmyfaraday)), 
  - Hakan Ergun, KU Leuven, ([@hakanergun](https://github.com/hakanergun)), and
  - Frederik Geth, GridQube, ([@frederikgeth](https://github.com/frederikgeth)).

## License
This code is provided under a BSD 3-Clause License.