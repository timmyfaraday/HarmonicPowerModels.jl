# HarmonicPowerModels.jl

## HarmonicHostingCapacity

### bus
- id::Int [1]
- type::Int [1] in {1,2,3}, where

- Pd::Real [1] - Fundamental active power
- Qd::Real [1] - Fundamental reactive power

- baseKV::Real [1] - voltage

- vre_ref::Real [h] if type == 3
- vim_ref::Real [h] if type == 3

- standard [1]

OR, alternatively -- preferably:

- vmin_rms::Real [1]
- vmax_rms::Real [1]
- vmax_ihd::Real [h/1]
- vmax_thd::Real [1]

### branch

- f_bus::Int [1]
- t_bus::Int [1]

- r::Real [h]
- x::Real [h]
- b::Real [h]

- cmax_rms [1]

### xfmr

- f_bus::Int [1]
- t_bus::Int [1]

- xsc::Real [h]
- gsh::Real [h]
- r1::Real [h]
- r2::Real [h]

- vg::String [1]
- gnd::Array{Bool} [1]
- re::Array{Real} [h]
- xe::Array{Real} [h]

- cmax_rms [1]

### generator

- Pmin::Real [1]
- Pmax::Real [1]
- Qmin::Real [1]
- Qmax::Real [1]

- rsc::Real [h]
- xsc::Real [h]

- cmax_rms [1]

### load

- cmax_rms [1]