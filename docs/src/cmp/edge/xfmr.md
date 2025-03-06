# Xfmr $x \in X$

The transformers (xfmr) $x \in X$ are a subset of the edges $e \in E$ of the extended graph, and group the edges over which a voltage transformation occurs. An xfmr is uniquely represented by a t-model, split into two parts: core and windings. 

[Illustration of the equivalent single-phase circuit diagram (t-model) of the core of an xfmr](figure/xfrm_core.JPG)
[Illustration of the equivalent single-phase circuit diagram (t-model) of the winding of an xfmr, for the positive and negative sequence components](figure/xfrm_winding_pos_seq.JPG)
[Illustration of the equivalent single-phase circuit diagram (t-model) of the winding of an xfmr, for the zero sequence components](figure/xfrm_winding_zero_seq.JPG)

## Parameters

| name          | symb.                 | unit  | type              | $\subset H$   | def.      | definition                                                            |
|---------------|-----------------------|-------|-------------------|---------------|-----------|-----------------------------------------------------------------------|
| index         | $x$                   | -     | Int               | 1             | -         | unique index of the xfmr                                              |
| bus           | $i$                   | -     | Vector{Int}       | 1             | -         | unique index of the connected buses of the xfmr                       |
| nw            | $N^{w}$               | -     | Int               | 1             | -         | number of windings of the xfmr                                        |
| cnf           | -                     | -     | Vector{String}    | 1             | -         | configuration of the xfmr windings                                    |
| gnd           | -                     | -     | Vector{Bool}      | 1             | -         | grounding of the xfmr windings                                        |
| x             | $x_{x,h}$             | pu    | Real              | H             | -         | core - series reactance of the xfmr                                   |
| b             | $b_{x,h}$             | pu    | Real              | H             | 0         | core - shunt reactance of the xfmr                                    |
| g             | $g_{x,h}$             | pu    | Real              | H             | -         | core - shunt resistance of the xfmr                                   |
| tr            | $t^{re}_{x,h}$        | pu    | Real              | H             | -         | core - real part of the phase shift of the xfmr                       |
| ti            | $t^{im}_{x,h}$        | pu    | Real              | H             | -         | core - imaginary part of the phase shift of the xfmr                  |
| r             | $r_{xi,h}$            | pu    | Vector{Real}      | H             | -         | winding - series resistance of the xfmr                               |
| r_gnd         | $r^{gnd}_{xi,h}$      | pu    | Vector{Real}      | H             | [0,...]   | winding - real part of the grounding impedance of the xfmr            |
| x_gnd         | $x^{gnd}_{xi,h}$      | pu    | Vector{Real}      | H             | [0,...]   | winding - imaginary part of the grounding impedance of the xfmr       |
| i_base_ka     | -                     | kA    | Vector{Real}      | 1             | -         | base xfmr winding current magnitude                                   |
| i_fund_magn   | $I^{fund,magn}_{xi}$  | pu    | Vector{Real}      | 1             | [0,...]   | fundamental xfmr winding current magnitude                            | 
| i_rms_max     | $I^{rms,max}_{xi}$    | pu    | Vector{Real}      | 1             | -         | maximum root-mean-square xfmr winding current                         |

## Variables 

### Voltages

| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                                        |
|---------------|-----------------------|-------|---------------|-----------------------|---------------------------------------------------------------------------------------------------|
| vxr           | $V^{re}_{xij,h}$      | pu    | H             | HarmonicPowerModel    | real part of the harmonic winding voltage vector at bus $i$ of xfmr $x$ for harmonic $h$          |
| vxi           | $V^{im}_{xij,h}$      | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic winding voltage vector at the bus $i$ of xfmr $x$ for harmonic $h$ |
| exr           | $E^{re}_{x,h}$        | pu    | H             | HarmonicPowerModel    | real part of the harmonic excitation voltage at the core of xfmr $x$ for harmonic $h$             |
| exi           | $E^{im}_{x,h}$        | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic excitation voltage at the core of xfmr $x$ for harmonic $h$        |

Additional remarks:
- the primal starting value is set to the relevant bus voltage limit, i.e., `v_rms_max` or `v_ihd_max`, and zero for the real and imaginary part of the harmonic winding voltage, respectively;
- if bounded, the real and imaginary part of the harmonic winding voltage is lower (-) and upper (+) bounded to the relevant bus voltage limit, i.e., `v_rms_max` or `v_ihd_max`;
- the primal starting value is set to the maximum of the relevant voltage limit, i.e., `v_rms_max` or `v_ihd_max`, of all connected buses, and zero for the real and imaginary part of the harmonic excitation voltage, respectively; and
- if bounded, the real and imaginary part of the harmonic excitation voltage is lower (-) and upper (+) bounded to the relevant voltage limit, i.e., `v_rms_max` or `v_ihd_max`, of all connected buses.

### Currents

| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                                        |
|---------------|-----------------------|-------|---------------|-----------------------|---------------------------------------------------------------------------------------------------|
| cxr           | $I^{re}_{xij,h}$      | pu    | H             | HarmonicPowerModel    | real part of the harmonic current vector from bus $i$ into xfmr $x$ for harmonic $h$              |
| cxi           | $I^{im}_{xij,h}$      | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic current vector from bus $i$ into xfmr $x$ for harmonic $h$         |
| cxsr          | $I^{s,re}_{x,h}$      | pu    | H             | HarmonicPowerModel    | real part of the harmonic series current vector of xfmr $x$ for harmonic $h$                      |
| cxsi          | $I^{s,im}_{x,h}$      | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic series current vector of xfmr $x$ for harmonic $h$                 |
| cxmr          | $I^{m,re}_{x,h}$      | pu    | H             | HarmonicPowerModel    | real part of the harmonic magnetizing current vector of xfmr $x$ for harmonic $h$               |
| cxmi          | $I^{m,im}_{x,h}$      | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic magnetizing current vector of xfmr $x$ for harmonic $h$          |

Additional remarks:
- the primal starting value is set to zero for the real and imaginary part of all harmonic xfmr currents; and
- if bounded, the real and imaginary part of any harmonic xmfr current is lower (-) and upper (+) bounded to the maximum root-mean-square current limit, i.e., `i_rms_max`.

## Constraints

### HarmonicPowerModel

Core voltage drop (Ohm's law) - $\forall xij \in T^{x,→}, h \in H$:
```math
\begin{align}
    V^{re}_{xij,h} &= E^{re}_{x,h} - x_{x,h} \cdot I^{s,im}_{x,h} \\
    V^{im}_{xij,h} &= E^{im}_{x,h} + x_{x,h} \cdot I^{s,re}_{x,h}
\end{align}
```

Core current balance (Kirchhoff's current law) - $\forall xij \in T^{x,→}, h \in H$:
```math
\begin{align}
    I^{s,re}_{xji,h} + t^{re}_{x,h} \cdot \Big( I^{s,re}_{xij,h} - I^{m,re}_{x,h} - g_{x,h} \cdot E^{re}_{x,h} \Big) 
                     + t^{im}_{x,h} \cdot \Big( I^{s,im}_{xij,h} - I^{m,im}_{x,h} - g_{x,h} \cdot E^{im}_{x,h} \Big) &= 0 \\
    I^{s,im}_{xji,h} + t^{re}_{x,h} \cdot \Big( I^{s,im}_{xij,h} - I^{m,im}_{x,h} - g_{x,h} \cdot E^{im}_{x,h} \Big) 
                     + t^{im}_{x,h} \cdot \Big( I^{s,re}_{xij,h} - I^{m,re}_{x,h} - g_{x,h} \cdot E^{re}_{x,h} \Big) &= 0
\end{align}
```

Core phase shift - $\forall xji \in T^{x,←}, h \in H$:
```math
\begin{align}
    E^{re}_{x,h} &= t^{re}_{x,h} \cdot V^{re}_{xji,h} - t^{im}_{x,h} \cdot V^{im}_{xji,h} \\
    E^{re}_{x,h} &= t^{re}_{x,h} \cdot V^{im}_{xji,h} + t^{im}_{x,h} \cdot V^{re}_{xji,h}
\end{align}
```

Winding current balance (Kirchhoff's current law) - $\forall xij \in T^{x}, h \in H$:
```math
\begin{align}
    I^{s,re}_{xij,h} &= I^{re}_{xij,h} + xxx
    I^{s,im}_{xij,h} &= I^{im}_{xij,h} + xxx
\end{align}
```

Winding pos./neg. seq. voltage drop (Ohm's law) - $\forall xij \in T^{x}, h \in H^{+} \cup H^{-}$:
```math
\begin{align}
    U^{re}_{i,h} &= V^{re}_{xij,h} + r_{xi,h} \cdot I^{re}_{xij,h} \\
    U^{im}_{i,h} &= V^{im}_{xij,h} + r_{xi,h} \cdot I^{im}_{xij,h}
\end{align}
```

Winding zero seq. voltage drop (Ohm's law) - $\forall xij \in T^{x,earthed}, h \in H^{0}$:
```math
\begin{align}
    U^{re}_{i,h} &= V^{re}_{xij,h} + (r_{xi,h} + 3 \cdot r^{gnd}_{xi,h}) \cdot I^{re}_{xij,h} - 3 \cdot x^{gnd}_{xi,h} \cdot I^{im}_{xij,h} \\
    U^{im}_{i,h} &= V^{im}_{xij,h} + (r_{xi,h} + 3 \cdot r^{gnd}_{xi,h}) \cdot I^{im}_{xij,h} + 3 \cdot x^{gnd}_{xi,h} \cdot I^{re}_{xij,h}
\end{align}
```

Winding zero seq. current blocking - $\forall xij \in T^{x,ne} \cup T^{x,delta}, h \in H^{0}$:
```math 
\begin{align}
    I^{re}_{xij,h} &= 0 \\
    I^{im}_{xij,h} &= 0
\end{align}
```

## Magnetizing Current

### Theoretical Background

In the time domain, the transformer excitation voltage $e(t)$ [V] and magnetizing current $i^m(t)$ [A] are related through the BH-curve:
```math
\begin{align}
    e(t) \rightarrow B(t) \stackrel{\mbox{BH-curve}}{\rightarrow} H(t) \rightarrow i^m(t)
\end{align}
```

[Illustration of the relationship between the excitation voltage~$e(t)$ and magnetizing current~$i^m(t)$](figure/BH_curve.JPG)

The transformer excitation voltage $e(t)$ relates to its magnetic flux density $B(t)$ [T]:
```math
\begin{align}
    B(t)    &= \sum_{h \in H} \frac{|E_h|}{A \cdot \omega_h} \cdot \cos(\omega_h \cdot t + \theta_h), \\
    \mbox{given:} & \nonumber \\
    e(t)    &= \sum_{h \in H} |E_h| \cdot \sin(\omega_h \cdot t + \theta_h), \\
    B(t)    &= \frac{1}{A} \int -e(t) \mathrm{d}t,
\end{align}
```
where $A$ [m^2], $\omega_h$ [rad/Hz], $|E_h|$ [V] and~$\theta_h$ [rad] denotes the core surface, harmonic angular frequency, harmonic excitation voltage magnitude and phase angle, respectively. 

The transformer magnetizing current $i^m(t)$ relates to its magnetic field intensity $H(t)$ [Ampere-turn/meter]:
```math
\begin{align}
    i^m(t)  &= H(t) \cdot l,
\end{align}
```
where $l$ [m] denotes the mean magnetic path. The frequency-domain magnetizing current $I^{e}$ [pu] is determined through a Fourrier transform of the time-domain magnetizing current $i^m(t)$, adjusted for current basis. 

### Implementation

```@docs
    HarmonicPowerModels.sample_magnetizing_current(hdata::Dict{String,<:Any}, xfmr_magn::Dict{String,<:Any})
```

All magnitizing data are stored in a dictionary `xfmr_magn` with:

- General input:
| key       | type          | description                                                       |
|-----------|---------------|-------------------------------------------------------------------|
| Hᴱ        | Vector{Int}   | set of relevant excitation voltage harmonics                      |
| Hᴵ        | Vector{Int}   | set of relevant magnetizing current harmonics                     |
| Eᵐᵃˣ      | Real          | maximum excitation root-mean-square voltage                       |
| IHD       | Vector{Real}  | individual harmonic voltage distortion limits for all $h \in Hᴱ$  |
| pcs       | Vector{Int}   | number of spline pieces for all $h \in Hᴱ$                        |

- Individual xfmr input, i.e, dictionary with key: xfmr id and values:
| key       | type          | description                                                       |
|-----------|---------------|-------------------------------------------------------------------|
| l         | Real          | mean magnetic path [m]                                            |
| A         | Real          | core surface [m^2]                                                |
| N         | Int           | nominal primary turns [-]                                         |
| BH        | Function      | anonymous function for the inversed BH-curve [T//A-turns/m]       |
| Vbase     | Real          | base voltage associated with the primary winding [V]              |
