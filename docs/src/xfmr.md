# Transformer Model

## 

## Transformer Excitation

### Theoretical Background

In the time domain, the transformer excitation voltage~$e(t)$\,[V] and magnetization current~$i^m(t)$\,[A] are related through the BH-curve:
```math
\begin{align}
    e(t) \rightarrow B(t) \stackrel{\mbox{BH-curve}}{\rightarrow} H(t) \rightarrow i^m(t)
\end{align}
```

[Illustration of the relationship between the excitation voltage~$e(t)$ and magnetization current~$i^m(t)$](figure/BH_curve.JPG)

The transformer excitation voltage~$e(t)$ relates to its magnetic flux density~$B(t)$\,[T]:
```math
\begin{align}
    B(t)    &= \sum_{h \in H} \frac{|E_h|}{A \cdot \omega_h} \cdot \cos(\omega_h \cdot t + \theta_h), \\
    \mbox{given:} & \nonumber \\
    e(t)    &= \sum_{h \in H} |E_h| \cdot \sin(\omega_h \cdot t + \theta_h), \\
    B(t)    &= \frac{1}{A} \int -e(t) \mathrm{d}t,
\end{align}
```
where~$A$\,[m^2], $\omega_h$\,[rad/Hz], $|E_h|$\,[V] and~$\theta_h$\,[rad] denotes the core surface, harmonic angular frequency, harmonic excitation voltage magnitude and phase angle, respectively. The link between the frequency-domain excitation current~$E$\,[pu] and the time-domain excitation voltage~$e(t)$, and consequently the magnetic flux density~$B(t)$, is provided through the following functions, depending on the chosen excitation voltage formulation `E_formulation`, i.e., `:polar` or `:rectangular`.
```@docs
HarmonicPowerModels.magnetic_flux_density_polar(E::Vector{<:Real}, θ::Vector{<:Real}, ω::Vector{<:Real}, t::Vector{<:Real}, A::Real, Vbase::Real)
```
```@docs
HarmonicPowerModels.magnetic_flux_density_rectangular(Ere::Vector{<:Real}, Eim::Vector{<:Real}, ω::Vector{<:Real}, t::Vector{<:Real}, A::Real, Vbase::Real)
```

The transformer magnetization current~$i^m(t)$ relates to its magnetic field intensity~$H(t)$\,[Ampere-turn/meter]:
```math
\begin{align}
    i^m(t)  &= H(t) \cdot l,
\end{align}
```
where~$l$\,[m] denotes the mean magnetic path. The frequency-domain magnetization current~$I^{e}$\,[pu] is determined through a Fourrier transform of the time-domain magnetization current~$i^m(t)$, adjusted for current basis. Depending on the chosen excitation voltage formulation `E_formulation`, the frequency-domain magnetization current is expressed in `:polar` or `:rectangular` coordinates.

### Implementation

```@docs
    HarmonicPowerModels.sample_xfmr_excitation(data::Dict{String, <:Any}, xfmr_exc::Dict{Int, Dict{String, <:Any})
```

All excitation data are stored in a dictionary `xfmr_exc` with:
- key = id of the xfmr [`Int`]
- val = a dictionary [`Dict{String,Any}`] consisting of the following input:
  - General input, including:
| key   | type          | description                                                     |
|-------|---------------|-----------------------------------------------------------------|
| "Hᴱ"  | Vector{Int}   | set of relevant excitation voltage harmonics                    |
| "Hᴵ"  | Vector{Int}   | set of relevant magnetizing current harmonics                   |
| "Fᴱ"  | Symbol        | excitation voltage formulation, i.e., :rectangular or :polar    |
| "Fᴵ"  | Symbol        | magnetization current formulation, i.e., :rectangular or :polar |
| "l"   | Real          | mean magnetic path [m]                                          |
| "A"   | Real          | core surface [m^2]                                              |
| "N"   | Int           | nominal primary turns [-]                                       |
| "BH"  | Function      | anonymous function for the inversed BH-curve [T//A-turns/m]     |
