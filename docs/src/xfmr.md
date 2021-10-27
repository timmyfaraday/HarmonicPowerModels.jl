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
- val = a dictionary consisting of three types of input:
  - General input, including:
| key           | type          | description                                                                       |
|---------------|---------------|-----------------------------------------------------------------------------------|
| E_formulation | Symbol        | excitation voltage formulation, i.e., :rectangular or :polar                      |
| I_formulation | Symbol        | magnetization current formulation, i.e., :rectangular or :polar                   |
| E_harmonics   | Vector{Int}   | set of excitation voltage harmonics                                               |
| I_harmonics   | Vector{Int}   | set of magnetization current harmonics                                            |
  - Magnetization model input, including:
| key           | type          | description                                                                       |
|---------------|---------------|-----------------------------------------------------------------------------------|
| BH_curve      | Function      | anonymous function for the inversed BH-curve [T//A-turns/m]                     |
| mean_path     | Real          | mean magnetic path [m]                                                            |
| core_surface  | Real          | core surface [m^2]                                                                |
  - Excitation voltage input, depending on the chosen E_formulation:
  if E_formulation == :rectangular
| key           | type          | description                                                                       |
|---------------|---------------|-----------------------------------------------------------------------------------|
| dEre          | Vector{Real}  | real excitation voltage step [pu] for all excitation voltage harmonics            |
| Ere_min       | Vector{Real}  | minimum real excitation voltage [pu] for all excitation voltage harmonics         |
| Ere_max       | Vector{Real}  | maximum real excitation voltage [pu] for all excitation voltage harmonics         |
| dEim          | Vector{Real}  | imaginary excitation voltage step [pu] for all excitation voltage harmonics       |
| Eim_min       | Vector{Real}  | minimum imaginary excitation voltage [pu] for all excitation voltage harmonics    |
| Eim_max       | Vector{Real}  | maximum imaginary excitation voltage [pu] for all excitation voltage harmonics    |    
  if E_formulation == :polar
| key           | type          | description                                                                       |
|---------------|---------------|-----------------------------------------------------------------------------------|
| dE            | Vector{Real}  | excitation voltage magnitude step [pu] for all excitation voltage harmonics       |
| dE_min        | Vector{Real}  | minimum excitation voltage magnitude [pu] for all excitation voltage harmonics    |
| dE_max        | Vector{Real}  | maximum excitation voltage magnitude [pu] for all excitation voltage harmonics    |
| dθ            | Vector{Real}  | excitation voltage phase angle step [rad] for all excitation voltage harmonics    |
| dθ_min        | Vector{Real}  | minimum excitation voltage phase angle [rad] for all excitation voltage harmonics |
| dθ_max        | Vector{Real}  | maximum excitation voltage phase angle [rad] for all excitation voltage harmonics |