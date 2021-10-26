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

<!-- ### Implementation

This function creates anonymous functions which wrap a spline model of the exitation current. As inputs it takes either the rectangular or polar coordinates of the excitation voltage of transformer x ∈ 𝓧:
    `E^{re}_{h,x}, E^{im}_{h,x}, ∀ h ∈ 𝓗ᵉ,`
        or
    `E_{h,x}, θ_{h,x}, ∀ h ∈ 𝓗ᵉ,`
and outputs either the rectangular or polar coordinates of the exictation 
current of tranformer x ∈ 𝓧:
    `I^{exc,re}_{h,x}, I^{exc,im}_{h,x}, ∀ h ∈ 𝓗ⁱ,`
        or
    `I^{exc}_{h,x}, φ^{exc}_{h,x}, ∀ h ∈ 𝓗ⁱ,`
where 𝓗ᵉ and 𝓗ⁱ denote the set of excitation voltage and current harmonics, 
respectively. -->