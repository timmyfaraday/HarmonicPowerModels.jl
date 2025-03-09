# Branches $b \in B$

The branches $b \in B$ are a subset of the edges $e \in E$ of the extended graph, and group the edges for which no voltage transformation occurs, e.g., cables and overhead lines. A branch is uniquely represented by a pi-model.

[Illustration of the equivalent single-phase circuit diagram (pi-model) of a branch](figure/branch_pi_model.JPG)

## Parameters

| name          | symb.                 | unit  | type          | $\subset H$   | def.  | definition                                                            |
|---------------|-----------------------|-------|---------------|---------------|-------|-----------------------------------------------------------------------|
| index         | $b$                   | -     | Int           | 1             | -     | unique index of the branch                                            |
| f_bus         | $i$                   | -     | Int           | 1             | -     | unique index of the from-bus of the branch                            |
| t_bus         | $j$                   | -     | Int           | 1             | -     | unique index of the to-bus of the branch                              |
| r             | $r_{b,h}$             | pu    | Real          | H             | -     | real part of the series impedance of the branch                       |
| x             | $x_{b,h}$             | pu    | Real          | H             | -     | imaginary part of the series impedance of the branch                  |
| g_fr          | $g_{bi,h}$            | pu    | Real          | H             | -     | real part of the shunt admittance at the from-bus of the branch       |
| b_fr          | $b_{bi,h}$            | pu    | Real          | H             | -     | imaginary part of the shunt admittance at the from-bus of the branch  |
| g_to          | $g_{bj,h}$            | pu    | Real          | H             | -     | real part of the shunt admittance at the to-bus of the branch         |
| b_to          | $b_{bj,h}$            | pu    | Real          | H             | -     | imaginary part of the shunt admittance at the to-bus of the branch    |
| i_base_ka     | -                     | kA    | Real          | 1             | -     | base branch current magnitude                                         |
| i_fund_magn   | $I^{fund,magn}_{b}$   | pu    | Vector{Real}  | 1             | [0,0] | fundamental from- and to-side branch current magnitude                | 
| i_rms_max     | $I^{rms,max}_{b}$     | pu    | Real          | 1             | -     | maximum root-mean-square branch current                               |

## Variables

| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                                    |
|---------------|-----------------------|-------|---------------|-----------------------|-----------------------------------------------------------------------------------------------|
| cbr           | $I^{re}_{bij,h}$      | pu    | H             | HarmonicPowerModel    | real part of the harmonic current vector from bus $i$ into branch $b$ for harmonic $h$        |
| cbi           | $I^{im}_{bij,h}$      | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic current vector from bus $i$ into branch $b$ for harmonic $h$   |
| cbsr          | $I^{s,re}_{b,h}$      | pu    | H             | HarmonicPowerModel    | real part of the harmonic series current vector of branch $b$ for harmonic $h$                |
| cbsi          | $I^{s,im}_{b,h}$      | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic series current vector of branch $b$ for harmonic $h$           |

Additional remarks:
- the primal starting value is set to zero for the real and imaginary part of all harmonic branch currents; and
- if bounded, the real and imaginary part of any harmonic branch current is lower (-) and upper (+) bounded to the maximum root-mean-square current limit, i.e., `i_rms_max`.

## Constraints

### HarmonicPowerModel

Branch voltage drop (Ohm's law) - $\forall bij \in T^{b,→}, h \in H$: 
```math
\begin{align}
    U^{re}_{j,h} &= U^{re}_{i,h} - r_{b,h} \cdot I^{s,re}_{bij,h} + x_{b,h} \cdot I^{s,im}_{bij,h} \\ 
    U^{im}_{j,h} &= U^{im}_{i,h} - r_{b,h} \cdot I^{s,im}_{bij,h} - x_{b,h} \cdot I^{s,re}_{bij,h}
\end{align}
```

Branch total current - $\forall bij \in T^{b}, h \in H$:
```math
\begin{align}
    I^{re}_{bij,h} &= I^{re,s}_{bij,h} + g_{bi,h} \cdot U^{re}_{i,h} - b_{bi,h} \cdot U^{im}_{i,h} \\ 
    I^{im}_{bij,h} &= I^{im,s}_{bij,h} + g_{bi,h} \cdot U^{im}_{i,h} + b_{bi,h} \cdot U^{re}_{i,h}
\end{align}
```

Branch root-mean-square current limit - $\forall bij \in T^{b}$:
```math
\begin{align}
    \sum_{h \in H} (I^{re}_{bij,h})^{2} + (I^{im}_{bij,h})^{2} &\leq (I^{rms,max}_{b})^{2}
\end{align}
```

### HarmonicPowerModel

Branch root-mean-square current limit - $\forall bij \in T^{b}$:
```math
\begin{align}
    \sum_{h \in H\backslash\{1\}} (I^{re}_{bij,h})^{2} + (I^{im}_{bij,h})^{2} &\leq (I^{rms,max}_{b})^{2} - (I^{fund,magn}_{b})^2
\end{align}
```