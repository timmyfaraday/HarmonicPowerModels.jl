# Harmonic sources $r \in \mathcal{R}$

The harmonic sources $r \in \mathcal{R}$ are the harmonic "pollution" sources used for harmonic hosting capacity calculation. They inject harmonic currents at the specified frequencies into the system.

## Parameters  - TODO

| name          | symb.                 | unit  | type      | $\subset H$   | def.  | definition                                                           |
|---------------|-----------------------|-------|-----------|---------------|-------|----------------------------------------------------------------------|
| index         | $r$                   | -     | Int       | 1             | -     | unique index of the harmonic source                                  |
| bus           | $i$                   | -     | Int       | 1             | -     | unique index of the bus to which the harmonic source is connected to |
| p_fund        | $P_{r,h=1}$           | pu    | Real      | 1             | -     | active power injection of the source at fundamental frequency        |
| q_fund        | $Q_{r,h=1}$           | pu    | Real      | 1             | -     | reactive power injection of the source at fundamental frequency      |
| i_base_ka     | $I_{r,h=1}$           | pu    | Real      | 1             | -     | base current magnitude of the source at fundamental frequency        |
| crar         | $\theta_{r,h=1}$      | rad   | Real      | 1             | -     | current angle reference of the source at fundamental frequency       |
| crhmax         | $\hat{I}^{max}_{r,h}$      | pu   | Real      | H             | -     | current angle reference of the source at fundamental frequency       |


## Variables

Includes variables for harmonic source current injections as well as auxillary variables needed to model different fairness criteria
| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                  |
|---------------|-----------------------|-------|---------------|-----------------------|-----------------------------------------------------------------------------|
| crr          | $I^{re}_{r,h}$        | pu    | H             | HarmonicPowerModel    | real part of the harmonic source current for harmonic source $u$ for harmonic $h$      |
| cri          | $I^{im}_{r,h}$        | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic source current for lharmonic source $u$  for harmonic $h$ |
| crm          | $I^{mag}_{r,h}$       | pu    | H             | HarmonicPowerModel    | magnitude of the harmonic source current for harmonic source $u$ for harmonic $h$      |
| cmh          | $I^{mag}_{h}$       | pu    | H             | HarmonicPowerModel    | magnitude of the upper bound harmonic source current for all harmonic sources for harmonic $h$      |
| fh          | $f_{h}$       | -    | H             | HarmonicPowerModel    | common ratio for decreasing harmonic current injection w.r.t. maximum individual limits over all harmonic sources for harmonic $h$      |

Additional remarks:
- the harmonic sources are defined for the set of loads
## Constraints

### HarmonicPowerModel

Linking harmonic current injection and fundamental frequency power injection - $\forall r \in R, h=1\$: 
```math
\begin{align}
P_{r,h} = v^{Re}_{i} \cdot I^{re}_{r,h} + v^{Im}_{i} \cdot I^{Im}_{r,h} \\
Q_{r,h} = v^{Im}_{i} \cdot I^{re}_{r,h} - v^{Re}_{i} \cdot I^{Im}_{r,h}
\end{align}
```

Linking harmonic current angle at fundamental frequency using power factor - $\forall r \in R, h=1\$: 
```math
\begin{align}
I^{re}_{r,h} = I^{mag}_{r,h} \cdot cos(\theta_{r,h}) \\
I^{im}_{r,h} = I^{mag}_{r,h} \cdot cos(\theta_{r,h})
\end{align}
```

Fairness principles for - $\forall r,t \in U, h \in H$:

#### Maximum efficiency principle
```math
\begin{align}
     max \sum_{r \in \mathcal{R}^{h}} {\sum_{h \in \mathcal{H} \backslash \{1\}}} I^{mag}_{r,h} 
\end{align}
```

#### Absolute equality principle
```math
\begin{align}
     max \sum_{r \in \mathcal{R}^{h}}{\sum_{h \in \mathcal{H} \backslash \{1\}}} I^{mag}_{r,h}  \\
     I^{mag}_{r,h} = I^{mag}_{t,h} 
\end{align}
```

#### Maximin principle
```math
\begin{align}
     max \sum_{h \in \mathcal{H} \backslash \{1\}} I^{mag}_{h} \\
     I^{mag}_{h} \leq I^{mag}_{r,h} 
\end{align}
```


#### Kalai-Smorodinsky bargaining
```math
\begin{align}
     max \sum_{h \in \mathcal{H} \backslash \{1\}} f_{h} \\
     I^{mag}_{r,h} = f_{h} \cdot \hat{I}^{max}_{r,h} \\
     0 \leq f_{h} \leq 1 
\end{align}
```