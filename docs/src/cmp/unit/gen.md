# Generators $g \in \mathcal{G}$

Generators $g \in \mathcal{G}$ can be see as shunt impedances to the system which conduct the harmonic currents to the ground. Furthermore, due to their active and reactive power injection, they provide the input for the fundamental frequency power flow which is used as the starting base for the harmonic hosting capacity calculation.

## Parameters  - TODO

| name          | symb.                 | unit  | type      | $\subset H$   | def.  | definition                                                          |
|---------------|-----------------------|-------|-----------|---------------|-------|---------------------------------------------------------------------|
| id            | $g$                   | -     | Int       | 1             | -     | unique index of generator                                           |
| bus           | $i$                   | -     | Int       | 1             | -     | unique index of the bus to which the generator is connected to      |
| bsc           | $b_{g,h}$             | pu    | Real      | H             | -     | short circuit susceptance of the generator                          |
| gsc           | $g_{g,h}$             | pu    | Real      | H             | -     | short circuit conductance of the generator                          |
| i_base_ka     | $I^{base}_{g,h=1}$    | pu    | Real      | 1             | -     | base current magnitude of the generator at fundamental frequency    |
| i_fund_magn   | $I_{g,h=1}$           | pu    | Real      | 1             | -     | fundamental frequency current magnitude of the generator            |
| i_rms_max     | $I^{max}_{g,rms}$     | pu    | Real      | 1             | -     | maximum allowable RMS current  of the generator                     |
| p_fund_max    | $P^{max}_{g,h=1}$     | pu    | Real      | 1             | -     | maximum allowable generator active power at fundamental frequency   |
| p_fund_min    | $P^{min}_{g,h=1}$     | pu    | Real      | 1             | -     | minimum allowable generator active power at fundamental frequency   |
| q_fund_max    | $Q^{max}_{g,h=1}$     | pu    | Real      | 1             | -     | maximum allowable generator reactive power at fundamental frequency |
| q_fund_min    | $Q^{min}_{g,h=1}$     | pu    | Real      | 1             | -     | minimum allowable generator reactive power at fundamental frequency |


## Variables

Includes variables for harmonic source current injections as well as auxillary variables needed to model different fairness criteria
| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                  |
|---------------|-----------------------|-------|---------------|-----------------------|-----------------------------------------------------------------------------|
| chsr          | $I^{re}_{u,h}$        | pu    | H             | HarmonicPowerModel    | real part of the harmonic source current for harmonic source $u$ for harmonic $h$      |
| chsi          | $I^{im}_{u,h}$        | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic source current for lharmonic source $u$  for harmonic $h$ |
| chsm          | $I^{mag}_{u,h}$       | pu    | H             | HarmonicPowerModel    | magnitude of the harmonic source current for harmonic source $u$ for harmonic $h$      |
| cmh          | $I^{mag}_{h}$       | pu    | H             | HarmonicPowerModel    | magnitude of the upper bound harmonic source current for all harmonic sources for harmonic $h$      |
| fh          | $f_{h}$       | -    | H             | HarmonicPowerModel    | common ratio for decreasing harmonic current injection w.r.t. maximum individual limits over all harmonic sources for harmonic $h$      |

Additional remarks:
- the harmonic sources are defined for the set of loads
## Constraints

### HarmonicPowerModel

Linking harmonic current injection and fundamental frequency power injection - $\forall u \in T^{b,U}, h=1\$: 
```math
\begin{align}
P_{u,h} = v^{Re}_{i} \cdot I^{re}_{u,h} + v^{Im}_{i} \cdot I^{Im}_{u,h} \\
Q_{u,h} = v^{Im}_{i} \cdot I^{re}_{u,h} - v^{Re}_{i} \cdot I^{Im}_{u,h}
\end{align}
```

Linking harmonic current angle at fundamental frequency using power factor - $\forall u \in T^{b,U}, h=1\$: 
```math
\begin{align}
I^{re}_{u,h} = I^{mag}_{u,h} \cdot cos(\theta_{u,h}) \\
I^{im}_{u,h} = I^{mag}_{u,h} \cdot cos(\theta_{u,h})
\end{align}
```

Fairness principles for - $\forall u,v \in U, h \in H$:

#### Maximum efficiency principle
```math
\begin{align}
     max \sum_{u \in \mathcal{U}^{h}} {\sum_{h \in \mathcal{H} \backslash \{1\}}} I^{mag}_{u,h} 
\end{align}
```

#### Absolute equality principle
```math
\begin{align}
     max \sum_{u \in \mathcal{U}^{h}}{\sum_{h \in \mathcal{H} \backslash \{1\}}} I^{mag}_{u,h}  \\
     I^{mag}_{u,h} = I^{mag}_{v,h} 
\end{align}
```

#### Maximin principle
```math
\begin{align}
     max \sum_{h \in \mathcal{H} \backslash \{1\}} I^{mag}_{h} \\
     I^{mag}_{h} \leq I^{mag}_{u,h} 
\end{align}
```


#### Kalai-Smorodinsky bargaining
```math
\begin{align}
     max \sum_{h \in \mathcal{H} \backslash \{1\}} f_{h} \\
     I^{mag}_{u,h} = f_{h} \cdot \hat{I}^{max}_{u,h} \\
     0 \leq f_{h} \leq 1 
\end{align}
```