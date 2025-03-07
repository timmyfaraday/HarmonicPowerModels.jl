# Buses $i \in I$

The buses $i \in I$ are the nodes of the extended graph. In general, a harmonic voltage is defined at each bus and limits are enforced on the (harmonic) bus voltage.

## Parameters  - TODO

| name          | symb.                 | unit  | type      | $\subset H$   | def.  | definition                                            |
|---------------|-----------------------|-------|-----------|---------------|-------|-------------------------------------------------------|
| index         | $i$                   | -     | Int       | 1             | -     | unique index of the bus                               |
| std           | -                     | -     | String    | 1             | -     | relevant standard for voltage quality                 |
| type          | -                     | -     | Int       | 1             | -     | bus type, see matpower manuel                         |
| v_base_kv     | -                     | kV    | Real      | 1             | -     | base bus voltage magnitude                            |
| v_fund_magn   | $U^{fund,magn}_{i}$   | pu    | Real      | 1             | 1.0   | fundamental bus voltage magnitude                     |
| v_ihd_max     | $U^{ihd,max}_{i,h}$   | pu    | Real      | H/1           | -     | maximum individual harmonic bus voltage distortion    |
| v_rms_min     | $U^{rms,min}_{i}$     | pu    | Real      | 1             | -     | minimum root-mean-square bus voltage                  |
| v_rms_max     | $U^{rms,max}_{i}$     | pu    | Real      | 1             | -     | maximum root-mean-square bus voltage                  |
| v_thd_max     | $U^{thd,max}_{i}$     | pu    | Real      | 1             | -     | maximum total harmonic bus voltage distortion         |

## Variables

| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                  |
|---------------|-----------------------|-------|---------------|-----------------------|-----------------------------------------------------------------------------|
| csr            | $I^{re}_{i,h}$       | pu    | H             | HarmonicPowerModel    | real part of the harmonic source current for load $l$ for harmonic $h$      |
| csi            | $I^{im}_{i,h}$       | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic source current for load $l$ for harmonic $h$ |
| csm            | $I^{mag}_{i,h}$      | pu    | H             | HarmonicPowerModel    | magnitude of the harmonic source current for load $l$ for harmonic $h$      |

Additional remarks:
- the harmonic sources are defined fot eh set of loads
## Constraints

### HarmonicPowerModel

Fairness principles for - $\forall l,m \in L, h \in H$:

#### Absolute equality principle
```math
\begin{align}
     I^{mag}_{l,h} = I^{mag}_{m,h} 
\end{align}
```