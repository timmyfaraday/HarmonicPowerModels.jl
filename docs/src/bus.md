# Buses $i \in I$

The buses $i \in I$ are the nodes of the extended graph. In general, a harmonic voltage is defined at each bus and limits are enforced on the (harmonic) bus voltage.

## Parameters

| name          | symb.                 | unit  | type      | $\subset H$   | def.  | definition                                            |
|---------------|-----------------------|-------|-----------|---------------|-------|-------------------------------------------------------|
| id            | $i$                   | -     | Int       | 1             | -     | unique index of the bus                               |
| std           | -                     | -     | String    | 1             | -     | relevant standard for voltage quality                 |
| type          | -                     | -     | Int       | 1             | -     | bus type, where ...                                   |
| v_base_kv     | -                     | kV    | Real      | 1             | -     | base bus voltage magnitude                            |
| v_fund_magn   | $U^{fund,magn}_{i}$   | pu    | Real      | 1             | 1.0   | fundamental bus voltage magnitude                     |
| v_ihd_max     | $U^{ihd,max}_{i,h}$   | pu    | Real      | H/1           | -     | maximum individual harmonic bus voltage distortion    |
| v_rms_min     | $U^{rms,min}_{i}$     | pu    | Real      | 1             | -     | minimum root-mean-square bus voltage                  |
| v_rms_max     | $U^{rms,max}_{i}$     | pu    | Real      | 1             | -     | maximum root-mean-square bus voltage                  |
| v_thd_max     | $U^{thd,max}_{i}$     | pu    | Real      | 1             | -     | maximum total harmonic bus voltage distortion         |

## Variables

| name          | symb.                 | unit  | $\subset H$   | formulation           | definition                                                                |
|---------------|-----------------------|-------|---------------|-----------------------|---------------------------------------------------------------------------|
| vr            | $U^{re}_{i,h}$        | pu    | H             | HarmonicPowerModel    | real part of the harmonic voltage vector at bus $i$ for harmonic $h$      |
| vi            | $U^{im}_{i,h}$        | pu    | H             | HarmonicPowerModel    | imaginary part of the harmonic voltage vector at bus $i$ for harmonic $h$ |

Additional remarks:
- the primal starting value is set to the relevant voltage limit, i.e., `v_rms_max` or `v_ihd_max`, and zero for the real and imaginary part of the harmonic bus voltage, respectively; and
- if bounded, the real and imaginary part of the harmonic bus voltage is lower (-) and upper (+) bounded to the relevant voltage limit, i.e., `v_rms_max` or `v_ihd_max`.

TODO add figure 

## Constraints

### HarmonicPowerModel

xxx

### dHHCPowerModel

xxx