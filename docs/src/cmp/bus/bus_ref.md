# Reference buses $i \in I^{ref}$

The reference buses $i \in I^{ref}$ are the nodes of the extended graph to which a reference frame is applied. To this end, the bus voltage angle of a reference bus is fixed by setting the imaginary part of the harmonic voltage vector to zero.

## Parameters

| name          | symb.                 | unit  | type      | $\subset H$   | def.  | definition                                            |
|---------------|-----------------------|-------|-----------|---------------|-------|-------------------------------------------------------|
| v_fund_ref    | $U^{fund,ref}_{i}$    | pu    | Real      | 1             | 1.0   | real part of the fundamental bus voltage reference    |

## Constraints

### HarmonicPowerModel

Reference bus real fundamental voltage - $\forall i \in I^{ref}$:
```math
\begin{align}
    U^{re}_{i,1} &= U^{fund,ref}_{i} 
\end{align}
```

Reference bus imaginary voltage - $\forall i \in I^{ref}, h \in H$:
```math
\begin{align}
    U^{im}_{i,h} &= 0 
\end{align}
```