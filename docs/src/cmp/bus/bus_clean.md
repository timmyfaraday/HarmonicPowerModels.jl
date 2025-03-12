# Clean buses $i \in I^{cln}$

The clean buses $i \in I^{cln}$ are the nodes of the extended graph to which a clean reference frame is applied. To this end, the bus harmonic voltage of a clean bus is fixed by setting the real and imaginary part of the harmonic voltage vector to zero. 

## Parameters

| name          | symb.                 | unit  | type      | $\subset H$   | def.  | definition                                            |
|---------------|-----------------------|-------|-----------|---------------|-------|-------------------------------------------------------|
| v_fund_ref    | $U^{fund,ref}_{i}$    | pu    | Real      | 1             | 1.0   | real part of the fundamental bus voltage reference    |

## Constraints

### HarmonicPowerModel

Clean bus real harmonic voltage - $\forall i \in I^{cln}$:
```math
\begin{align}
    U^{re}_{i,1} &= U^{fund,ref}_{i} 
    U^{im}_{i,1} &= 0 
\end{align}
```

Clean bus harmonic voltage - $\forall i \in I^{cln}, h \in H$:
```math
\begin{align}
    U^{re}_{i,h} &= 0 
    U^{im}_{i,h} &= 0 
\end{align}
```