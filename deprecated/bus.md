# Bus

A bus (`i ∈ 𝓘`) of a circuit represents the vertices of the graph to which 
edges and units are connected. A bus has a number of terminals (`t ∈ 𝓣ᵢ`), 
depending on the physical representation, e.g., balanced `i → t ∈ {1}`, 
four-wire `i → t ∈ {a, b, c, n, g}`, etc.

## Variables

### Voltage Phasor

#### Rectangular Current-Voltage
- real voltage: ure_n,i vector of cardinality `|𝓣ᵢ|`
- imag voltage: uim_n,i vector of cardinality `|𝓣ᵢ|`

## Input

- min RMS voltage
- max RMS voltage
- base RMS voltage

## Constraints

### Power Flow

### Optimal Power Flow
