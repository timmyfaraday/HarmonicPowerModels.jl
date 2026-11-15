# Saved reduction evidence

These tables retain the accepted `ieee30-full-exact` derivation at characteristic
harmonics 5 through 25, with exact moments and the Gaussian factor. The run used
HPM revision `4aec30b05762160dac375a747807c2cf079b6a19` with the uncommitted level
port and LCM revision `e180b435e39495eeb7eb16dbdaa5193d374a9f60`.

The CSV files are unchanged copies of the numerical search tables written by
`reduce.jl`. The hull and cell tables record the bounds and subdivisions.
The other tables record admission, coordinate sharing, source dependencies,
joint extremes and equivalent containment results.
`factor-verification.toml` retains the full factor containment and screen-loss
results from that run's generated surrogate data. Its empty escaped-index lists
correspond to the escape count of 0. The script now writes this table with
each full derivation.

Default tests read these inputs to the generators. They compare the generated
cases, data and reduction reports with the stored deliverables byte for byte.
Separate model evaluations check containment at 256 fixed Sobol points per design.
The full searches and 4096-point checks remain available through
`RHHC_REDERIVE=1` or `reduce.jl`. See
[reproduction and tests](../../../REDUCTION.md#reproduction-and-tests).
