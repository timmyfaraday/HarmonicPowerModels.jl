# Physical-input register

This register lists the physical inputs of the IEEE 30-bus cables and overhead lines. Each
input has a nominal value with its source and an interval with its basis. The studies of
`examples/rhhc-variants/` cite its rows by their `id`. An audit checked the register against
licensed copies of the standards it cites. A value enters a spec only after the register holds
it.

| File | Content |
| --- | --- |
| `inputs.csv` | one row per input, with the sixteen columns below |
| `derive.py` | the arithmetic behind every derived value, in the Python standard library |
| `derived.json` | the output of `derive.py` |

Run `python3 -I examples/rhhc-variants/data/register/derive.py` from the repository root to
reproduce `derived.json`.

## Columns

| Column | Content |
| --- | --- |
| `id` | the row's name, which the specs cite |
| `construction` | the spec or the family the row describes |
| `component` | the part of the construction |
| `quantity` | the physical quantity |
| `lcm_input` | the LineCableModels input it sets, or its role as a check value |
| `unit` | the unit of `nominal`, `lower` and `upper` |
| `nominal` | the nominal value |
| `definition_class` | the class of the nominal value |
| `definition_source` | the source of the nominal value |
| `lower`, `upper` | the interval |
| `interval_class` | the class of the interval |
| `interval_basis` | the basis of the interval |
| `sensitivity` | the effect of the input on the harmonic network |
| `coordinate` | `yes` for an established coordinate, `sweep` for an interval awaiting the sensitivity study, `no` for a fixed input, an absorbed quantity or a value outside LCM. Rows outside the full IEEE30 audit retain their earlier labels |
| `notes` | remarks and the history of the value |

## Classes

| Nominal value | Definition |
| --- | --- |
| STD | a value that a standard fixes, such as a tabulated thickness or a maximum resistance |
| CAT | a value in the datasheet of a commercial product, with its page |
| DER | a value that `derive.py` computes from cited input rows, including stated assumptions |
| ASM | an engineering assumption, with its basis |

| Interval | Definition |
| --- | --- |
| N | the acceptance limit of a standard, often one-sided |
| P | the tolerance that the manufacturer prints for the product |
| E | production evidence, the values that type-test reports measured |
| M | the range of a material family, across compounds or products |
| S | a study interval inside the N, P and E limits, with its reason |
| ENV | installation, environment or operation, never a manufacturing tolerance |
| FIX | a discrete choice of construction, or an input that the result does not depend on |
| DER | an interval that follows from other rows |

An asterisk, as in `N*`, `P*` or `CAT*`, marks a value taken by analogy from another product or
standard. It also marks a clause that the audit did not check against a licensed copy.

The candidate study is the physical-candidate study of the IEEE 30-bus case that preceded this
register. It is not part of the repository. A row that names it records a choice or a value of
that study, and the row's class says how far the audit supports it.

## Catalogue values used as checks

The studies compare their U33 specs with these values of the Olex catalogue [OLEX, printed pp.
12, 13 and 44]. Resistance and reactance are in Ω/km, capacitance in µF/km and the rating in A.
The values hold for three cables in contact in trefoil, with the screens bonded at both ends.

| Spec | R at 90 °C | X at 50 Hz | C | Diameter over insulation mm | Overall diameter mm | Buried rating |
| --- | --- | --- | --- | --- | --- | --- |
| U33-400 | 0.101 | 0.116 | 0.267 | 40.7 | 50.9 | 487 |
| U33-500 | 0.0797 | 0.111 | 0.297 | 44.2 | 54.6 | 554 |
| U33-630 | 0.0629 | 0.108 | 0.329 | 48.0 | 58.6 | 625 |
| U33-800 | 0.0507 | 0.104 | 0.366 | 52.3 | 63.4 | 699 |

The ratings hold for a soil at 25 °C with a thermal resistivity of 1.2 K·m/W. The center of the
group lies 0.8 m deep [OLEX, p. 44].

## IEEE 30-bus assignments

Each 33 kV branch takes the smallest U33 size whose buried rating is at least the branch current
at its `rateA` at 33 kV. The rule gives U33-400 on branches 20, 31, 33 and 34, U33-630 on branch
26, U33-800 on branch 21, and U33-500 on every other cable.

| Branch | Buses | rateA MVA | Current A | Spec |
| --- | --- | --- | --- | --- |
| 17, 18, 22, 23, 24, 28, 29, 30, 32 | as the case | 29 | 507 | U33-500 |
| 19, 25, 27 | as the case | 30 | 525 | U33-500 |
| 35, 37, 38, 39 | as the case | 28 | 490 | U33-500 |
| 20 | 14–15 | 20 | 350 | U33-400 |
| 31 | 22–24 | 26 | 455 | U33-400 |
| 33 | 24–25 | 27 | 472 | U33-400 |
| 34 | 25–26 | 25 | 437 | U33-400 |
| 26 | 10–17 | 33 | 577 | U33-630 |
| 21 | 16–17 | 38 | 665 | U33-800 |

The 132 kV lines take the ACSR codes of the candidate study. That study rated each code with the
steady heat balance of IEEE 738. Among the codes whose rating meets the branch current, it took
the one whose resistance matches the case best. The audit checked the code data against the
table of the manufacturer, not the selection.

| Branch | Buses | Code |
| --- | --- | --- |
| 1, 4, 6 | 1–2, 3–4, 2–6 | Lion |
| 2, 5, 41 | 1–3, 2–5, 6–28 | Goat |
| 3, 7, 9, 10, 40 | 2–4, 4–6, 6–7, 6–8, 8–28 | Bear |
| 8 | 5–7 | Panther |

## Hello-world assignments

Branches 2 and 3 use U33-630 in touching trefoil. Branch 4 uses Goat on O33-H.
The routes come from `examples/rhhc/data/hello_world.m:38-40`.

| Branch | Buses | Spec | Formation | Length m |
| --- | --- | --- | --- | --- |
| 2 | 1–2 | U33-630 | trefoil | 15000 |
| 3 | 2–3 | U33-630 | trefoil | 10000 |
| 4 | 1–4 | Goat | O33-H | 10000 |

## References

The audit read the licensed standards in the copies of the KU Leuven license.

| Key | Standard | Parts cited |
| --- | --- | --- |
| IEC60228 | NBN EN IEC 60228:2024 (IEC 60228:2023, Edition 4.0) | 3.2, 6.3, 6.4, Table 4, Annex B, Table C.2 |
| IEC60502-2 | NBN IEC 60502-2:2023 (IEC 60502-2:2014) | 4.2, 14.3, 17.5.2, 17.5.3, 17.6, 18.2.6, 18.2.10.3, Tables 3, 6 and 15 |
| IEC60840 | IEC 60840:2011, Edition 4.0 | 10.5, 10.6.2, 10.6.3, 10.7, 10.10, 12.4.5, 12.4.9.2, Table 3 |
| IEC62067 | IEC 62067:2011, Edition 2.0 | 10.6.2, 10.6.3, 10.7.1, 10.10 |
| IEC60287-1-1 | NBN IEC 60287-1-1:2023 (IEC 60287-1-1:2006, Edition 2.0) | 2.3.1, Tables 1 and 3 |
| IEC62641 | NBN EN IEC 62641:2022 | 6.4.2, 6.4.7, Tables 1 and 2 |
| IEC63248 | NBN EN IEC 63248:2022 | Tables A.4 and B.1 |
| IEC60865-1 | NBN EN 60865-1:2015 (IEC 60865-1:2011, Edition 3.0) | Clause 1, Annex A.11 |
| EN50182 | BS EN 50182:2001 | 5.4, 5.5.2, 5.5.4 to 5.5.6, 5.7.3, 5.10, Tables 1, 2 and 4, Annex F Table F.3 |
| HD620 | DIN VDE 0276-620:2009-05 (HD 620 S1:1996/A3:2007) | Parts 1 and 5-C |
| NGED | National Grid Electricity Distribution EE SPEC 77/5 (2024) | 5.0, 6.0, 7.0, 10.0, 11.0 |

Table 4 of IEC 60228:2023 holds the values of Table 2 of IEC 60228:2004. IEC 60287-1-1 has a
current Edition 3.0 of 2023. The licensed set lacks IEC 60028, IEC 60949 and IEC 60986. It also
lacks AS/NZS 1429.1. Table 1 of IEC 60287-1-1 gives the copper value of IEC 60028,
1.7241 × 10⁻⁸ Ω·m.

| Key | Paper | Parts cited |
| --- | --- | --- |
| TB531 | CIGRE WG B1.30, Technical Brochure 531 (2013) | 4.2.2.3 and its Table 6, 4.2.2.5, 5.3, 5.4 |
| GUSTAVSEN2005 | Gustavsen, Martinez and Durbak, IEEE Trans. Power Delivery 20(3), 2005, doi:10.1109/TPWRD.2005.848774 | III-A |
| GUSTAVSEN2016 | Gustavsen, Høyer-Hansen, Triverio and Patel, IEEE Trans. Power Delivery 31(6), 2016 | conclusion 4 |
| KASAS2015 | Kasas-Lazetic and coauthors, Acta Electrotehnica 56(3), 2015 | p. 3 |
| MEYBERG2021 | Meyberg and coauthors, IET Science, Measurement and Technology 15, 2021, doi:10.1049/smt2.12016 | pp. 143 to 153 |
| KING2011 | King, PhD thesis, Cardiff University, 2011 | chapter 4 |

| Key | Document | Link |
| --- | --- | --- |
| OLEX | Nexans Olex High Voltage catalogue OLC12641, printed pp. 11 to 13 and 44 | <https://www.nexans.com.au/dam/jcr:c332abc1-1637-4cbd-9fae-6a65a5485d97/OLC12641_HighVoltageCat.pdf> |
| TFK | TF Kable catalogue of high and extra high voltage cables, printed p. 40 | <https://www.tfkabely.cz/wp-content/uploads/High-and-extra-high-voltage-cables.pdf> |
| TREFINASA | Trefilados de Navarra overhead conductors (2016), p. 6 British-size masses. French-size 30/7 mechanical table. Material table for galvanized steel | <https://www.trefinasa.com/wp-content/uploads/2017/04/overhead_conductors_trefinasa_2016web.pdf> |
| APAR | APAR Industries ACSR table to BS 215 Part 2 | <https://apar.com/wp-content/uploads/2021/02/1.-ACSR.pdf> |
| KEMA | KEMA Laboratories type test certificate 1533-16 (2016) | <http://nci.ae/ttc/5.%20KEMA%201533-16%20%20(3x300%20CU-XLPE-LAT-SWA-PE%2011%20kV).pdf> |
| KIN | Kinectrics type test report K-422899-RC-0001-R00 (2013) | <https://www.nci.ae/ttc/28.%20TT%20REPORT%203X300%20MM2%2015KV%20KINECTRICS%20CANADA%20-%202013.pdf> |
| DEKRA | DEKRA test report 2236358.51-CAT (2019) | <https://www.nci.ae/ttc/60.%20DEKRA-TYPE-TEST-REPORT-IEC-60502-1-4x150-CU-XLPE-SWA-MDPE.pdf> |
| ENEL | Enel Global Standard GSC003 Rev. 2 (2016) | <https://www.eneldistribuicao.com.br/ce/documentos/GSC-003_R-02.pdf> |
| ESKOM | Eskom 240-185000201 Rev. 1, section 3.3.1 and Appendix D Schedule A2 | <https://www.etenders.gov.za/home/Download/?blobName=a435f872-aa2f-4993-a7d3-91c2c6e0678a.pdf&downloadedFileName=240-185000201+-+Rev.1+%281%29+Published.pdf> |
| BOREALIS | Borealis data sheets Borlink LS4201R and LE0595 (2024), typical values | the manufacturer's site |
| HFDK-4201 | Dow ENDURANCE HFDK-4201 SC data sheet, typical values | <https://www.dow.com/en-us/pdp.endurance-hfdk-4201-sc-compound-for-cable-systems.127277z.html> |
| NKT | a third-party copy of the NKT 132 kV brochure, which the audit did not check | <https://manuals.plus/m/72406bd69b69737a8f8c0f48714be3adea8787f4cee0890238bc0c237f976a2b.pdf> |
| POWERGRID | Power Grid Corporation of India technical specification, Volume II Part II | <https://apps.powergrid.in/pgciltenders/Files/034a2eab-3098-45ba-9667-3d6f53cc9f4b/Volume-II_Part-II_TS-1.pdf> |
| ETHEKWINI | eThekwini Electricity, Rev. 04 (2022-08-30), 4.5.2, 4.5.3, 4.6.1, 4.6.2, 4.7.2 and 8.7.8.14 | <https://www.etenders.gov.za/home/Download/?blobName=d8f1d52e-319c-4122-8f67-bd084afe9e91.pdf&downloadedFileName=3-GenTechSpec-Construct-Steel-Lattice+Power+Lines.pdf> |
| UKPN | UK Power Networks EDS 08-4000 | <https://media.umbraco.io/uk-power-networks/iunhkfsn/eds-08-4000-ehv-network-design.pdf> |
| SPEN | SP Energy Networks OHL-03-099 | <https://www.spenergynetworks.co.uk/userfiles/file/OHL-03-099.pdf> |
| AUSGRID | Ausgrid NS116, Annex C | <https://www.ausgrid.com.au/-/media/Documents/Technical-Documentation/NS/NS116.pdf> |
| COPPERFYI | European Copper Institute note on copper conductivity | <https://help.copper.fyi/hc/en-us/articles/4410229175442-How-does-copper-set-the-standard-for-electrical-conductivity> |
| CALEDONIAN | Caledonian Cables table of thermoplastic properties | <https://www.caledonian-cables.net/Material-Properties/Thermoplastrc.html> |
| SUDA | Wuxi Suda copolymer-coated aluminium tape | <https://www.copolymercoatedaluminumtape.com/copolymer-coated-tape/copolymer-coated-aluminum-tape/plastic-coated-aluminum-tape.html> |
| VIGEPD | Vig Metal EPD-IES-0028115:003 (2026) | <https://www.environdec.com/library/epd28115> |
| IEEE738 | IEEE 738-2023, 4.4 and 4.6.1.1, heat balance and illustrative rating weather | local corpus copy |
| ASTM B3 | the ASTM standard for annealed copper wire, by analogy only | not read |

## Specs pending audit

The alternative-system candidate study defines the designs below. They await audit and
have no rows in `inputs.csv`. The pipeline cases retain their source electrical data with
`ScaledBranchModel`. Industrial and Antwerp vary R, L and C. IEEE33, IEEE9, IEEE39
and the stars vary R and L, with capacitance fixed. The ranges come from
`examples/rhhc/ieee30_robust.jl:58-60`. These ranges describe network parameters. They are
not manufacturing tolerances. Source equivalents, transformers and reactors retain their
source data. The stars retain their identical branches.

The definitions below come from an unaudited candidate study that is not part of the
repository. The table gives the conductor and cable constructions. The following paragraph
gives the structure geometry. These definitions await audit.

| System | Cable specs pending audit | Overhead specs pending audit |
| --- | --- | --- |
| IEEE33 | U22-Al95, U22-Al150, U22-Al240, U22-Al300, U22-Al400, U22-Al500, U22-Al630 and U22-Al1200: 12.7/22 kV aluminum/XLPE, copper wire screens, PVC jackets, touching trefoil at 0.8 m depth | Rabbit, Dog, Bison, Weasel, Tiger and Wolf ACSR on O13-F |
| Industrial | U36-Al400: 20.8/36 kV aluminum/XLPE, 35 mm² copper screen, PE jacket. U11-Al240: 6.35/11 kV aluminum/XLPE, copper screen, PVC jacket. Both in touching trefoil at 0.8 m depth. NYY-Cu240-120: three 240 mm² copper phases and one 120 mm² core, PVC insulation and jacket | O36-F with Goat and O10-F with Panther. The conductors are audited, the structures await audit. The elevated NYY installation at 6 m also awaits audit |
| Antwerp | U155-Cu1600: 87/150–161 kV, five-sector copper core, 19 mm XLPE, 95 mm² screen of copper wires, aluminum barrier and PE jacket, 1 m depth | Twin Falcon ACSR on O155-D2, 0.45 m bundle spacing. Steel earth wire |
| IEEE9 | U230-Cu1200: 127/220–230 kV, five-sector copper core, 21 mm XLPE, 150 mm² screen of copper wires, aluminum barrier and PE jacket, 1 m depth | Cardinal ACSR on O230-F, with ACS 20SA earth wires, 19 wires of 2.5 mm diameter |
| IEEE39 | U345-Cu1600: 200/345 kV, five-sector copper core, 28 mm XLPE, 2.5 mm lead sheath and PE jacket, flat 0.45 m spacing at 1.2 m depth | Cardinal ACSR on O345-F with the same ACS earth wires as O230-F |
| Three- and five-node stars | U155-Cu1600, as above | Falcon ACSR on O100-D, with steel earth wire |

The source defines O13-F and O10-F with phase centers at horizontal coordinates
−0.8, 0 and 0.8 m, mean height 10 m, with no earth wire. O36-F has coordinates
−1.5, 0 and 1.5 m, mean height 10.5 m, with no earth wire. O100-D and O155-D2
have phase centers (−2.5, 14), (0, 18.330) and (2.5, 14) m, and an earth wire at
(0, 25.4) m. O230-F and O345-F have phase centers at −10, 0 and 10 m, mean
height 23.333 m, and proposed earth wires at ±8 m, mean height 30.833 m.
All of these structure definitions await audit, including their sag and support data.
The O33-H formation used by hello-world instead cites audited rows `O33.geom` and
`O33.sag`. Its Goat conductor and U33-630 cables already have audited specs.

The candidate HV cable ratings assume cross-bonding. The pipeline cases assume solid grounding of screens and earth wires at both ends.
They use ideal transposition and omit cross-bonding. Their electrical parameters come
from the source files. The candidate compensation and thermal ratings also await audit. These cases keep the
source ratings and compensation.

## Full-uncertainty IEEE30 inputs

The audit covers the four U33 sizes in touching trefoil and the Bear, Goat, Lion and
Panther conductors on O132-D. It also covers the 7 × 3.15 mm steel earthwire and earth.
Solid grounding and ideal transposition form the fixed premise. `inputs.csv` holds the
nominal values, intervals, classes, bases, sensitivities and coordinate decisions.

Each physical input has one coordinate per family, shared by every branch of that family.
A `yes` marks an established coordinate. A `sweep` marks another nonzero interval for the
next study at 50 and 650 Hz. An input enters the box if any sequence z or y changes by
at least 0.1 percent. Fixed inputs and quantities derived from other inputs have `no`.
Mass, metal area, resistance and lay must satisfy their joint construction constraints.

The intervals below use the cited corpus. The span is fixed at 300 m, and the mechanical intervals use the spread of calculation
values for comparable constructions. The support dimensions remain nominal because clearance limits and survey
reserves do not define their uncertainty. Earth permittivity stays at 10 because the
corpus does not bound it. `EARTH.eps` and `O132.gw.R20` each have one row.
The earthwire resistance interval includes wire diameter, resistivity and lay.
Its nominal remains 3.55 Ω/km, as the cases specify.

The review keeps two earlier U33 discrepancies for a later decision. The 400 mm² diameter
interval allows +0.40 mm, but its comparison-product citation states +0.30 mm. The insulation
interval starts at 7.6 mm, but its cited HD 620 mean-thickness limit is 8.0 mm.

### Wire realization and aggregate limits

The full case selects the ACSR wire inputs. LCM computes their resistance, and the
catalogue R20 remains a check value. The derived resistance envelope is absorbed by
those wire inputs and has no additional coordinate.

EN 50182:2001 5.4 limits conductor diameter to ±1 percent at these sizes, and 5.7.3
limits total mass to ±2 percent. The mass includes the helical wire-length increments.
For a requested point, the builder first applies 5.5.6 with
`lay_outer_eff = min(lay_outer, lay_inner)`. The inner aluminum range ends at 16 and
the steel range starts at 16. Their ordering cannot bind within the registered ranges
except at equality.

A feasible wire realization maps to itself. Otherwise, the builder finds the largest
t in [0, 1] for which both aggregate limits hold along
`d_eff = d_nom + t * (d - d_nom)`. Bisection reduces both diameter deviations along
the same direction. The realized mass and metal area also enter the sag equation.
The input coordinates retain their rectangular bounds. This mapping couples the
physical wires inside the builder, so HPM can keep its parameter box without an
additional constrained-domain protocol.

### Temperature and sag

The phase range is −5 to 75 °C. The cold end transfers eThekwini's minimum conductor
and earthwire temperature to the unloaded ambient state. Its normal 132 kV clearance
basis gives the hot end. Ambient and passive earthwire temperatures range from −5 to
45 °C, with the hot end taken from APAR's ACSR rating conditions. This combines utility
and supplier design states by analogy.
The phase temperature must be at least ambient. The study equates earthwire temperature
to ambient and omits solar heating.

The nominal phase sag of 6.0 m and earthwire sag of 5.4 m belong to the assumed 20 °C
reference state. The utility's setting allowance gives the intervals [6.0, 6.06] m and
[5.4, 5.454] m. For a level span of length L, mass per length m, gravity g, metal area A,
final modulus E and thermal expansion coefficient α, horizontal tension is H = mgL²/(8f).
At the reference state, H₀ = mgL²/(8f₀). The parabolic arc length is L + 8f²/(3L).
Equal unstressed lengths with small elastic and thermal strains give

```math
\frac{8(f^2-f_0^2)}{3L^2}
= \frac{H-H_0}{EA} + \alpha(\theta-\theta_0).
```

`derive.py` solves this equation by bisection. Each temperature sets the corresponding
sag, and mean height is support height minus two thirds of sag. The final modulus
represents a settled conductor under the elastic, still-air approximation. The model
omits ice, wind loads, creep and a change in which metal carries the tension.
The reference temperature and theoretical sags remain assumptions.

The table gives temperature endpoints with nominal mechanical inputs. The last column
encloses the intervals of span, reference sag, mechanical properties, metal area and mass.
The phase enclosure uses independent area and mass limits. A realized construction must
obey their joint mass and wire-geometry constraints. The earthwire calculation keeps
wire diameter, helical mass, area and tension coupled. The span remains 300 m. A span range requires its attachment heights and tensions.
The structural capacity factor does not define a span distribution.

| Code | Mass kg/m | H₀ kN at 20 °C | Sag at −5 °C m | Sag at 75 °C m | Joint enclosing sag interval m |
| --- | --- | --- | --- | --- | --- |
| Bear | 1.213 | 22.3040 | 5.1022 | 7.8853 | 5.0120–8.1012 |
| Goat | 1.489 | 27.3789 | 5.1024 | 7.8850 | 5.0122–8.1011 |
| Lion | 1.095 | 20.1343 | 5.1027 | 7.8847 | 5.0125–8.1007 |
| Panther | 0.974 | 17.9094 | 5.1026 | 7.8849 | 5.0122–8.1010 |

The earthwire has nominal mass 0.428461 kg/m and H₀ of 8.7537 kN.
Its sags at −5 and 45 °C are 4.7918 and 6.0063 m.
The joint interval is 4.7911–6.0600 m.

### Assumptions and transferred values

An asterisk marks an analogy or an unchecked clause, as defined above. The mechanical
spreads describe calculation models and comparable products. They do not give production
tolerances. APAR supplies rating weather but does not tabulate modulus or expansion.
The licensed EN and IEC tables and the Eskom schedules supply the mechanical comparisons.

| Row ids | Assumption, inference or extrapolation |
| --- | --- |
| `U33.cond.d.400`, `.500`, `.630`, `.800` | The earlier audit transfers TFK's comparison-product diameter tolerances to Olex aluminum cables. |
| `U33.cond.R20.400`, `.500`, `.630`, `.800` | The earlier audit transfers the lower end of copper type-test resistance evidence to the aluminum family. |
| `U33.cond.rho.400`, `.500`, `.630`, `.800` | Resistance and circular diameter jointly define equivalent core resistivity. |
| `U33.sc1.t.400`, `.500`, `.630`, `.800`, `U33.sc2.t` | Catalogue radial differences and comparison evidence define the earlier layer-thickness assumptions. |
| `U33.ins.t`, `U33.ins.eps`, `U33.ins.rho`, `U33.ins.tand` | The earlier intervals transfer production evidence and material-family properties to the cable construction. |
| `U33.sc.rho.cond`, `U33.sc.rho.ins`, `U33.sc.eps` | Comparison measurements and the assumed large permittivity describe the semiconducting layers. |
| `U33.scr.d`, `U33.scr.rho`, `U33.scr.R20`, `U33.scr.lay` | The earlier audit transfers wire tolerances and copper-family data to an equivalent tube with the catalogue resistance. |
| `U33.jkt.t.400`, `.500`, `.630`, `.800`, `U33.jkt.eps` | The earlier audit transfers comparison sheath measurements and PVC compound data to U33. |
| `U33.od`, `U33.scr.tube.400`, `.500`, `.630`, `.800`, `U33.jkt.model.400`, `.500`, `.630`, `.800` | Layer sums set diameter, equal resistance sets tube area, and the jacket correction preserves outside diameter. |
| `U33.inst.depth`, `U33.op.theta`, `U33.jkt.rho` | The earlier assumptions cover burial depth, cold soil temperature and the transfer of the insulating-polymer leakage interval to PVC. |
| `O132.span` | eThekwini 4.5.2 gives the fixed 300 m standard span. Its 1.7 structural capacity factor does not describe span uncertainty. |
| `O132.phase.reference`, `O132.gw.reference` | eThekwini 8.7.8.14 permits a one-sided setting allowance of min(1 percent of theoretical sag, 0.100 m). Transfer gives [6, 6.06] m and [5.4, 5.454] m. |
| `O132.reference.theta` | The case's nominal 20 °C defines the reference state of the assumed theoretical sags. |
| `O132.geom` | The dimensions stay at nominal because the corpus lacks their uncertainty. A survey or dimensioned drawing with erection tolerances would supply a basis. |
| `ACSR.E` | Trefinasa's 30/7 value converts to 80.41453 GPa. EN 50182 Table F.3 gives 82 GPa for equal-wire 30/7 AL1/ST1A, defining the transferred model spread. |
| `ACSR.expansion` | The same 30/7 tables give 17.8 and 19.2 × 10⁻⁶/K, defining the transferred model spread. |
| `O132.gw.E` | Eskom Appendix D gives 190 GPa for a 7/3.51 S1A strand. Trefinasa's galvanized-steel value converts to 196.133 GPa, defining the strand-to-material analogy. |
| `O132.gw.expansion` | Trefinasa, Eskom and IEC 63248 Table B.1 agree on 11.5 × 10⁻⁶/K. The coefficient remains nominal because the corpus does not bound its spread. |
| `O132.ambient.theta`, `ACSR.op.theta` | eThekwini 4.7.2 supplies the −5 °C cold state and APAR's rating notes supply 45 °C ambient. eThekwini 4.6.1.3 supplies the phase hot end of 75 °C. |
| `ACSR.Bear.area`, `ACSR.Goat.area`, `ACSR.Lion.area`, `ACSR.Panther.area` | The 37 wire circles define metal area. |
| `O132.phase.H0.Bear`, `.Goat`, `.Lion`, `.Panther`, `O132.phase.sag.Bear`, `.Goat`, `.Lion`, `.Panther` | The calculation obtains tension from parabolic equilibrium and sag from the elastic thermal equation. |
| `O132.gw.area`, `O132.gw.mass`, `O132.gw.H0`, `O132.gw.sag`, `O132.gw.R20` | Wire geometry, helical length, steel properties and the state equation determine the earthwire quantities jointly. |
| `ACSR.al.rho`, `ACSR.st.rho`, `ACSR.st.mu` | The earlier material-family conductivity, steel resistivity and magnetic screening intervals remain. |
| `EARTH.eps` | The value stays at 10 because Gustavsen 2005 bounds earth resistivity and cable-insulation permittivity, while the local papers do not bound soil permittivity. |
| `EARTH.rho`, `EARTH.mu` | The earlier environmental resistivity interval and nonmagnetic earth approximation remain. |
| `MATERIAL.T0`, `MATERIAL.metal.eps`, `MATERIAL.nonmagnetic.mu`, `MATERIAL.dielectric.alpha`, `MATERIAL.tand.default`, `MATERIAL.rho_thermal`, `MATERIAL.theta_max`, `MATERIAL.sigma_solar` | The material calls keep fixed reference conventions and defaults that the electrical calculation does not use. |
| `SYSTEM.length`, `SYSTEM.connections` | The per-meter normalization and the grounding and transposition premise remain fixed. |

Family fields name physical quantities. Case coordinates append the family name, such as
`R20_U33`, `t_ins_U33` and `rho_EARTH`. The `registered` scale identifies a normalized
position: zero is nominal, minus one is the registered lower bound and plus one is the
registered upper bound. One-sided intervals retain the applicable half of that interval.

## Falcon and the physical star cases

The twelve `ACSR.Falcon.*` rows and `STAR3.length`, `STAR5.length` supply the two
physical star cases. The conductor library and strand table in
`local/codex/alternative-test-cases/physical-candidates.md` give the 54 aluminum
wires of 4.36 mm and the 19 steel wires of 2.62 mm, arranged as 12/18/24 aluminum
courses over a steel center and 6/12 courses. The component masses are 2243 and
799 kg/km, from Trefinasa printed p. 4. Wire dimensions give
39.26 mm outside diameter. The candidate library gives 39.24 mm.

The diameter intervals come from the Falcon rows of
`local/codex/manufacturing-tolerances/physical-input-tolerances.csv`: aluminum
4.3164 to 4.4036 mm, steel 2.58 to 2.66 mm. Aggregate diameter and mass limits use
EN 50182 by analogy. The middle aluminum lay adopts the inner-course interval
10 to 16, midpoint 13. The second steel course uses `ACSR.lay.st12`, 14 to 22,
midpoint 18. Outward lays are capped by the preceding course. Every course
contributes to the mass and lay calculation.

Material assumption: the conductor uses the registered ACSR material laws and intervals,
including steel relative permeability 50 to 300, nominal 100.

Assumption: the Trefinasa printed p. 3 calculation values for 54+19 construction,
66.68522 GPa modulus and 0.0000194 per kelvin expansion, are fixed for Falcon.
The 30/7 mechanical coordinates do not act on Falcon.

Installation assumption: O132-D serves the 100 kV stars. It retains the existing
supports and the 300 m span with the reference sags. GW-7x3.15 has both-end grounding,
and each circuit has ideal transposition. Phase temperature is at least ambient temperature.

The candidate study's `scratch/branch-candidates.csv` supplies the shared fitted
lengths: 25690.406433487706 m for star3 and 513.8081286697541 m for star5.
`derive.py` reproduces its joint relative R/X fit. These fixed study lengths
have no surveyed route interval. The physical cases retain LCM charging.

The current rating is an imported result of the candidate study's heat balance.
The same candidate table supplies the physical current rating, 1372.8443238788311 A,
or 237.783612 MVA at 100 kV. The source rating is 100 MVA, or 577.350269 A.
Its IEEE 738 reference conditions are 75 °C conductor, 25 °C ambient, 0.6 m/s
perpendicular wind, sea level, 1000 W/m² solar flux, emissivity and absorptivity
0.5. Assumption: this fixed reference ceiling serves the electrical study without
recomputing a weather-dependent rating at each point.

`sweep.jl` admits eleven coordinates at 0.1 percent for each physical case.
One coordinate per input and family acts on every spoke. Both lay coordinates
remain separate. The mechanical values are fixed.
The source cases remain in `cases/star3.toml` and `cases/star5.toml`.

## Cardinal and the physical transmission cases

The twelve `ACSR.Cardinal.*` rows and `O230.geom`, `O230.phase.reference`,
`O230.gw.reference` supply the IEEE9 and IEEE39 physical cases. The candidate study's
conductor library and strand table specify 54 aluminum and seven steel wires,
all 3.38 mm in diameter. Aluminum courses contain 12, 18 and 24 wires around
the steel center and six-wire course. The resulting diameter is 30.42 mm.
the catalog DC resistance is 0.05992 ohm/km.

The manufacturing tolerance study gives aluminum wire diameters 3.3462 to
3.4138 mm and steel wire diameters 3.33 to 3.43 mm. The existing ACSR wire
realization applies the shared material intervals, lay ordering, aggregate
one percent diameter and two percent mass limits. The middle aluminum course
uses the inner-course lay interval 10 to 16, nominal 13, by analogy.
The steel lay is 16 to 26, nominal 21, and the outer aluminum lay is
10 to 14, nominal 12. There is one coordinate per input and family.

Assumption: aluminum density 2700 kg/m3 and steel density 7780 kg/m3 determine
nominal component masses from wire circles and helical lengths, giving
1.3487363340027083 and 0.49331435514222144 kg/m, respectively.
For wire area A, each course contributes density times A times its wire count
times sqrt(1 + (pi/lay)^2). The straight steel center contributes density times A.
Assumption: the fixed final modulus is 70 GPa and thermal expansion is
0.0000193 per kelvin for this 54/7 construction. The 30/7 mechanical values
are replaced by these values for Cardinal.

O230-F and O345-F share the O230 family. The source geometry gives phase
supports at x = -10, 0, 10 m and height 30 m, earth-wire supports at
x = -3.25, 3.25 m and height 37.5 m, and reference sags of 10 m for both.
The source's 6.5 m earth-wire separation is retained. The candidate study's
16 m separation is superseded by the source geometry. Mean height subtracts
two thirds of the sag. Assumption: two GW-7x3.15 earth wires replace the
source's unspecified earth wires, with both-end grounding and ideal transposition.
Assumption: the towers share the existing span of 300 m and reference temperature.
The ambient interval and earth-wire mechanical and material inputs also apply.
The existing setting allowance gives reference sag intervals 10 to 10.1 m.
The phase temperature remains at least ambient temperature.

The physical cases use every supplied `length_m` in their network JSON.
Each branch assignment cites its fixed source route. The lengths have no
uncertainty coordinate. LCM supplies charging.
The candidate study's fixed thermal current is 994.5173367327438 A, equivalent
to 396.1875479327156 MVA at 230 kV and 594.2813218990733 MVA at 345 kV.
Both source systems rate their lines at 1000 MVA.
Assumption: this fixed current ceiling uses the candidate study's IEEE 738
reference weather without recalculating ampacity at each electrical point.

IEEE39 branch 31 (27-29) keeps its source 1000 MVA rating. Its nominal
fundamental terminal currents are 0.676359489 and 0.677167628 per unit,
above the physical current limit 0.660312580 per unit that HPM derives from
594.2813218990733 MVA and the source minimum voltage of 0.9 per unit.
The source current limit is 1.111111111 per unit. Every other IEEE39 line
keeps the physical rating. The fundamental power flow selects this exception.

The transformer representation stays in the source branch records. IEEE9's
case notes propose Dy11 in generator-to-grid order with a solid grid neutral.
IEEE39's notes propose Yd11 on generator-connected records and Yy0 on internal
records, with solid grid-side grounding. These assumed vector groups remain
notes because the source branch representation does not apply winding groups.
The IEEE39 generator source impedances and capacitive load components remain.

Both sweeps admit ten coordinates at 0.1 percent: phase temperature, earth
resistivity, aluminum resistivity and diameter, inner and middle aluminum lays,
and ambient temperature, wire diameter, steel resistivity and permeability
of the earth-wire family. The remaining registered coordinates stay nominal.
