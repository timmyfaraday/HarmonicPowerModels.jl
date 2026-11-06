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
| `coordinate` | `yes` for a recommended coordinate of the parameter box, `candidate` for a possible one, `derived` or `no` |
| `notes` | remarks and the history of the value |

## Classes

| Nominal value | Definition |
| --- | --- |
| STD | a value that a standard fixes, such as a tabulated thickness or a maximum resistance |
| CAT | a value in the datasheet of a commercial product, with its page |
| DER | a value that `derive.py` computes from STD or CAT values |
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
| EN50182 | BS EN 50182:2001 | 5.4, 5.5.2, 5.5.4 to 5.5.6, 5.7.3, 5.10, Tables 1, 2 and 4 |
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
| TREFINASA | Trefilados de Navarra overhead conductors (2016), p. 6, ACSR to BS 215 Part 2 | <https://www.trefinasa.com/wp-content/uploads/2017/04/overhead_conductors_trefinasa_2016web.pdf> |
| APAR | APAR Industries ACSR table to BS 215 Part 2 | <https://apar.com/wp-content/uploads/2021/02/1.-ACSR.pdf> |
| KEMA | KEMA Laboratories type test certificate 1533-16 (2016) | <http://nci.ae/ttc/5.%20KEMA%201533-16%20%20(3x300%20CU-XLPE-LAT-SWA-PE%2011%20kV).pdf> |
| KIN | Kinectrics type test report K-422899-RC-0001-R00 (2013) | <https://www.nci.ae/ttc/28.%20TT%20REPORT%203X300%20MM2%2015KV%20KINECTRICS%20CANADA%20-%202013.pdf> |
| DEKRA | DEKRA test report 2236358.51-CAT (2019) | <https://www.nci.ae/ttc/60.%20DEKRA-TYPE-TEST-REPORT-IEC-60502-1-4x150-CU-XLPE-SWA-MDPE.pdf> |
| ENEL | Enel Global Standard GSC003 Rev. 2 (2016) | <https://www.eneldistribuicao.com.br/ce/documentos/GSC-003_R-02.pdf> |
| ESKOM | Eskom 240-185000201 Rev. 1, section 3.3.1 | <https://www.etenders.gov.za/home/Download/?blobName=a435f872-aa2f-4993-a7d3-91c2c6e0678a.pdf&downloadedFileName=240-185000201+-+Rev.1+%281%29+Published.pdf> |
| BOREALIS | Borealis data sheets Borlink LS4201R and LE0595 (2024), typical values | the manufacturer's site |
| HFDK-4201 | Dow ENDURANCE HFDK-4201 SC data sheet, typical values | <https://www.dow.com/en-us/pdp.endurance-hfdk-4201-sc-compound-for-cable-systems.127277z.html> |
| NKT | a third-party copy of the NKT 132 kV brochure, which the audit did not check | <https://manuals.plus/m/72406bd69b69737a8f8c0f48714be3adea8787f4cee0890238bc0c237f976a2b.pdf> |
| POWERGRID | Power Grid Corporation of India technical specification, Volume II Part II | <https://apps.powergrid.in/pgciltenders/Files/034a2eab-3098-45ba-9667-3d6f53cc9f4b/Volume-II_Part-II_TS-1.pdf> |
| ETHEKWINI | eThekwini Electricity, the specification for lines on steel lattice towers | <https://www.etenders.gov.za/home/Download/?blobName=d8f1d52e-319c-4122-8f67-bd084afe9e91.pdf&downloadedFileName=3-GenTechSpec-Construct-Steel-Lattice+Power+Lines.pdf> |
| UKPN | UK Power Networks EDS 08-4000 | <https://media.umbraco.io/uk-power-networks/iunhkfsn/eds-08-4000-ehv-network-design.pdf> |
| SPEN | SP Energy Networks OHL-03-099 | <https://www.spenergynetworks.co.uk/userfiles/file/OHL-03-099.pdf> |
| AUSGRID | Ausgrid NS116, Annex C | <https://www.ausgrid.com.au/-/media/Documents/Technical-Documentation/NS/NS116.pdf> |
| COPPERFYI | European Copper Institute note on copper conductivity | <https://help.copper.fyi/hc/en-us/articles/4410229175442-How-does-copper-set-the-standard-for-electrical-conductivity> |
| CALEDONIAN | Caledonian Cables table of thermoplastic properties | <https://www.caledonian-cables.net/Material-Properties/Thermoplastrc.html> |
| SUDA | Wuxi Suda copolymer-coated aluminium tape | <https://www.copolymercoatedaluminumtape.com/copolymer-coated-tape/copolymer-coated-aluminum-tape/plastic-coated-aluminum-tape.html> |
| VIGEPD | Vig Metal EPD-IES-0028115:003 (2026) | <https://www.environdec.com/library/epd28115> |
| ASTM B3 | the ASTM standard for annealed copper wire, by analogy only | not read |
