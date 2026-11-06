"""Reproduce the arithmetic of the physical-input register.

From the repository root, run this script with `python3 -I`.
The script uses only the Python standard library.
It writes `derived.json` beside this file and prints a summary.
Every value of class DER in `inputs.csv` comes from this script.
"""

import csv
import itertools
import json
import math
import os

EPS0 = 8.8541878128e-12  # F/m
MU0 = 4e-7 * math.pi  # H/m. Classical value. The 2019 SI value differs by ~1e-10.
RHO_CU = 1.7241e-8  # ohm*m at 20 C, annealed copper. IEC 60287-1-1:2006 Table 1, IEC 60028 value.
RHO_AL = 2.8264e-8  # ohm*m at 20 C, A1/1350 maximum (IEC 62641 / IEC 60889)
RHO_ST = 1.9157e-7  # ohm*m at 20 C, S1A coated steel calculation value (IEC 63248)
OUT = {}


def coax_c(eps_r, r1, r2):
    """Capacitance per meter of a coaxial dielectric [F/m]."""
    return 2 * math.pi * EPS0 * eps_r / math.log(r2 / r1)


def lay_factor(ratio):
    """Length of a helical wire per unit axial length, lay ratio = pitch / layer diameter."""
    return math.sqrt(1 + (math.pi / ratio) ** 2)


# 1. Equivalent conductor resistivity of compacted cable conductors
# The candidate study tabulated (d_c mm, R20 ohm/km, metal, nominal mm2) from the catalogs.
# R20 is the catalog value. It equals the IEC 60228:2023 Table 4 maximum for every size below.
# The 2004 edition gives this maximum in Table 2.
CABLES = {
    "U33-400": (23.0, 0.0778, "Al", 400), "U33-500": (26.5, 0.0605, "Al", 500),
    "U33-630": (29.9, 0.0469, "Al", 630), "U33-800": (34.2, 0.0367, "Al", 800),
    "U132-800": (34.2, 0.0221, "Cu", 800),
    "U22-Al95": (11.4, 0.320, "Al", 95), "U22-Al150": (14.2, 0.206, "Al", 150),
    "U22-Al240": (18.0, 0.125, "Al", 240), "U22-Al300": (20.1, 0.100, "Al", 300),
    "U22-Al400": (23.0, 0.0778, "Al", 400), "U22-Al500": (26.5, 0.0605, "Al", 500),
    "U22-Al630": (29.9, 0.0469, "Al", 630), "U22-Al1200": (43.8, 0.0247, "Al", 1200),
    "U36-Al400": (23.4, 0.0778, "Al", 400), "U11-Al240": (18.0, 0.125, "Al", 240),
    "U155-Cu1600": (49.6, 0.0113, "Cu", 1600), "U230-Cu1200": (41.6, 0.0151, "Cu", 1200),
    "U345-Cu1600": (49.6, 0.0113, "Cu", 1600),
}
rho_eq = {}
for name, (d, r20, metal, a_nom) in CABLES.items():
    rho_m = RHO_AL if metal == "Al" else RHO_CU
    r = d / 2 * 1e-3
    area_circle = math.pi * r * r  # m2
    R = r20 * 1e-3  # ohm/m
    rho_e = R * area_circle
    a_equiv = rho_m / R * 1e6  # mm2 of straight metal that gives R20
    rho_eq[name] = {
        "d_c_mm": d, "R20_ohm_per_km": r20, "metal": metal, "nominal_mm2": a_nom,
        "rho_eq_ohm_m": rho_e,
        "R_if_solid_core_with_metal_rho_ohm_per_km": rho_m / area_circle * 1e3,
        "error_solid_core_pct": 100 * (rho_m / area_circle / R - 1),
        "R_if_nominal_area_straight_ohm_per_km": rho_m / (a_nom * 1e-6) * 1e3,
        "error_nominal_area_pct": 100 * (rho_m / (a_nom * 1e-6) / R - 1),
        "equivalent_straight_metal_mm2": a_equiv,
        "fill_factor_equivalent": a_equiv / (area_circle * 1e6),
        "fill_factor_if_nominal_area": a_nom / (area_circle * 1e6),
    }
OUT["conductor_equivalent_resistivity"] = rho_eq

# Production evidence for R20 relative to the IEC 60228 maximum (type-test reports).
OUT["R20_evidence_ratio_to_cap"] = {
    "Kinectrics K-422899 (Cu 300, cap 0.0601)": [0.0552 / 0.0601, 0.0554 / 0.0601],
    "KEMA 1533-16 (Cu 300, cap 0.06010)": [0.059820 / 0.06010, 0.060049 / 0.06010],
    "DEKRA 2236358.51-CAT (Cu 150 shaped, cap 0.124)": [0.116 / 0.124, 0.117 / 0.124],
}

# 2. Catalog-implied XLPE permittivity, Olex 19/33 kV Al (pp. 12-13)
# The tuple is (d_c, t_sc1, t_ins, C uF/km).
# t_sc1 = (diameter over insulation - d_c)/2 - t_ins.
U33 = {"U33-400": (23.0, 40.7, 8.0, 0.267), "U33-500": (26.5, 44.2, 8.0, 0.297),
       "U33-630": (29.9, 48.0, 8.0, 0.329), "U33-800": (34.2, 52.3, 8.0, 0.366)}
eps_implied = {}
for name, (d, d_ins, t, c) in U33.items():
    t_sc1 = (d_ins - d) / 2 - t
    r1 = (d / 2 + t_sc1) * 1e-3
    r2 = r1 + t * 1e-3
    eps = c * 1e-9 * math.log(r2 / r1) / (2 * math.pi * EPS0)
    eps_implied[name] = {"t_sc1_mm": round(t_sc1, 3), "r1_mm": r1 * 1e3, "r2_mm": r2 * 1e3,
                         "eps_r_implied": eps}
OUT["xlpe_eps_r_implied_by_catalogue"] = eps_implied

# Catalog insulation resistance (Olex 400 mm2: 9100 MOhm*km) -> volume resistivity.
r1, r2 = eps_implied["U33-400"]["r1_mm"], eps_implied["U33-400"]["r2_mm"]
OUT["xlpe_rho_implied_by_catalogue_ohm_m"] = 9100e6 * 1e3 * 2 * math.pi / math.log(r2 / r1)

# 3. Copper wire screen: catalog-implied lay factor and lay ratio
screens = {}
for label, (n, dw, r_s_cat) in {"48x1.35 (68.7 mm2), Olex R_s 0.263": (48, 1.35, 0.263),
                                "43x1.35 (61.5 mm2), Olex R_s 0.293": (43, 1.35, 0.293)}.items():
    area = n * math.pi * dw * dw / 4  # mm2
    k = r_s_cat * 1e-3 * area * 1e-6 / RHO_CU
    screens[label] = {"area_mm2": area, "straight_R_ohm_per_km": RHO_CU / (area * 1e-6) * 1e3,
                      "lay_factor_implied": k, "lay_ratio_implied": math.pi / math.sqrt(k * k - 1)}
OUT["screen_lay_implied_by_catalogue"] = screens
# The candidate study's lay ratio 20 [16, 24] against the catalog resistance.
OUT["screen_R_with_candidate_lay_ratio_20"] = RHO_CU / (48 * math.pi * 1.35 ** 2 / 4 * 1e-6) * lay_factor(20) * 1e3

# 4. Sensitivities of the principal 33 kV cable, U33-500
d_c, t_sc1, t_ins, t_sc2, d_w, t_j = 26.5, 0.85, 8.0, 0.90, 1.35, 2.95
eps = 2.40


def u33_c(d_c=d_c, t_sc1=t_sc1, t_ins=t_ins, eps=eps):
    r1 = (d_c / 2 + t_sc1) * 1e-3
    return coax_c(eps, r1, r1 + t_ins * 1e-3)


c0 = u33_c()
sens = {"C_nominal_uF_per_km": c0 * 1e9,
        "dC_pct_t_ins_minus5": 100 * (u33_c(t_ins=0.95 * t_ins) / c0 - 1),
        "dC_pct_t_ins_plus5": 100 * (u33_c(t_ins=1.05 * t_ins) / c0 - 1),
        "dC_pct_t_ins_plus10_starter": 100 * (u33_c(t_ins=1.10 * t_ins) / c0 - 1),
        "dC_pct_eps_2.25": 100 * (u33_c(eps=2.25) / c0 - 1),
        "dC_pct_eps_2.50": 100 * (u33_c(eps=2.50) / c0 - 1),
        "dC_pct_d_c_plus_0.4mm": 100 * (u33_c(d_c=d_c + 0.4) / c0 - 1),
        "dC_pct_t_sc1_plus10": 100 * (u33_c(t_sc1=1.1 * t_sc1) / c0 - 1)}
# Series share of the reciprocal capacitance for the semiconducting layer.
# This uses the lossless law, LCM's default.
r_c = d_c / 2
r_sc = r_c + t_sc1
share = (math.log(r_sc / r_c) / 1000) / (math.log((r_sc + t_ins) / r_sc) / eps)
sens["semicon_lossless_series_share_eps1000"] = share
# Radial resistance relative to insulation reactance at 650 Hz under the lossy law.
# This uses the worst-case IEC limit of 1000 ohm*m.
R_sc = 1000 * math.log(r_sc / r_c) / (2 * math.pi)  # ohm*m
X_c = 1 / (2 * math.pi * 650 * c0)  # ohm*m
sens["semicon_lossy_R_over_Xc_650Hz_rho1000"] = R_sc / X_c
# Jacket -> touching-trefoil spacing -> loop inductance (screening estimate, screens open).
OD = d_c + 2 * (t_sc1 + t_ins + t_sc2 + d_w + t_j)
gmr = 0.7788 * r_c
x0 = 2 * math.pi * 50 * 2e-7 * math.log(OD / gmr) * 1e3  # ohm/km
OD_hi = OD + 2 * 0.35 * t_j
sens["OD_nominal_mm"] = OD
sens["X50_screening_ohm_per_km"] = x0
sens["dX_pct_jacket_plus35"] = 100 * (math.log(OD_hi / gmr) / math.log(OD / gmr) - 1)
# Solenoid term of helically laid screen wires relative to the core-screen coaxial term.
r_s = (OD / 2 - t_j - d_w / 2)  # mm, wire-center radius
for lam in (10, 16, 20):
    pitch = lam * 2 * r_s
    ratio = 2 * math.pi ** 2 * (r_s ** 2 - r_c ** 2) / (pitch ** 2 * math.log(r_s / r_c))
    sens[f"solenoid_over_coaxial_inductance_lay{lam}"] = ratio
OUT["sensitivity_U33_500"] = sens

# 5. Steel earthwire 7/3.15 mm with the IEC 63248 calculation resistivity
# The lay ratio of the 6-wire layer of a 7-wire steel strand is 16-26, nominal 21.
# EN 50182 Table 1 and Eskom 240-185000201 3.3.1 give this range.
# The interval covers rho [1.79e-7, 2.155e-7] (register ACSR.st.rho) and the lay range.
a1 = math.pi * 3.15 ** 2 / 4 * 1e-6
gw_R = lambda rho, lam: rho / (a1 * (1 + 6 / lay_factor(lam))) * 1e3
OUT["earthwire_7x3.15_R20_ohm_per_km_IEC_rho"] = gw_R(RHO_ST, 21)
OUT["earthwire_7x3.15_R20_ohm_per_km_candidate_rho_1.8e-7"] = gw_R(1.8e-7, 21)
OUT["earthwire_7x3.15_R20_interval_ohm_per_km"] = [gw_R(1.79e-7, 26), gw_R(2.155e-7, 16)]
OUT["eskom_10-S1A-7/3.51_implied_rho_incl_lay_ohm_m"] = 2.86e-3 * 7 * math.pi * 3.51 ** 2 / 4 * 1e-6

# 6. ACSR resistance interval from documented component limits
# Conductor mass +/-2 % (Enel GSC003 5.7) bounds the metal area.
# Al rho is in [61.8, 61.0] % IACS.
rho_lo = 1.7241e-8 / 0.618
OUT["ACSR_R20_interval_factor"] = [rho_lo / RHO_AL / 1.02, 1 / 0.98]

# 7. Alternatives. 345 kV single Cardinal surface gradient and 132 kV 50 mm2 screen
r = 30.42 / 2 / 10  # cm
deq = (10 * 10 * 20) ** (1 / 3) * 100  # cm
u_ph = 362 / math.sqrt(3)  # kV rms at Um
OUT["Cardinal_345kV_avg_surface_gradient_kV_per_cm_at_Um"] = u_ph / (r * math.log(deq / r))
# Adiabatic short-circuit capacity of a copper screen, 80 -> 250 C, 1 s:
#   I^2 t = K^2 S^2 ln((theta_f + beta) / (theta_i + beta)),  K = sqrt(sigma (beta + 20) / rho_20).
# IEC 60949 tabulates K = 226 A*s^0.5/mm2 and beta = 234.5 K for copper (not in the licensed set).
# The same K follows from licensed tables alone.
# IEC 60865-1:2011 Annex A gives sigma = c * density of copper (390 J/(kg K), 8900 kg/m3).
# IEC 60228:2023 Annex B gives beta through k = 254.5 / (234.5 + t).
# This is 1 / (1 + 0.00393 (t - 20)). IEC 60287-1-1:2006 Table 1 gives rho_20.
# IEC 60865-1's own bare-conductor data gives the bound below.
# It uses kappa_20 = 56e6 S/m and alpha_20 = 0.0039/K.
# The temperatures 80 C (initial) and 250 C (final) are assumptions.
# IEC 60502-2 refers the short-circuit limits to IEC 60986, which is not in the licensed set.
beta = 234.5
sigma_cu = 390 * 8900  # J/(K m3), IEC 60865-1:2011 Annex A
K = math.sqrt(sigma_cu * (beta + 20) / RHO_CU) / 1e6  # A*s^0.5/mm2
K_60865 = math.sqrt(56e6 * 390 * 8900 / 0.0039) / 1e6  # A*s^0.5/mm2, IEC 60865-1 copper data
OUT["Cu_adiabatic_K_licensed_A_s05_per_mm2"] = K
OUT["Cu_adiabatic_K_over_IEC60949_226"] = K / 226
OUT["Cu_adiabatic_K_IEC60865_data_A_s05_per_mm2"] = K_60865
OUT["Cu_screen_50mm2_adiabatic_1s_kA"] = K * 50 * math.sqrt(math.log((250 + beta) / (80 + beta))) / 1e3
OUT["Cu_screen_50mm2_adiabatic_1s_kA_IEC60865_data"] = (
    K_60865 * 50 * math.sqrt(math.log((1 + 0.0039 * 230) / (1 + 0.0039 * 60))) / 1e3)

# 8. Checks against the licensed sources (bibliography)
# CIGRE TB 531 Table 6 gives the XLPE IEC 60287-1-1 value of 2.5.
# It gives the EPRI typical value of 2.2 and range of 2.1-2.3.
OUT["dC_pct_U33_500_eps_2.10_vs_2.40"] = 100 * (u33_c(eps=2.10) / c0 - 1)
OUT["resonance_frequency_factor_eps_2.10_vs_2.40"] = math.sqrt(2.40 / 2.10)
# NGED EE SPEC 77/5: design stress <= 8 kV/mm at the conductor screen, <= 4 kV/mm at the
# insulation screen. U132-800 with the NGED semiconducting thickness of 1.0 mm.
u0 = 132 / math.sqrt(3)
for t_sc in (1.0, 1.2):
    r1 = 34.2 / 2 + t_sc
    r2 = r1 + 16.0
    OUT[f"U132_stress_kV_per_mm_sc{t_sc}"] = [u0 / (r1 * math.log(r2 / r1)), u0 / (r2 * math.log(r2 / r1))]
    OUT[f"U132_C_uF_per_km_eps2.4_sc{t_sc}"] = coax_c(2.4, r1 * 1e-3, r2 * 1e-3) * 1e9
# Copper screen area for the NGED lead-sheath fault ratings, 1 s, adiabatic law of section 7.
fac = K * math.sqrt(math.log((250 + beta) / (80 + beta)))
OUT["Cu_screen_mm2_for_1s_fault_kA"] = {f"{i}": i * 1e3 / fac for i in (21.5, 25, 31.5)}

# IEC 60840 10.10 (and IEC 62067 10.10): measured capacitance <= 1.08 x declared nominal.
# U132-800 has a declared capacitance of 0.212 uF/km.
# Box corners have eps_r in [2.10, 2.50] and t_ins in [15.2, 16.8] mm.
r1 = 34.2 / 2 + 1.0
corners = {f"eps{e}_t{t}": coax_c(e, r1 * 1e-3, (r1 + t) * 1e-3) * 1e9 / 0.212
           for e in (2.10, 2.50) for t in (15.2, 16.8)}
OUT["U132_C_over_declared_at_box_corners"] = corners

# 9. Full IEEE30 intervals and elastic sag, using the register as the input record.
# The parabolic arc length is L + 8*f**2/(3*L), H = m*g*L**2/(8*f).
# Equal unstressed lengths at two states, with small elastic and thermal strains, give
# 8*(f**2-f0**2)/(3*L**2) = (H-H0)/(E*A) + alpha*(theta-theta0).
# This study uses level supports in still air and an elastic conductor state.
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "inputs.csv"), newline="") as source:
    rows = {row["id"]: row for row in csv.DictReader(source)}


def number(name, column="nominal"):
    """Read a numerical register value in its stated unit."""
    return float(rows[name][column])


def sag(span, mass, area, modulus, expansion, reference, theta, theta0):
    """Solve the elastic state equation for positive sag in meters."""
    load = mass * number("O132.gravity")
    tension0 = load * span**2 / (8 * reference)
    low, high = 1e-9, span
    for _ in range(100):
        value = (low + high) / 2
        tension = load * span**2 / (8 * value)
        residual = (8 * (value**2 - reference**2) / (3 * span**2)
                    - (tension - tension0) / (modulus * area)
                    - expansion * (theta - theta0))
        if residual > 0:
            high = value
        else:
            low = value
    return (low + high) / 2


full = {"core_resistivity_ohm_m": {}, "cable_diameter_mm": {},
        "acsr": {}, "earthwire": {}, "screen_tube_mm": {}, "jacket_model_mm": {}}
for size in (400, 500, 630, 800):
    diameter = "U33.cond.d." + str(size)
    resistance = "U33.cond.R20." + str(size)
    full["core_resistivity_ohm_m"][str(size)] = [
        number(resistance, end) * 1e-3 * math.pi * (number(diameter, end) * 1e-3 / 2)**2
        for end in ("lower", "upper")]
    full["cable_diameter_mm"][str(size)] = [
        number(diameter, end) + 2 * (number("U33.sc1.t." + str(size), end)
        + number("U33.ins.t", end) + number("U33.sc2.t", end)
        + number("U33.scr.d", end) + number("U33.jkt.t." + str(size), end))
        for end in ("lower", "upper")]

    tubes, jackets = [], []
    for column in ("nominal", "lower", "upper"):
        dimensions = []
        for ends in itertools.product((column,) if column == "nominal" else ("lower", "upper"), repeat=8):
            core, inner, ins, outer, copper, screen_r, jacket, wire = [
                number(name, end) for name, end in zip(
                    (diameter, "U33.sc1.t." + str(size), "U33.ins.t", "U33.sc2.t",
                     "U33.scr.rho", "U33.scr.R20", "U33.jkt.t." + str(size), "U33.scr.d"), ends)]
            radius = core / 2 + inner + ins + outer
            tube = math.sqrt(radius**2 + copper * 1e9 / (math.pi * screen_r)) - radius
            dimensions.append((tube, jacket + wire - tube))
        if column == "nominal":
            tubes.append(dimensions[0][0])
            jackets.append(dimensions[0][1])
        else:
            select = min if column == "lower" else max
            tubes.append(select(item[0] for item in dimensions))
            jackets.append(select(item[1] for item in dimensions))
    full["screen_tube_mm"][str(size)] = tubes
    full["jacket_model_mm"][str(size)] = jackets


def states(prefix, mass, area, reference, temperature, areas):
    """Return reference tension, temperature endpoints and the joint interval."""
    span = number("O132.span")
    modulus = number(prefix + ".E") * 1e9
    expansion = number(prefix + ".expansion")
    f0 = number(reference)
    t0 = number("O132.reference.theta")
    nominal = [sag(span, mass, area, modulus, expansion, f0,
                   number(temperature, end), t0) for end in ("lower", "upper")]
    bounds = []
    # Enclose mass and area independently here. The result is a conservative bound,
    # not an instruction to split their shared physical causes into coordinates.
    for length, factor, section, elastic, thermal, initial, theta in itertools.product(
            (number("O132.span", end) for end in ("lower", "upper")),
            (0.98, 1.02), areas,
            (number(prefix + ".E", end) * 1e9 for end in ("lower", "upper")),
            (number(prefix + ".expansion", end) for end in ("lower", "upper")),
            (number(reference, end) for end in ("lower", "upper")),
            (number(temperature, end) for end in ("lower", "upper"))):
        bounds.append(sag(length, mass * factor, section, elastic, thermal,
                          initial, theta, t0))
    return {"area_m2": area, "mass_kg_per_m": mass,
            "H0_N": mass * number("O132.gravity") * span**2 / (8 * f0),
            "H0_interval_N": [mass * factor * number("O132.gravity")
                * number("O132.span", end)**2 / (8 * number(reference, other))
                for factor, end, other in ((0.98, "lower", "upper"), (1.02, "upper", "lower"))],
            "temperature_end_sag_m": nominal, "joint_sag_interval_m": [min(bounds), max(bounds)]}


for code in ("Bear", "Goat", "Lion", "Panther"):
    d = float(rows["ACSR." + code + ".al"]["nominal"].split(" x ")[1]) * 1e-3
    area = 37 * math.pi * d**2 / 4
    areas = [sum(int(rows["ACSR." + code + "." + metal][end].split(" x ")[0])
        * math.pi * (float(rows["ACSR." + code + "." + metal][end].split(" x ")[1]) * 1e-3)**2 / 4
        for metal in ("al", "st")) for end in ("lower", "upper")]
    full["acsr"][code] = states("ACSR", number("ACSR." + code + ".mass"), area,
                                "O132.phase.reference", "ACSR.op.theta", areas)
    full["acsr"][code]["area_interval_m2"] = areas
# Earthwire area and mass use the same wire diameter and lay, not an extra mass factor.
gw = full["earthwire"]
d = number("O132.gw.d") * 1e-3
area = 7 * math.pi * d**2 / 4
mass = number("O132.gw.density") * math.pi * d**2 / 4 * (1 + 6 * lay_factor(21))
gw.update(states("O132.gw", mass, area, "O132.gw.reference", "O132.ambient.theta", (area,)))
# Recompute its joint interval with diameter and lay tied to its actual steel area and mass.
values = []
weights = []
tensions = []
for length, diameter, lay, elastic, thermal, initial, theta in itertools.product(
        (number("O132.span", end) for end in ("lower", "upper")),
        (number("O132.gw.d", end) * 1e-3 for end in ("lower", "upper")),
        (16, 26),
        (number("O132.gw.E", end) * 1e9 for end in ("lower", "upper")),
        (number("O132.gw.expansion", end) for end in ("lower", "upper")),
        (number("O132.gw.reference", end) for end in ("lower", "upper")),
        (number("O132.ambient.theta", end) for end in ("lower", "upper"))):
    wire_area = math.pi * diameter**2 / 4
    weight = number("O132.gw.density") * wire_area * (1 + 6 * lay_factor(lay))
    weights.append(weight)
    tensions.append(weight * number("O132.gravity") * length**2 / (8 * initial))
    values.append(sag(length, weight, 7 * wire_area, elastic, thermal, initial, theta,
                      number("O132.reference.theta")))
gw["joint_sag_interval_m"] = [min(values), max(values)]
gw["mass_interval_kg_per_m"] = [min(weights), max(weights)]
gw["H0_interval_N"] = [min(tensions), max(tensions)]
gw["R20_interval_ohm_per_km"] = [
    number("ACSR.st.rho", end) / (math.pi * (number("O132.gw.d", opposite) * 1e-3)**2 / 4
    * (1 + 6 / lay_factor(number("O132.gw.lay", opposite)))) * 1e3
    for end, opposite in (("lower", "upper"), ("upper", "lower"))]
# Reproduce the arithmetic of the corpus transfers, distinct from production limits.
full["corpus_interval_derivations"] = {
    "span_m": [300, 300],  # ETHEKWINI 4.5.2 standard span.
    "phase_reference_sag_m": [6.0, 6.0 + min(0.01 * 6.0, 0.100)],
    "earthwire_reference_sag_m": [5.4, 5.4 + min(0.01 * 5.4, 0.100)],
    "ACSR_modulus_GPa": [8200 * number("O132.gravity") / 1000, 82000 / 1000],
    "ACSR_expansion_per_K": [17.8e-6, 19.2e-6],
    "earthwire_modulus_GPa": [190000 / 1000, 20000 * number("O132.gravity") / 1000],
    "earthwire_expansion_per_K": [11.5e-6, 11.5e-6],
    "ambient_C": [-5, 45],  # ETHEKWINI 4.7.2 and APAR rating-table notes.
}
OUT["full_uncertainty_ieee30"] = full

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "derived.json"), "w") as f:
    json.dump(OUT, f, indent=2)

print("rho_eq and resistance error of a solid core at catalogue diameter with metal rho:")
for k, v in rho_eq.items():
    print(f"  {k:13s} rho_eq={v['rho_eq_ohm_m']:.4e}  solid-core err={v['error_solid_core_pct']:+6.1f}%"
          f"  nominal-area err={v['error_nominal_area_pct']:+5.1f}%  fill(eq)={v['fill_factor_equivalent']:.3f}"
          f"  fill(nominal)={v['fill_factor_if_nominal_area']:.3f}")
for k, v in eps_implied.items():
    print(f"  {k}: t_sc1={v['t_sc1_mm']} mm, eps_r implied={v['eps_r_implied']:.3f}")
print("  XLPE rho implied by 9100 MOhm*km:", f"{OUT['xlpe_rho_implied_by_catalogue_ohm_m']:.3e}")
for k, v in screens.items():
    print(f"  screen {k}: k={v['lay_factor_implied']:.4f}, lay ratio={v['lay_ratio_implied']:.1f}")
print("  screen R with lay ratio 20:", f"{OUT['screen_R_with_candidate_lay_ratio_20']:.4f} ohm/km")
for k, v in sens.items():
    print(f"  {k}: {v:.4g}")
for k in ("earthwire_7x3.15_R20_ohm_per_km_IEC_rho", "earthwire_7x3.15_R20_ohm_per_km_candidate_rho_1.8e-7",
          "earthwire_7x3.15_R20_interval_ohm_per_km",
          "eskom_10-S1A-7/3.51_implied_rho_incl_lay_ohm_m", "ACSR_R20_interval_factor",
          "Cardinal_345kV_avg_surface_gradient_kV_per_cm_at_Um", "Cu_adiabatic_K_licensed_A_s05_per_mm2",
          "Cu_adiabatic_K_over_IEC60949_226", "Cu_adiabatic_K_IEC60865_data_A_s05_per_mm2",
          "Cu_screen_50mm2_adiabatic_1s_kA", "Cu_screen_50mm2_adiabatic_1s_kA_IEC60865_data",
          "Cu_screen_mm2_for_1s_fault_kA"):
    print(f"  {k}: {OUT[k]}")
for k, v in OUT["R20_evidence_ratio_to_cap"].items():
    print(f"  {k}: {v[0]:.3f}-{v[1]:.3f}")
