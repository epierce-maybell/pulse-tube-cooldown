"""Reduce each validation system's BoM to node masses by material class."""

import collections
import csv
import re
from pathlib import Path

DATA = Path(__file__).parent
FILES = {"SN09": DATA / "sn09/bom_910-000018.csv", "FP1": DATA / "fp1/bom_910-000040.csv"}
OUTDIR = {"SN09": DATA / "sn09", "FP1": DATA / "fp1"}
GG_PLATE_WEIGHT = {"PRP": 1.0, "RGP": 2.0, "CFP": 2.5, "STP": 2.5, "ICP": 2.0, "MXP": 1.0}
N_SWITCHES = 11
NODES = ["PRP", "RGP", "CFP", "STP", "ICP", "MXP", "DHX", "UNASSIGNED"]
MATS = ["Cu", "SS", "Pb", "Al", "Brass", "G10", "Other"]


def matclass(s):
    """Classifies a BoM material description as one of `MATS`."""
    u = s.upper()
    if re.search(r"COPPER|C10100|C11000", u):
        return "Cu"
    if re.search(r"STAINLESS|304|316|17-4|18-8|\bSS\b", u):
        return "SS"
    if re.match(r"AL\b|ALUMINUM|ALUMINIUM", u):
        return "Al"
    if "G10" in u:
        return "G10"
    if "BRASS" in u:
        return "Brass"
    if "LEAD" in u:
        return "Pb"
    return "Other"


R = []


def add(rx, node, dm=None, flag=""):
    """Registers a name pattern with its node, default material class, and flag; earlier patterns take precedence."""
    R.append((re.compile(rx), node, dm, flag))


add(r"SHIPPING", "EXCLUDED", "Al", "shipping hardware excluded")
add(r"PULSE TUBE, STRIPPED", "EXCLUDED", "Other", "pulse tube body excluded (in cooler curves)")
add(r"PORT, WIRING INSERT|TOP WIRING INSERT", "EXCLUDED", "Other", "room-temperature wiring-port blank excluded")
add(r"SLAB MXC|BASE, 6STK", "EXCLUDED", "Cu", "6STK/slab MXC design not installed")
add(r"^FLANGE, (SHORT|LONG), PRP PT", "PRP", "Cu", "PT-to-PRP flange, N/A assumed Cu")
add(r"^PLATE, PRIMARY|^PRP, 2\.0", "PRP", "Cu", "plate material N/A, OFHC copper (decision)")
add(r"^PLATE, REGENERATOR|^RGP, 2\.0", "RGP", "Cu", "plate material N/A, OFHC copper (decision)")
add(r"^PLATE, COLD FOOT|^CFP, 2\.0", "CFP", "Cu", "plate material N/A, OFHC copper (decision)")
add(r"^PLATE, STILL|^STP, 2\.0", "STP", "Cu", "plate material N/A, OFHC copper (decision)")
add(r"^PLATE, INTERMEDIATE COLD|^ICP, 2\.0", "ICP", "Cu", "plate material N/A, OFHC copper (decision)")
add(r"^PLATE, MIXING CHAMBER|^MXP, 2\.0", "MXP", "Cu", "plate material N/A, OFHC copper (decision)")
add(r"^TT COVER, (\w+)", None, "Cu", "")
add(r"GAS GAP|GGHS", "GGHS", "Cu", "material N/A, assumed Cu")
add(
    r"^SHELL, (BASE|CAP) (UPPER|LOWER)|^FOIL, PARTIAL ETCH|^SINTER, LOWER|^SINTERED FOIL, UPPER",
    "DHX",
    "Other",
    "DHX identification inferred; sinter/foil N/A assumed silver",
)
add(r"^INSERT, UPPER STAGE WIRING", "SPREAD5", "Cu", "spread over PRP,RGP,CFP,STP,ICP (inferred)")
add(r"^BACK PORT COVER, UPPER STAGES", "SPREAD5", "Cu", "spread over 5 upper plates (inferred); N/A assumed Cu")
add(r"^INSERT, BLANK, MXP|^BACK PORT COVER, MXP", "MXP", "Cu", "")
add(
    r"BASE, STILL|CHAMBER, STILL|HX LOWER, STILL|HX UPPER, STILL|INNER TUBE, STILL|OUTER TUBE, STILL|"
    r"HOUSING, HX PUMPING|TUBE, HX PUMPING|LID, STILL|COND HX, STILL|LIGHT BLOCKER.*STILL|"
    r"CLAMP.*HEATER CARTRIDGE, STILL|ADAPTER, STILL LINE",
    "STP",
    "Cu",
    "",
)
add(r"(PANEL|SHIELD).*\bSTP\b|BRACKET, SHIELD.*, STP", "STP", "Cu", "shield material N/A assumed Cu")
add(r"LIGHT BLOCKING.*CFP", "CFP", "Cu", "material N/A assumed Cu")
add(r"LIGHT BLOCKING.*PRP", "PRP", "Al", "")
add(
    r"LIGHT BLOCKING|SHIELD, (DOOR|BOTTOM|FRONT|REAR|SIDE),? ?MXP|PANEL, SIDE MXP|BRACKET, SHIELD.*MXP",
    "MXP",
    "Cu",
    "shield material N/A assumed Cu",
)
add(r"SINTER, (BASE|CAP) MXC", "MXP", "Other", "silver")
add(r"MIXING CHAMBER|MXC|\bMXP\b", "MXP", "Cu", "material N/A assumed Cu")
add(r"COLD FOOT|\bCFP\b", "CFP", "Cu", "material N/A assumed Cu")
add(r"\bPRP\b|PRIMARY PLATE|PRP SHLD", "PRP", "Al", "shield/bracket N/A assumed Al")
add(r"HX ICP PLATE", "ICP", "Cu", "")
add(r"REGENERATOR|REGNERATOR|\bRGP\b", "RGP", "Cu", "material N/A assumed Cu")
add(r"^SUPPORT, (PRIMARY|TOP-PRP)", "PRP", "G10", "")
add(r"^SUPPORT, MIXING", "MXP", "G10", "")
add(r"THERMAL STRAP", "UNASSIGNED", "Cu", "strap end plates not stated")
add(r"MLI SPACER", "UNASSIGNED", "Al", "node not stated")
add(r"HANDLE, SHIELD", "UNASSIGNED", "Cu", "shield stage not stated")
add(r"SHIPPING SUPPORT", "EXCLUDED", "Al", "shipping hardware")
add(r"TRAP, CONC|ISOTHERM, CONC TRAP", "UNASSIGNED", "Other", "concentric trap stage not stated")
add(r"PLHX", "UNASSIGNED", "Other", "PLHX stage not stated")
add(r"MOUNT, PT[12]", "UNASSIGNED", "Other", "PT mount stage not stated")

SKIP = re.compile(
    r"SCREW|BOLT|WASHER|GASKET|O-RING|NUT\b|STUD|PIN,|RETAINER|PLUG|CLAMP, KF|VALVE|HOSE|CASTER|ELEVATOR|ISOLATOR|"
    r"FRAME|VACUUM CHAMBER|PUMP|MANIFOLD|TURBO|BRACKET|HINGE|LATCH|DOOR, VAC|HEAVY|"
    r"PANEL, (FRT|SIDE FRT|FRONT\.|RIGHT SIDE|LEFT SIDE|LOWER REAR|REAR E)|COVER|TRACK|SHELF|PORT, (THERM|VACUUM)|"
    r"PLATE, STAGE|CARRIAGE|ARM,|SPACER, ELEVATOR|GUARD|FLANGE|ADAPTER|TUBE, BUFFER|PLATE, (BOTTOM|TOP), BUFFER|"
    r"LINE|SPRING|MOUNT|HOIST|CLAMP|BLOCK, (90I|CASTER)|WATER|CHILLER|SCANNER|HANDLE|SUPPORT, (EARTH|RAD)|HANGER|"
    r"BEARING|LEAD|GEAR|HANDWHEEL|BELLOWS|TEE|CONNECTOR|GLAND|FITTING|CAP, TRAP|CATCH|SLIDE|RAIL, (ELEV|OFFSET)|"
    r"PANEL, (ELECTRICAL)|PLATE, (LEFT TOP|CENTER TOP|TURBO)|SHAFT|AEROQUIP|TUBE STUB|UNION|MEDIA|TOP PLATE|"
    r"PLATE, CASTER|BLOCK, REAR|NUT PLATE|HOUSING|CONTROL DISC|MOTOR|HEAD|DISC",
    re.I,
)


def parse_mass(s):
    """Returns the mass in kg from a BoM mass field such as "0.25 kg", or None when it holds no mass."""
    m = re.match(r"([\d.]+)\s*kg", s)
    return float(m.group(1)) if m else None


def classify(sysname):
    """Assigns each massed BoM row of a system to a node and material class.

    Returns:
        The assigned rows, and the unmatched rows of at least 0.15 kg that are not skipped hardware.
    """
    seen_gg = set()
    rows = list(csv.DictReader(open(FILES[sysname], encoding="utf-8-sig")))
    out, unm = [], []
    for r in rows:
        m = parse_mass(r["Mass"])
        if not m:
            continue
        q = float(r["Quantity"] or 1)
        name = r["Name"].upper()
        hit = None
        for rx, node, dm, flag in R:
            mm = rx.search(name)
            if mm:
                hit = (node, dm, flag, mm)
                break
        if not hit:
            if m * q >= 0.15 and not SKIP.search(name):
                unm.append((r["Item"], r["Name"], q, r["Material"], m, round(m * q, 3)))
            continue
        node, dm, flag, mm = hit
        if node is None:
            node = mm.group(1).upper()
            if node not in NODES:
                continue
        if sysname == "SN09" and re.search(r"^PANEL, (LEFT|RIGHT) PRP", name):
            q = 1
            flag = "listed qty 8 treated as 1 (implausible multiplier)"
        if node == "GGHS":
            key = r["Part number"]
            if key in seen_gg:
                node = "EXCLUDED"
                flag = "duplicate gas-gap row (first row per part number is the total)"
            else:
                seen_gg.add(key)
                flag = "gas-gap mass allocated to plates per 11-switch layout; " + flag
        fm = r["Material"]
        mat = dm if re.search(r"N/A|None", fm, re.I) else matclass(fm)
        if "SILVER" in fm.upper():
            mat = "Other"
        out.append(
            dict(
                item=r["Item"],
                pn=r["Part number"],
                name=r["Name"],
                qty=q,
                unit_kg=m,
                mass_kg=round(m * q, 4),
                file_material=fm,
                material=mat,
                node=node,
                flag=flag,
            )
        )
    return out, unm


def table(out):
    """Sums assigned masses by node and material class, spreading shared rows across their plates."""
    t = collections.defaultdict(lambda: collections.defaultdict(float))
    for o in out:
        if o["node"] == "SPREAD5":
            for n in ["PRP", "RGP", "CFP", "STP", "ICP"]:
                t[n][o["material"]] += o["mass_kg"] / 5
        elif o["node"] == "GGHS":
            for n, wt in GG_PLATE_WEIGHT.items():
                t[n][o["material"]] += o["mass_kg"] * wt / N_SWITCHES
        elif o["node"] != "EXCLUDED":
            t[o["node"]][o["material"]] += o["mass_kg"]
    return t


if __name__ == "__main__":
    for s in FILES:
        out, unm = classify(s)
        t = table(out)
        od = OUTDIR[s]
        with open(od / "bom_by_node.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["node"] + [m + "_kg" for m in MATS] + ["total_kg"])
            tot = collections.defaultdict(float)
            for n in NODES:
                row = [round(t[n][m], 3) for m in MATS]
                for m, v in zip(MATS, row, strict=True):
                    tot[m] += v
                w.writerow([n] + row + [round(sum(row), 3)])
            w.writerow(["TOTAL"] + [round(tot[m], 3) for m in MATS] + [round(sum(tot.values()), 3)])
        with open(od / "bom_line_assignments.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                ["item", "part_number", "name", "quantity", "extended_mass_kg", "node", "material_class", "flag"]
            )
            for o in out:
                w.writerow([o["item"], o["pn"], o["name"], o["qty"], o["mass_kg"], o["node"], o["material"], o["flag"]])
        print("=====", s)
        print("node", *MATS, "tot", sep="	")
        for n in NODES:
            row = [t[n][m] for m in MATS]
            print(n, *[f"{v:.2f}" for v in row], f"{sum(row):.2f}", sep="	")
        print("TOTAL", *[f"{sum(t[n][m] for n in NODES):.2f}" for m in MATS], sep="	")
