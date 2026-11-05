"""Geometry for `data/<system>/geometry.csv` from a STEP export of the system's CAD assembly.

Two stages, each a subcommand:

    uv run --with cadquery python data/scripts/step_geometry.py extract <system> <assembly.step> <out_dir>
    uv run --with cadquery python data/scripts/step_geometry.py analyze <system> <out_dir>

`extract` loads the assembly once, writes every solid's name, volume, and bounding box to
`<out_dir>/inventory.csv`, and saves the solids named in the system's table as BREP files under
`<out_dir>/parts/`. `analyze` works from those files: shield envelopes and plate positions from the
inventory, gas cavities as the non-solid pieces inside a part's footprint, and tube bores from
cylindrical faces. The report is read by a person and transcribed into `geometry.csv`.

A cavity is found by fusing a group of solids, subtracting the result from a padded box, and keeping the
pieces that fall inside a footprint box: the plan outline of a named part over a stated height. Pieces are
listed largest first; the enclosed cavity is usually the largest, and any piece that spans the full
footprint at a lid's height is the open space around that lid.
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Cylinder
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.gp import gp_Pnt
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TCollection import TCollection_ExtendedString
from OCP.TDataStd import TDataStd_Name
from OCP.TDF import TDF_Label, TDF_LabelSequence
from OCP.TDocStd import TDocStd_Document
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Shape
from OCP.TopTools import TopTools_ListOfShape
from OCP.XCAFDoc import XCAFDoc_DocumentTool

Box = tuple[float, float, float, float, float, float]

# Per system: which solids to save, and how to read the cavities and tubes from them.
#   select: solids whose name matches are saved by `extract`.
#   cavities: label -> (group pattern, footprint pattern, lower pattern, upper pattern); the footprint
#     is the plan outline of the first solid matching the footprint pattern, over the full height of
#     the solids matching the lower and upper patterns.
#   tubes: solids whose cylindrical faces are reported.
SYSTEMS = {
    "fp1": {
        "select": r"STILL, FP|PLHX|CHX|HX ICP|SHELL, (BASE|CAP)|SINTER|FOIL|MIXING CHAMBER|YURT|MXC TO|DIL INLET"
        r"|BELLOWS, \.75|GAS GAP|, PT$|ADAPTER, \.125|FITTING, HX",
        "cavities": {
            "still": (
                r"^(BASE|LID|COND HX), STILL, FP$|^SHELL, PLHX$",
                r"^LID, STILL, FP$",
                r"^BASE, STILL, FP$",
                r"^LID, STILL, FP$",
            ),
            "PLHX shell": (r"^(SHELL|TUBE|HOUSING), PLHX$", r"^SHELL, PLHX$", r"^SHELL, PLHX$", r"^SHELL, PLHX$"),
            "ICP heat exchanger": (
                r"^(BASE|INSERT|FITTING), HX ICP",
                r"^INSERT, HX ICP",
                r"^INSERT, HX ICP",
                r"^BASE, HX ICP",
            ),
            "mixing chamber": (
                r"MIXING CHAMBER, YURT|INTERNAL YURT|ADAPTER, MXC",
                r"^CAP, MIXING CHAMBER",
                r"^BASE, MIXING CHAMBER",
                r"^BASE, MIXING CHAMBER",
            ),
            "DHX lower": (
                r"^SHELL, (BASE|CAP) LOWER$",
                r"^SHELL, CAP LOWER$",
                r"^SHELL, CAP LOWER$",
                r"^SHELL, BASE LOWER$",
            ),
            "DHX upper": (
                r"^SHELL, (BASE|CAP) UPPER$",
                r"^SHELL, CAP UPPER$",
                r"^SHELL, BASE UPPER$",
                r"^SHELL, CAP UPPER$",
            ),
            "gas-gap switch": (
                r"^(CORE|OUTER CORE|SHELL), GAS GAP$",
                r"^OUTER CORE, GAS GAP$",
                r"^OUTER CORE, GAS GAP$",
                r"^OUTER CORE, GAS GAP$",
            ),
        },
        "tubes": r"CHX|TUBE|BELLOWS|SINTER|FOIL|COND HX|SHELL, PLHX",
    },
    "sn09": {
        "select": r"STILL|CHX|HX ICP|SHELL, (BASE|CAP)|SINTER|FOIL|MIXING CHAMBER|WELDED FLAG|BELLOWS, \.75"
        r"|GAS GAP|PULSE TUBE|REGEN|CLAMP, RGP|FITTING, HX",
        "cavities": {
            "still": (
                r"^(BASE|CHAMBER|INNER TUBE|OUTER TUBE|HX UPPER|HX LOWER), STILL$",
                r"^CHAMBER, STILL$",
                r"^CHAMBER, STILL$",
                r"^INNER TUBE, STILL$",
            ),
            "ICP heat exchanger": (
                r"^(BASE|INSERT|FITTING), HX ICP",
                r"^INSERT, HX ICP",
                r"^INSERT, HX ICP",
                r"^BASE, HX ICP",
            ),
            "mixing chamber, welded flag": (
                r"^(BASE|TOP CAP|BOTTOM CAP|FLAG), MXC, WELDED FLAG$",
                r"^TOP CAP, MXC, WELDED FLAG$",
                r"^BOTTOM CAP, MXC, WELDED FLAG$",
                r"^TOP CAP, MXC, WELDED FLAG$",
            ),
            "mixing chamber, 6STK": (
                r"MIXING CHAMBER$",
                r"^CAP, MIXING CHAMBER$",
                r"^BASE, 6STK, MIXING CHAMBER$",
                r"^BASE, 6STK, MIXING CHAMBER$",
            ),
            "DHX lower": (
                r"^SHELL, (BASE|CAP) LOWER$",
                r"^SHELL, CAP LOWER$",
                r"^SHELL, CAP LOWER$",
                r"^SHELL, BASE LOWER$",
            ),
            "DHX upper": (
                r"^SHELL, (BASE|CAP) UPPER$",
                r"^SHELL, CAP UPPER$",
                r"^SHELL, BASE UPPER$",
                r"^SHELL, CAP UPPER$",
            ),
            "gas-gap switch": (
                r"^(CORE|OUTER CORE|SHELL), GAS GAP$",
                r"^OUTER CORE, GAS GAP$",
                r"^OUTER CORE, GAS GAP$",
                r"^OUTER CORE, GAS GAP$",
            ),
        },
        "tubes": r"CHX|TUBE|BELLOWS|SINTER|FOIL",
    },
}

PLATE = re.compile(r"^(PRP|RGP|CFP|STP|ICP|MXP)\b|TOP PLATE")
SHIELD = re.compile(r"^(SHIELD|PANEL), (REAR|FRONT|SIDE|LEFT|RIGHT|BOTTOM|DOOR)|^SHIELD DOOR|^SHIELD ASSY, BOTTOM")


def read_step(path: str) -> tuple[TDocStd_Document, XCAFDoc_DocumentTool, TDF_LabelSequence]:
    """Returns the document with its shape tool and root labels; the document must outlive both."""
    doc = TDocStd_Document(TCollection_ExtendedString("doc"))
    reader = STEPCAFControl_Reader()
    reader.SetNameMode(True)
    if reader.ReadFile(path) != IFSelect_RetDone:
        sys.exit(f"could not read {path}")
    reader.Transfer(doc)
    shape_tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    roots = TDF_LabelSequence()
    shape_tool.GetFreeShapes(roots)
    return doc, shape_tool, roots


def label_name(label: TDF_Label) -> str:
    """The product name stored on an assembly label, or an empty string."""
    name = TDataStd_Name()
    if label.FindAttribute(TDataStd_Name.GetID_s(), name):
        return name.Get().ToExtString()
    return ""


def solids(shape_tool, roots: TDF_LabelSequence):
    """Yields (path, name, solid) for every solid, placed in assembly coordinates."""

    def visit(label: TDF_Label, location: TopLoc_Location, path: list[str]):
        if shape_tool.IsReference_s(label):
            ref = TDF_Label()
            shape_tool.GetReferredShape_s(label, ref)
            yield from visit(ref, location.Multiplied(shape_tool.GetLocation_s(label)), path + [label_name(label)])
            return
        name = label_name(label)
        if shape_tool.IsAssembly_s(label):
            children = TDF_LabelSequence()
            shape_tool.GetComponents_s(label, children)
            for i in range(1, children.Length() + 1):
                yield from visit(children.Value(i), location, path + [name])
            return
        placed = shape_tool.GetShape_s(label).Moved(location)
        explorer = TopExp_Explorer(placed, TopAbs_SOLID)
        while explorer.More():
            yield "/".join(path), name, explorer.Current()
            explorer.Next()

    for i in range(1, roots.Length() + 1):
        yield from visit(roots.Value(i), TopLoc_Location(), [])


def volume(shape: TopoDS_Shape) -> float:
    """Solid volume in mm^3."""
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return props.Mass()


def bbox(shape: TopoDS_Shape) -> Box:
    """Axis-aligned bounding box as (xmin, ymin, zmin, xmax, ymax, zmax) in mm."""
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box, False)
    return box.Get()


def describe(box: Box) -> str:
    """Box extents and z range for the report, in mm."""
    return f"{box[3] - box[0]:6.1f} x {box[4] - box[1]:6.1f} x {box[5] - box[2]:6.1f}  z[{box[2]:.1f}, {box[5]:.1f}]"


def fuse(shapes: list[TopoDS_Shape]) -> TopoDS_Shape:
    """One solid from the union of the shapes, with a small fuzzy tolerance for touching faces."""
    if len(shapes) == 1:
        return shapes[0]
    arguments, tools = TopTools_ListOfShape(), TopTools_ListOfShape()
    arguments.Append(shapes[0])
    for shape in shapes[1:]:
        tools.Append(shape)
    op = BRepAlgoAPI_Fuse()
    op.SetArguments(arguments)
    op.SetTools(tools)
    op.SetFuzzyValue(0.01)
    op.Build()
    return op.Shape()


def footprint_pieces(group: list[TopoDS_Shape], footprint: Box, z_lo: float, z_hi: float) -> list[tuple[float, Box]]:
    """Non-solid pieces inside the footprint box, as (volume in mm^3, bounding box), largest first."""
    fused = fuse(group)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox(fused)
    pad = 5.0
    outer = BRepPrimAPI_MakeBox(
        gp_Pnt(xmin - pad, ymin - pad, zmin - pad), gp_Pnt(xmax + pad, ymax + pad, zmax + pad)
    ).Shape()
    outside = BRepAlgoAPI_Cut(outer, fused).Shape()
    box = BRepPrimAPI_MakeBox(
        gp_Pnt(footprint[0], footprint[1], z_lo), gp_Pnt(footprint[3], footprint[4], z_hi)
    ).Shape()
    common = BRepAlgoAPI_Common(outside, box).Shape()
    pieces = []
    explorer = TopExp_Explorer(common, TopAbs_SOLID)
    while explorer.More():
        pieces.append((volume(explorer.Current()), bbox(explorer.Current())))
        explorer.Next()
    return sorted(pieces, key=lambda piece: -piece[0])


def cylinder_faces(shape: TopoDS_Shape) -> dict[float, float]:
    """Total face area in mm^2 of each cylindrical radius in mm."""
    areas: dict[float, float] = defaultdict(float)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face_s(explorer.Current())
        surface = GeomAdaptor_Surface(BRep_Tool.Surface_s(face))
        if surface.GetType() == GeomAbs_Cylinder:
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, props)
            areas[round(surface.Cylinder().Radius(), 3)] += props.Mass()
        explorer.Next()
    return dict(sorted(areas.items()))


def nearest(boxes: list[Box], anchor: Box) -> Box:
    """The box that overlaps the anchor in plan and whose z center is closest to the anchor's."""
    anchor_z = (anchor[2] + anchor[5]) / 2
    return min((b for b in boxes if overlaps(b, anchor, 1e9)), key=lambda b: abs((b[2] + b[5]) / 2 - anchor_z))


def overlaps(a: Box, b: Box, slack: float = 2.0) -> bool:
    """Whether two boxes overlap in plan and come within `slack` mm of each other in z."""
    return not (a[3] < b[0] or a[0] > b[3] or a[4] < b[1] or a[1] > b[4] or a[5] < b[2] - slack or a[2] > b[5] + slack)


def extract(system: str, step_path: str, out_dir: Path) -> None:
    """Writes the inventory of every solid and saves the system's selected solids as BREP files."""
    select = re.compile(SYSTEMS[system]["select"], re.I)
    parts_dir = out_dir / "parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    doc, shape_tool, roots = read_step(step_path)
    inventory, parts = [], []
    for path, name, solid in solids(shape_tool, roots):
        inventory.append([path, name, volume(solid), *bbox(solid)])
        if select.search(name):
            file_name = f"{len(parts):05d}.brep"
            BRepTools.Write_s(solid, str(parts_dir / file_name))
            parts.append([file_name, name, path])
    with open(out_dir / "inventory.csv", "w", newline="") as f:
        csv.writer(f).writerows(
            [["path", "name", "volume_mm3", "xmin", "ymin", "zmin", "xmax", "ymax", "zmax"]] + inventory
        )
    with open(parts_dir / "index.csv", "w", newline="") as f:
        csv.writer(f).writerows([["file", "name", "path"]] + parts)
    print(f"{len(inventory)} solids inventoried, {len(parts)} saved")


def analyze(system: str, out_dir: Path) -> None:
    """Prints the plates, shield panels, cavity pieces, and tube radii read from an extraction."""
    config = SYSTEMS[system]
    with open(out_dir / "inventory.csv", newline="") as f:
        inventory = [
            (
                r["name"],
                float(r["volume_mm3"]),
                tuple(float(r[k]) for k in ("xmin", "ymin", "zmin", "xmax", "ymax", "zmax")),
            )
            for r in csv.DictReader(f)
        ]

    print("## plates: name, volume cm^3, plan mm, z range")
    for name, vol, box in sorted(inventory, key=lambda r: r[2][2]):
        if PLATE.search(name) and vol > 3e5 and box[5] - box[2] < 70:
            print(f"  {name:28s} {vol / 1e3:8.1f}  {describe(box)}")
    print("## shield panels: name, volume cm^3, box mm, z range")
    for name, vol, box in sorted(inventory, key=lambda r: r[0]):
        if SHIELD.search(name) and vol > 1e5:
            print(f"  {name:36s} {vol / 1e3:8.1f}  {describe(box)}")

    with open(out_dir / "parts" / "index.csv", newline="") as f:
        index = list(csv.DictReader(f))
    shapes: dict[str, TopoDS_Shape] = {}

    def load(file_name: str) -> TopoDS_Shape:
        if file_name not in shapes:
            shape = TopoDS_Shape()
            BRepTools.Read_s(shape, str(out_dir / "parts" / file_name), BRep_Builder())
            shapes[file_name] = shape
        return shapes[file_name]

    def matching(pattern: str) -> list[tuple[str, TopoDS_Shape]]:
        rx = re.compile(pattern, re.I)
        return [(r["name"], load(r["file"])) for r in index if rx.search(r["name"])]

    print("## cavities: non-solid pieces inside each footprint, cm^3")
    for label, (group, footprint, z_from, z_to) in config["cavities"].items():
        anchor_box = bbox(matching(footprint)[0][1])
        spanning = [nearest([bbox(s) for _, s in matching(p)], anchor_box) for p in (z_from, z_to)]
        z_lo, z_hi = min(b[2] for b in spanning), max(b[5] for b in spanning)
        span_box = (anchor_box[0], anchor_box[1], z_lo, anchor_box[3], anchor_box[4], z_hi)
        group_shapes = [s for _, s in matching(group) if overlaps(bbox(s), span_box)]
        pieces = footprint_pieces(group_shapes, anchor_box, z_lo, z_hi)
        metal = sum(volume(s) for s in group_shapes) / 1e3
        plan = f"{anchor_box[3] - anchor_box[0]:.1f} x {anchor_box[4] - anchor_box[1]:.1f}"
        print(f"  {label}: {len(group_shapes)} solids, metal {metal:.2f}, footprint {plan}, z[{z_lo:.1f}, {z_hi:.1f}]")
        for vol, box in pieces[:5]:
            print(f"    piece {vol / 1e3:8.3f}  box {describe(box)}")
        inside = [(n, volume(s) / 1e3) for n, s in matching(r"SINTER|FOIL") if overlaps(bbox(s), span_box, 0.0)]
        if inside:
            print(f"    porous solids inside: {[(n, round(v, 2)) for n, v in inside]}")

    print("## tubes: name, solid cm^3, box mm, cylindrical radius mm -> face area mm^2")
    seen = set()
    for name, shape in matching(config["tubes"]):
        key = (name, round(volume(shape)))
        if key in seen:
            continue
        seen.add(key)
        box = bbox(shape)
        radii = {r: round(a) for r, a in cylinder_faces(shape).items() if a > 30}
        print(f"  {name:48s} {volume(shape) / 1e3:7.3f}  {describe(box)}  {radii}")


if __name__ == "__main__":
    command, system = sys.argv[1], sys.argv[2]
    if command == "extract":
        extract(system, sys.argv[3], Path(sys.argv[4]))
    elif command == "analyze":
        analyze(system, Path(sys.argv[3]))
    else:
        sys.exit(f"unknown command {command}")
