"""Every load-map run in `runs.csv`, with its steady-state points.

Q1 and Q2 are the stage-1 and stage-2 heater powers in W. T1 and T2 are the stage temperatures in K, and T3 is the
Denver stand's pulse-tube thermometer, logged as TPR; Münster figures carry no T3 and the S0028 sweep carries no T1,
so those are NaN. Every run is its own dataset.
"""

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[2] / "data/test_stand"


@dataclass(frozen=True, slots=True)
class Run:
    """One load map: a head on one stand on one date.

    Attributes:
        head: Serial, as the overview names it, or the archive's name for heads the overview lacks.
        stand: "Denver" or "Münster".
        date: ISO date of the map.
        variant: Head build: "regular 1.5", "regular 1.5 HX1 prototype", "IBM", or "unknown".
        compressors: Number of F-70 compressors on the head, or None if unrecorded.
        valve_frequency: Shipped valve frequency in Hz, or None if unknown.
        kind: "grid" for a map over pairs of heater powers, "sweep" for one stage-1 power and a sweep of stage 2.
        note: Provenance and caveats from `runs.csv`.
        q1, q2, t1, t2, t3: One value per steady-state point, in measurement order for Denver runs.
    """

    head: str
    stand: str
    date: str
    variant: str
    compressors: int | None
    valve_frequency: float | None
    kind: str
    note: str
    q1: NDArray[np.float64] = field(repr=False)
    q2: NDArray[np.float64] = field(repr=False)
    t1: NDArray[np.float64] = field(repr=False)
    t2: NDArray[np.float64] = field(repr=False)
    t3: NDArray[np.float64] = field(repr=False)

    @property
    def label(self) -> str:
        """The run's name in figures, such as "S0009 Denver 2026-05-14"."""
        suffix = " sweep" if self.kind == "sweep" else ""
        return f"{self.head} {self.stand} {self.date}{suffix}"


def read_archive_map(path: Path) -> NDArray[np.float64]:
    """Reads the Q1, Q2, T1, T3, and T2 columns of a test-stand load-map file, skipping its header lines."""
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 6 and fields[0].replace(".", "", 1).isdigit():
            rows.append([float(fields[i]) for i in (0, 1, 3, 4, 5)])
    return np.array(rows)


def load_runs(
    directory: Path = DATA, archive: Path = HERE.parents[4] / "data" / "test_stand" / "loadmap_data"
) -> list[Run]:
    """Loads every run in `runs.csv`, plus the S0028 sweep, from the stand archive and `muenster_maps.csv`.

    The 2026-04-13 S0007 map logs its last point at 0.5 W on stage 2; it is read as 0 W, the power its temperature
    and the map's order give.

    Args:
        directory: The folder holding `runs.csv` and `muenster_maps.csv`.
        archive: The stand's load-map archive, which `runs.csv` sources prefixed `archive/` name a file in.

    Returns:
        The runs in `runs.csv` order, with the sweep last.
    """
    with open(directory / "muenster_maps.csv", encoding="utf-8") as file:
        digitized = list(csv.DictReader(file))
    with open(directory / "runs.csv", encoding="utf-8") as file:
        catalog = list(csv.DictReader(file))

    runs = []
    for entry in catalog:
        meta = {
            "head": entry["head"],
            "stand": entry["stand"],
            "date": entry["date"],
            "variant": entry["variant"],
            "compressors": int(entry["compressors"]) if entry["compressors"] else None,
            "valve_frequency": float(entry["valve_frequency"]) if entry["valve_frequency"] else None,
            "note": entry["note"],
        }
        if entry["source"].startswith("archive/"):
            q1, q2, t1, t3, t2 = read_archive_map(archive / entry["source"].removeprefix("archive/")).T
            if entry["source"].endswith("SN007-04-13-2026.dat"):
                q2[-1] = 0.0
            if entry["stand"] == "Münster":
                t3 = np.full_like(t1, np.nan)
            runs.append(Run(**meta, kind="grid", q1=q1, q2=q2, t1=t1, t2=t2, t3=t3))
        else:
            points = [r for r in digitized if r["head"] == entry["head"] and r["date"] == entry["date"]]
            for kind in ("grid", "sweep"):
                chosen = [r for r in points if r["kind"] == kind]
                if not chosen:
                    continue
                q1, q2, t1, t2 = (
                    np.array([float(r[k]) if r[k] else np.nan for r in chosen]) for k in ("Q1", "Q2", "T1", "T2")
                )
                runs.append(Run(**meta, kind=kind, q1=q1, q2=q2, t1=t1, t2=t2, t3=np.full_like(t1, np.nan)))
    return sorted(runs, key=lambda run: run.kind == "sweep")
