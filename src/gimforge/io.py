"""Small, dependency-free TSV/CSV helpers for compact GIM artefacts."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    return gzip.open(path, mode, newline="") if path.suffix == ".gz" else path.open(mode, newline="")


def detect_delimiter(path: str | Path) -> str:
    suffixes = Path(path).suffixes
    base_suffix = suffixes[-2] if suffixes and suffixes[-1] == ".gz" and len(suffixes) > 1 else Path(path).suffix
    if base_suffix.lower() == ".csv":
        return ","
    return "\t"


def read_table(path: str | Path) -> list[dict[str, str]]:
    return list(iter_table(path))


def iter_table(path: str | Path):
    """Stream rows without materialising a potentially large input table."""
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter=detect_delimiter(path))


def table_header(path: str | Path) -> list[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=detect_delimiter(path))
        return next(reader, [])


def write_table(path: str | Path, rows: Iterable[Mapping[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    materialised = list(rows)
    if fieldnames is None:
        fieldnames = list(materialised[0]) if materialised else []
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(materialised)


def write_json(path: str | Path, value: object) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_json(path: str | Path) -> object:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def as_float(value: object) -> float | None:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    return result


def as_int(value: object) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None
