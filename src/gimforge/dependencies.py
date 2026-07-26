"""Dependency checks with actionable Linux installation guidance."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from .plink import version


PLINK2_INSTALL = "Download PLINK2 from https://www.cog-genomics.org/plink/2.0/ and place the executable in PATH, or pass --plink2 /absolute/path/to/plink2."


def resolve_plink2(requested: str | None = None) -> Path:
    candidate = requested or shutil.which("plink2")
    if not candidate:
        raise ValueError("Required dependency 'plink2' was not found. " + PLINK2_INSTALL)
    path = Path(candidate)
    if not path.is_file() or not path.stat().st_mode & 0o111:
        raise ValueError(f"PLINK2 is not executable: {path}. " + PLINK2_INSTALL)
    return path


def doctor(plink2: str | None = None) -> tuple[bool, list[tuple[str, str, str]]]:
    rows: list[tuple[str, str, str]] = []
    python_ok = sys.version_info >= (3, 10)
    rows.append(("Python", "ok" if python_ok else "missing", sys.version.split()[0]))
    try:
        executable = resolve_plink2(plink2)
        rows.append(("PLINK2", "ok", f"{executable} — {version(executable)}"))
        plink_ok = True
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        rows.append(("PLINK2", "missing", str(error)))
        plink_ok = False
    rows.append(("Python packages", "ok", "none required outside the standard library"))
    return python_ok and plink_ok, rows
