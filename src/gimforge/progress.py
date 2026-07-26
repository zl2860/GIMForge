"""Concise console progress messages mirrored to a run log."""

from __future__ import annotations

import sys
import platform
import shlex
from datetime import datetime
from pathlib import Path
from typing import Iterable


_log_path: Path | None = None


def configure_progress_log(path: str | Path | None) -> None:
    """Send subsequent progress messages to ``path`` as well as stderr."""

    global _log_path
    _log_path = Path(path) if path is not None else None
    if _log_path is not None:
        _log_path.parent.mkdir(parents=True, exist_ok=True)


def progress(message: str) -> None:
    """Print one timestamped progress line and append it to the active log."""

    line = f"{datetime.now().astimezone().isoformat(timespec='seconds')} [GIMForge] {message}"
    print(line, file=sys.stderr, flush=True)
    if _log_path is not None:
        with _log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")


def format_bytes(size: int) -> str:
    """Return a compact binary file-size label."""

    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TiB"


def progress_file(label: str, path: str | Path) -> None:
    """Log a resolved file path and its current size without reading contents."""

    resolved = Path(path).resolve()
    state = format_bytes(resolved.stat().st_size) if resolved.is_file() else "missing"
    progress(f"{label}: {resolved} [{state}]")


def start_progress_session(*, version: str, command: Iterable[str], log_path: str | Path) -> None:
    """Write a PLINK-like run header to console and the configured log."""

    progress("=" * 72)
    progress(
        f"GIMForge {version} | Python {platform.python_version()} | "
        f"{platform.system()} {platform.machine()}"
    )
    progress(f"Command: {shlex.join(list(command))}")
    progress(f"Log: {Path(log_path).resolve()}")
    progress("=" * 72)
