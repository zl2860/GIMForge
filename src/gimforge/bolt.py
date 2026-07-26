"""Narrow wrapper around the optional BOLT-LMM mixed-model backend."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence


class BoltError(RuntimeError):
    pass


def require_executable(path: str | Path) -> Path:
    path = Path(path)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"BOLT-LMM executable is not executable: {path}")
    return path


def version(bolt: str | Path) -> str:
    completed = subprocess.run(
        [str(bolt), "-h"],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    identifying = next((line for line in lines if "BOLT-LMM" in line), None)
    if identifying:
        return identifying
    if lines:
        return lines[0]
    raise BoltError("BOLT-LMM did not print version/help information.")


def run(
    bolt: str | Path,
    args: Sequence[str | Path],
    *,
    context: str,
    quiet: bool = True,
) -> subprocess.CompletedProcess[str]:
    command = [str(bolt), *map(str, args)]
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if not quiet and completed.stdout:
        print(
            completed.stdout,
            file=sys.stderr,
            end="" if completed.stdout.endswith("\n") else "\n",
        )
    if completed.returncode != 0:
        detail = ""
        if not quiet and completed.stdout:
            detail = "\nBOLT-LMM output:\n" + completed.stdout[-4000:]
        raise BoltError(
            f"BOLT-LMM failed during {context}; rerun with --verbose for "
            f"diagnostic output.{detail}"
        )
    return completed
