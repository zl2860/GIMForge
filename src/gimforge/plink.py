"""Narrow, audited wrappers around established PLINK 2 operations."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence


class PlinkError(RuntimeError):
    pass


def require_bfile(prefix: str | Path) -> Path:
    prefix = Path(prefix)
    missing = [str(prefix.with_suffix(suffix)) for suffix in (".bed", ".bim", ".fam") if not prefix.with_suffix(suffix).is_file()]
    if missing:
        raise FileNotFoundError("PLINK binary input is incomplete: " + ", ".join(missing))
    return prefix


def require_executable(path: str | Path) -> Path:
    path = Path(path)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"PLINK2 executable is not executable: {path}")
    return path


def version(plink2: str | Path) -> str:
    completed = subprocess.run([str(plink2), "--version"], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return completed.stdout.strip()


def run(plink2: str | Path, args: Sequence[str | Path], *, context: str, quiet: bool = True) -> subprocess.CompletedProcess[str]:
    command = [str(plink2), *map(str, args)]
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.DEVNULL if quiet else subprocess.PIPE,
        stderr=subprocess.DEVNULL if quiet else subprocess.PIPE,
    )
    if not quiet:
        if completed.stdout:
            print(completed.stdout, file=sys.stderr, end="" if completed.stdout.endswith("\n") else "\n")
        if completed.stderr:
            print(completed.stderr, file=sys.stderr, end="" if completed.stderr.endswith("\n") else "\n")
    if completed.returncode != 0:
        detail = ""
        diagnostic = (completed.stdout or "") + (completed.stderr or "")
        if not quiet and diagnostic:
            detail = "\nPLINK2 output:\n" + diagnostic[-4000:]
        raise PlinkError(f"PLINK2 failed during {context}; rerun with --verbose for diagnostic output.{detail}")
    return completed


def common_filters(*, mac_min: int, geno_missing_max: float | None, chromosome: str | None = "1-22") -> list[str]:
    arguments = ["--snps-only", "just-acgt", "--max-alleles", "2", "--mac", str(mac_min)]
    if chromosome is not None:
        arguments = ["--chr", chromosome, *arguments]
    if geno_missing_max is not None:
        arguments.extend(["--geno", str(geno_missing_max)])
    return arguments


def read_fam_ids(prefix: str | Path) -> set[tuple[str, str]]:
    identifiers: set[tuple[str, str]] = set()
    with Path(prefix).with_suffix(".fam").open() as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                identifiers.add((parts[0], parts[1]))
    return identifiers


def write_ids(path: str | Path, identifiers: Iterable[tuple[str, str]]) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for fid, iid in identifiers:
            handle.write(f"{fid}\t{iid}\n")
