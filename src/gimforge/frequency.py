"""Minor-allele-frequency normalisation and reference-study categories."""

from __future__ import annotations

import math

from .io import as_float


def normalise_allele_frequency(value: object) -> float | None:
    frequency = as_float(value)
    if frequency is None or not math.isfinite(frequency) or not 0 <= frequency <= 1:
        return None
    return frequency


def normalise_maf(value: object) -> float | None:
    maf = normalise_allele_frequency(value)
    if maf is None or maf > 0.5:
        return None
    return maf


def maf_from_a1_frequency(value: object) -> float | None:
    frequency = normalise_allele_frequency(value)
    return min(frequency, 1 - frequency) if frequency is not None else None


def classify_maf(value: object) -> str:
    maf = normalise_maf(value)
    if maf is None:
        return "unknown"
    if maf <= 0.01:
        return "rare"
    if maf <= 0.05:
        return "low_frequency"
    return "common"
