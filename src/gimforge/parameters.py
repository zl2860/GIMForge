"""Default GIM parameters and explicit, auditable overrides."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from typing import Any, Mapping


@dataclass(frozen=True)
class GIMParameters:
    """Parameters used by the GIM workflow.

    Reference locus and conditional-analysis thresholds are exposed as
    auditable defaults. Genotype QC necessarily depends on the supplied
    representation. The PLINK1 BED backend can apply MAF, MAC, HWE, and
    missingness filters, but it cannot recover imputation INFO from a hard-call
    BED file.
    """

    sentinel_p: float = 1.25e-11
    conditional_p: float = 1.24741348813236e-8
    sentinel_clump_r2: float = 0.1
    sentinel_clump_window_kb: int = 1_000_000
    ld_span_r2: float = 0.1
    cross_metabolite_merge_r2: float = 0.6
    ld_window_kb: int = 1_000_000
    no_ld_half_width_kb: int = 500
    region_padding_kb: int = 250
    mac_min: int = 10
    maf_min: float | None = None
    hwe_p_min: float | None = None
    geno_missing_max: float | None = None
    threads: int = 1
    metabolite_batch_size: int = 1
    force_single_forward_lead: bool = True
    genetic_model: str = "additive"
    regression_model: str = "linear"
    mixed_backend: str = "bolt-lmm"

    def validate(self) -> "GIMParameters":
        for name in (
            "sentinel_p",
            "conditional_p",
            "sentinel_clump_r2",
            "ld_span_r2",
            "cross_metabolite_merge_r2",
        ):
            value = getattr(self, name)
            if not 0 < value < 1:
                raise ValueError(f"{name} must lie strictly between 0 and 1.")
        if self.cross_metabolite_merge_r2 <= self.ld_span_r2:
            raise ValueError("cross_metabolite_merge_r2 must be stricter than ld_span_r2.")
        for name in (
            "sentinel_clump_window_kb",
            "ld_window_kb",
            "no_ld_half_width_kb",
            "region_padding_kb",
            "mac_min",
            "threads",
            "metabolite_batch_size",
        ):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be positive.")
        if self.geno_missing_max is not None and not 0 <= self.geno_missing_max < 1:
            raise ValueError("geno_missing_max must be in [0, 1), or None.")
        if self.maf_min is not None and not 0 < self.maf_min <= 0.5:
            raise ValueError("maf_min must be in (0, 0.5], or None.")
        if self.hwe_p_min is not None and not 0 < self.hwe_p_min <= 1:
            raise ValueError("hwe_p_min must be in (0, 1], or None.")
        if self.genetic_model not in {"additive", "dominant", "recessive"}:
            raise ValueError("genetic_model must be additive, dominant, or recessive.")
        if self.regression_model not in {"linear", "mixed"}:
            raise ValueError("regression_model must be linear or mixed.")
        if self.mixed_backend != "bolt-lmm":
            raise ValueError("mixed_backend currently supports only bolt-lmm.")
        if self.regression_model == "mixed" and self.genetic_model != "additive":
            raise ValueError(
                "BOLT-LMM tests additive allele dosage only. Use regression_model='linear' "
                "for dominant or recessive conditional analysis."
            )
        return self

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def parameters_from_args(overrides: Mapping[str, Any] | None = None) -> GIMParameters:
    """Build validated GIM parameters, rejecting misspelled overrides."""

    overrides = dict(overrides or {})
    permitted = {item.name for item in fields(GIMParameters)}
    unknown = set(overrides).difference(permitted)
    if unknown:
        raise ValueError(f"Unknown GIM parameter(s): {', '.join(sorted(unknown))}")
    return GIMParameters(**overrides).validate()
