#!/usr/bin/env python3
"""Shared AlphaGenome biological-scope and reproducibility settings."""

from __future__ import annotations

from collections.abc import Mapping


SCORING_SEQUENCE_LENGTH = "1MB"
SCORING_MODEL_VERSION = "ALL_FOLDS"
SCORING_SCOPE_FILTER_VERSION = 2
TARGET_SCOPE_NAMES = ("gastric_tissue", "gastric_cancer", "immune")

GASTRIC_TISSUE_TERMS = (
    "stomach",
    "gastric",
    "uberon:0000945",
    "uberon:0001199",
    "uberon:0004222",
)

# AlphaGenome metadata uses biosample names for cell lines, so retain both
# disease phrases and the common gastric-cancer line names that may be returned.
GASTRIC_CANCER_TERMS = (
    "gastric cancer",
    "gastric carcinoma",
    "gastric adenocarcinoma",
    "stomach cancer",
    "stomach carcinoma",
    "stomach adenocarcinoma",
    "ags",
    "kato iii",
    "kato-iii",
    "mkn1",
    "mkn-1",
    "mkn7",
    "mkn-7",
    "mkn28",
    "mkn-28",
    "mkn45",
    "mkn-45",
    "mkn74",
    "mkn-74",
    "nci-n87",
    "nci n87",
    "snu-1",
    "snu-5",
    "snu-16",
    "hs 746t",
    "hgc-27",
    "hgc27",
    "mgc-803",
    "mgc803",
    "nugc",
    "ocum",
    "ycc-",
)

IMMUNE_TERMS = (
    "immune cell",
    "leukocyte",
    "lymphocyte",
    "lymphoblast",
    "t cell",
    "t-cell",
    "t lymphocyte",
    "b cell",
    "b-cell",
    "b lymphocyte",
    "natural killer",
    "nk cell",
    "monocyte",
    "macrophage",
    "neutrophil",
    "eosinophil",
    "basophil",
    "dendritic cell",
    "plasma cell",
    "mast cell",
    "granulocyte",
    "thymocyte",
    "bone marrow",
    "spleen",
    "jurkat",
    "gm12878",
    "gm12891",
    "gm12892",
)

# These recommended scorers have a tissue-independent component that is
# necessary for the official merged splicing prioritisation score.
SHARED_OUTPUT_TYPES = frozenset({"SPLICE_SITES"})


def _text(row: Mapping, *keys: str) -> str:
    return " ".join(str(row.get(key) or "") for key in keys).lower()


def target_scopes(row: Mapping) -> list[str]:
    """Return the allowed portal scopes represented by one tidy score row."""

    output_type = str(row.get("outputType") or row.get("output_type") or "")
    if output_type in SHARED_OUTPUT_TYPES:
        return ["shared"]

    entity = _text(
        row,
        "track",
        "track_name",
        "biosample",
        "biosample_name",
        "biosampleType",
        "biosample_type",
        "ontologyCurie",
        "ontology_curie",
        "gtexTissue",
        "gtex_tissue",
    )
    cancer = any(term in entity for term in GASTRIC_CANCER_TERMS)
    scopes: list[str] = []
    if cancer:
        scopes.append("gastric_cancer")
    elif any(term in entity for term in GASTRIC_TISSUE_TERMS):
        scopes.append("gastric_tissue")
    if any(term in entity for term in IMMUNE_TERMS):
        scopes.append("immune")
    return scopes


def row_is_in_scope(row: Mapping, scope: str) -> bool:
    scopes = row.get("scopes")
    if not isinstance(scopes, list):
        scopes = target_scopes(row)
    if scope == "target":
        return bool(scopes)
    return scope in scopes or "shared" in scopes

