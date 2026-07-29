#!/usr/bin/env python3
"""Create lightweight AlphaGenome lookup data plus one score file per GIM SNP.

The complete AlphaGenome score table is intentionally retained as a single
audit file. The portal must not download that large table before a user has
selected a variant, so this script writes a small index and lazily-loadable
per-SNP score files beside it.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = str(Path(__file__).resolve().parent)
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from alphagenome_scopes import row_is_in_scope, target_scopes


PORTAL_ROOT = Path(__file__).resolve().parents[1]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_scores_summary.json")
    parser.add_argument("--index-out", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_index.json")
    parser.add_argument("--scores-dir", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_scores")
    return parser.parse_args()


def numeric(row: dict, *keys: str) -> float | None:
    for key in keys:
        try:
            value = float(row.get(key))
        except (TypeError, ValueError):
            continue
        if value == value:
            return value
    return None


def strongest_by_absolute_ranking(rows: list[dict]) -> tuple[dict | None, float | None]:
    candidates = [
        (row, value)
        for row in rows
        if (value := numeric(row, "rankingScore")) is not None
    ]
    if not candidates:
        return None, None
    row, value = max(candidates, key=lambda item: abs(item[1]))
    return row, value


def summarise_splicing(rows: list[dict]) -> dict | None:
    component_names = {
        "SPLICE_SITES": "spliceSites",
        "SPLICE_SITE_USAGE": "spliceSiteUsage",
        "SPLICE_JUNCTIONS": "spliceJunctions",
    }
    grouped: dict[tuple[str, str], dict] = defaultdict(dict)
    for row in rows:
        component = component_names.get(str(row.get("outputType") or ""))
        raw = numeric(row, "rawScore", "score")
        if component is None or raw is None:
            continue
        gene = str(row.get("geneId") or row.get("geneName") or "unassigned")
        alternate = str(row.get("alternate") or "")
        key = (alternate, gene)
        previous = grouped[key].get(component)
        if previous is None or abs(raw) > abs(previous):
            grouped[key][component] = raw
            grouped[key]["geneId"] = row.get("geneId")
            grouped[key]["geneName"] = row.get("geneName")
            grouped[key]["alternate"] = alternate
    candidates = []
    for components in grouped.values():
        splice_sites = abs(float(components.get("spliceSites") or 0))
        splice_usage = abs(float(components.get("spliceSiteUsage") or 0))
        splice_junctions = abs(float(components.get("spliceJunctions") or 0))
        merged = splice_sites + splice_usage + splice_junctions / 5
        if merged:
            candidates.append(
                {
                    **components,
                    "score": merged,
                    "formula": (
                        "max |splice sites| + max |splice-site usage| + "
                        "max |splice junctions| / 5"
                    ),
                }
            )
    return max(candidates, key=lambda item: item["score"]) if candidates else None


def summarise_scope(rows: list[dict]) -> dict:
    by_modality: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        modality = str(row.get("outputType") or "UNSPECIFIED")
        by_modality[modality].append(row)
    modality_stats = []
    for modality, modality_rows in sorted(by_modality.items()):
        raw_values = [
            abs(value)
            for row in modality_rows
            if (value := numeric(row, "rawScore", "score")) is not None
        ]
        strongest, strongest_ranking = strongest_by_absolute_ranking(
            modality_rows
        )
        modality_stats.append(
            {
                "outputType": modality,
                "nTracks": len(modality_rows),
                "maxAbsRankingScore": (
                    abs(strongest_ranking)
                    if strongest_ranking is not None
                    else None
                ),
                "maxRankingScore": (
                    abs(strongest_ranking)
                    if strongest_ranking is not None
                    else None
                ),
                "rankingScoreAtMaximum": strongest_ranking,
                "rankingDirection": (
                    "negative"
                    if strongest_ranking is not None and strongest_ranking < 0
                    else "positive_or_unsigned"
                    if strongest_ranking is not None
                    else None
                ),
                "maxAbsRawScore": max(raw_values, default=None),
                "topTrack": strongest.get("track") if strongest else None,
                "topBiosample": (
                    strongest.get("biosample") if strongest else None
                ),
                "topGeneId": strongest.get("geneId") if strongest else None,
                "topGeneName": strongest.get("geneName") if strongest else None,
            }
        )
    modality_stats.sort(
        key=lambda item: (
            -(
                item["maxAbsRankingScore"]
                if item["maxAbsRankingScore"] is not None
                else -1
            ),
            item["outputType"],
        )
    )
    strongest, strongest_ranking = strongest_by_absolute_ranking(rows)
    return {
        "nTracks": len(rows),
        "nModalities": len(by_modality),
        "modalities": sorted(by_modality),
        "maxAbsRankingScore": (
            abs(strongest_ranking) if strongest_ranking is not None else None
        ),
        "maxRankingScore": (
            abs(strongest_ranking) if strongest_ranking is not None else None
        ),
        "rankingScoreAtMaximum": strongest_ranking,
        "rankingDirection": (
            "negative"
            if strongest_ranking is not None and strongest_ranking < 0
            else "positive_or_unsigned"
            if strongest_ranking is not None
            else None
        ),
        "topTrack": strongest.get("track") if strongest else None,
        "topBiosample": strongest.get("biosample") if strongest else None,
        "splicingCombined": summarise_splicing(rows),
        "modalityStats": modality_stats,
    }


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    args = parse_args()
    rows = json.loads(args.input.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise SystemExit("AlphaGenome score input must be a JSON array.")
    by_snp: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        snp_id = row.get("snpId")
        scopes = (
            row.get("scopes")
            if isinstance(row.get("scopes"), list)
            else target_scopes(row)
        )
        if snp_id and scopes:
            by_snp[str(snp_id)].append({**row, "scopes": scopes})
    args.scores_dir.mkdir(parents=True, exist_ok=True)
    for stale in args.scores_dir.glob("*.json"):
        if stale.stem not in by_snp:
            stale.unlink()
    index = []
    for snp_id, score_rows in sorted(by_snp.items()):
        atomic_json(args.scores_dir / f"{snp_id}.json", score_rows)
        scopes = {
            scope: summarise_scope(
                [row for row in score_rows if row_is_in_scope(row, scope)]
            )
            for scope in (
                "target",
                "gastric_tissue",
                "gastric_cancer",
                "immune",
            )
        }
        index.append(
            {
                "snpId": snp_id,
                "nTracks": len(score_rows),
                "nGastricTissueTracks": sum(
                    "gastric_tissue" in row.get("scopes", [])
                    for row in score_rows
                ),
                "nGastricCancerTracks": sum(
                    "gastric_cancer" in row.get("scopes", [])
                    for row in score_rows
                ),
                "nImmuneTracks": sum(
                    "immune" in row.get("scopes", []) for row in score_rows
                ),
                "modalities": sorted({row.get("outputType") for row in score_rows if row.get("outputType")}),
                "nModalities": scopes["target"]["nModalities"],
                "maxAbsRankingScore": scopes["target"][
                    "maxAbsRankingScore"
                ],
                "rankingScoreAtMaximum": scopes["target"][
                    "rankingScoreAtMaximum"
                ],
                "splicingCombined": scopes["target"]["splicingCombined"],
                "scopes": scopes,
            }
        )
    atomic_json(args.index_out, index)
    print(f"Published lazy AlphaGenome files for {len(index)} SNPs ({len(rows)} tracks).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
