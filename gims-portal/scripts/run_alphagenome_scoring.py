#!/usr/bin/env python3
"""Score every study ALT with the authorized AlphaGenome API.

The browser never receives an API key. Run this local job with an authorized key:

    python3.12 -m venv .venv-alphagenome
    .venv-alphagenome/bin/pip install "git+https://github.com/google-deepmind/alphagenome.git" pandas
    export ALPHAGENOME_API_KEY='...'
    .venv-alphagenome/bin/python scripts/run_alphagenome_scoring.py

It uses the official 1 MB context and recommended variant scorers, pins the
ALL_FOLDS model, retains only gastric tissue, gastric-cancer, immune-cell and
shared splicing outputs, resumes after interruption, and writes a reproducible
coverage/status manifest after each checkpoint.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import os
import time
from datetime import UTC, datetime
from pathlib import Path

from alphagenome_scopes import (
    SCORING_MODEL_VERSION,
    SCORING_SCOPE_FILTER_VERSION,
    SCORING_SEQUENCE_LENGTH,
    TARGET_SCOPE_NAMES,
    target_scopes,
)


PORTAL_ROOT = Path(__file__).resolve().parents[1]


def required_packages():
    try:
        import pandas as pd
        from alphagenome.data import genome
        from alphagenome.models import dna_client, dna_model, variant_scorers
    except ImportError as error:
        raise SystemExit(
            "AlphaGenome scoring requires its official SDK and pandas. "
            "Install with: pip install 'git+https://github.com/google-deepmind/alphagenome.git' "
            f"pandas numpy. Import error: {error}"
        ) from error
    return (
        pd,
        genome,
        dna_client,
        dna_model,
        variant_scorers,
        importlib.metadata.version("alphagenome"),
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_input.tsv")
    parser.add_argument("--out", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_scores_summary.json")
    parser.add_argument("--status-out", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_coverage.json")
    parser.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY", help="Environment variable containing the API key.")
    parser.add_argument(
        "--sequence-length",
        choices=["16KB", "100KB", "500KB", "1MB"],
        default=SCORING_SEQUENCE_LENGTH,
    )
    parser.add_argument(
        "--model-version",
        choices=["ALL_FOLDS", "FOLD_0", "FOLD_1", "FOLD_2", "FOLD_3"],
        default=SCORING_MODEL_VERSION,
    )
    parser.add_argument(
        "--top-tracks-per-modality",
        type=int,
        default=0,
        help=(
            "Optional cap per scorer after target-scope filtering. "
            "The default 0 preserves every gastric/immune target track."
        ),
    )
    parser.add_argument("--snp-id", action="append", default=[], help="Score only this rsID; repeat option for multiple SNPs.")
    parser.add_argument("--max-variants", type=int, default=0, help="Score at most this many unscored SNPs (0 means all).")
    parser.add_argument("--workers", type=int, default=5, help="Concurrent AlphaGenome requests. The official SDK uses 5 by default.")
    parser.add_argument("--checkpoint-every", type=int, default=5, help="Atomically publish results after this many completed SNPs.")
    parser.add_argument("--overwrite", action="store_true", help="Discard all existing summaries before starting.")
    parser.add_argument("--retries", type=int, default=3, help="Retry a transient API failure this many times before recording a failure.")
    return parser.parse_args()


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(f"{path.suffix}.tmp")
    temporary_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    temporary_path.replace(path)


def read_json_list(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []
    return value if isinstance(value, list) else []


def read_json_object(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return value if isinstance(value, dict) else {}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_rows(path: Path, pd) -> list[dict]:
    if not path.exists():
        raise SystemExit(f"Input file not found: {path}. Run build_gims_portal_data.py first.")
    frame = pd.read_csv(path, sep="\t", dtype=str).dropna(subset=["chromosome", "position", "reference_bases", "alternate_bases"])
    return frame.to_dict("records")


def metadata_value(selected, pd, *names: str) -> str | None:
    for name in names:
        item = selected.get(name)
        if pd.notna(item):
            return str(item)
    return None


def summarise_scores(
    snp_id: str,
    alternate: str,
    tidy,
    pd,
    top_tracks: int,
) -> list[dict]:
    # The raw score is the scientific quantity shown in the portal. AlphaGenome
    # may also return a quantile layer; that normalized value is only used to
    # select the strongest tracks of a scorer, never as the displayed effect.
    ranking_column = "quantile_score" if "quantile_score" in tidy and tidy["quantile_score"].notna().any() else "raw_score"
    raw_numeric = pd.to_numeric(tidy["raw_score"], errors="coerce")
    ranking_numeric = pd.to_numeric(tidy[ranking_column], errors="coerce")
    working = tidy.assign(_raw_score=raw_numeric, _ranking_score=ranking_numeric).dropna(subset=["_raw_score", "_ranking_score"])
    if working.empty:
        return []

    summaries: list[dict] = []
    for (output_type, variant_scorer), group in working.groupby(["output_type", "variant_scorer"], dropna=False):
        ranked = group.assign(_absolute_score=group["_ranking_score"].abs()).sort_values("_absolute_score", ascending=False)
        seen_tracks = set()
        distinct_rank = 0
        for _, selected in ranked.iterrows():
            selected_scopes = target_scopes(selected)
            if not selected_scopes:
                continue
            track = metadata_value(selected, pd, "track_name", "assay_title", "title")
            biosample = metadata_value(selected, pd, "biosample_name")
            biosample_type = metadata_value(selected, pd, "biosample_type")
            ontology_curie = metadata_value(selected, pd, "ontology_curie")
            transcription_factor = metadata_value(selected, pd, "transcription_factor")
            histone_mark = metadata_value(selected, pd, "histone_mark")
            gene_id = metadata_value(selected, pd, "gene_id")
            gene_name = metadata_value(selected, pd, "gene_name")
            gene_strand = metadata_value(selected, pd, "gene_strand")
            track_strand = metadata_value(selected, pd, "track_strand")
            signature = (
                track,
                biosample,
                transcription_factor,
                histone_mark,
                gene_id,
                gene_strand,
                track_strand,
            )
            if signature in seen_tracks:
                continue
            seen_tracks.add(signature)
            distinct_rank += 1
            if top_tracks and distinct_rank > top_tracks:
                break
            summaries.append(
                {
                    "snpId": snp_id,
                    "alternate": alternate,
                    "outputType": str(output_type),
                    "variantScorer": str(variant_scorer),
                    "rankInOutput": distinct_rank,
                    "selection": "target_scope",
                    "scopes": selected_scopes,
                    "scoreType": "raw_score",
                    "score": float(selected["_raw_score"]),
                    "rawScore": float(selected["_raw_score"]),
                    "rankingScoreType": ranking_column,
                    "rankingScore": float(selected["_ranking_score"]),
                    "rankingAbsScore": abs(float(selected["_ranking_score"])),
                    "track": track,
                    "biosample": biosample,
                    "biosampleType": biosample_type,
                    "ontologyCurie": ontology_curie,
                    "assayTitle": metadata_value(
                        selected, pd, "Assay title", "assay_title"
                    ),
                    "gtexTissue": metadata_value(selected, pd, "gtex_tissue"),
                    "geneId": gene_id,
                    "geneName": gene_name,
                    "geneType": metadata_value(selected, pd, "gene_type"),
                    "geneStrand": gene_strand,
                    "trackStrand": track_strand,
                    "dataSource": metadata_value(selected, pd, "data_source"),
                    "transcriptionFactor": transcription_factor,
                    "histoneMark": histone_mark,
                }
            )
    return summaries


def coverage_manifest(
    rows: list[dict],
    results: list[dict],
    failed_ids: set[str],
    args,
    *,
    input_sha256: str,
    sdk_version: str,
    sequence_length_bases: int,
    scorer_names: list[str],
) -> dict:
    input_ids = {row["snp_id"] for row in rows}
    scored_ids = {row["snpId"] for row in results}
    scored_alleles = {
        (row["snpId"], row.get("alternate")) for row in results
    }
    return {
        "schemaVersion": 2,
        "source": "AlphaGenome SDK RECOMMENDED_VARIANT_SCORERS",
        "sdkVersion": sdk_version,
        "modelVersion": args.model_version,
        "organism": "HOMO_SAPIENS",
        "sequenceLength": args.sequence_length,
        "sequenceLengthBases": sequence_length_bases,
        "inputSha256": input_sha256,
        "scopeFilterVersion": SCORING_SCOPE_FILTER_VERSION,
        "targetScopes": list(TARGET_SCOPE_NAMES),
        "retainsOnlyTargetScopes": True,
        "topTracksPerModality": args.top_tracks_per_modality,
        "variantScorers": scorer_names,
        "totalInputSnps": len(input_ids),
        "totalInputAlleles": sum(
            len([alt for alt in str(row["alternate_bases"]).split(",") if alt])
            for row in rows
        ),
        "scoredSnps": len(scored_ids),
        "scoredSnpIds": sorted(scored_ids),
        "scoredAlleles": len(scored_alleles),
        "remainingSnps": max(0, len(input_ids - scored_ids)),
        "failedSnps": len(failed_ids),
        "summaryRows": len(results),
        "scopeRows": {
            scope: sum(scope in row.get("scopes", []) for row in results)
            for scope in (*TARGET_SCOPE_NAMES, "shared")
        },
        "modalities": sorted({row["outputType"] for row in results}),
        "completed": len(input_ids - scored_ids) == 0 and not failed_ids,
        "updatedAt": datetime.now(UTC).isoformat(),
    }


def compatible_cache(
    manifest: dict,
    args,
    *,
    input_sha256: str,
    sdk_version: str,
) -> bool:
    return (
        manifest.get("schemaVersion") == 2
        and manifest.get("sequenceLength") == args.sequence_length
        and manifest.get("modelVersion") == args.model_version
        and manifest.get("sdkVersion") == sdk_version
        and manifest.get("inputSha256") == input_sha256
        and manifest.get("scopeFilterVersion") == SCORING_SCOPE_FILTER_VERSION
        and manifest.get("retainsOnlyTargetScopes") is True
    )


def score_variant(model, dna_client, genome, variant_scorers, row: dict, alternate: str, args):
    variant = genome.Variant(
        chromosome=row["chromosome"],
        position=int(row["position"]),
        reference_bases=row["reference_bases"],
        alternate_bases=alternate,
    )
    interval = variant.reference_interval.resize(dna_client.SUPPORTED_SEQUENCE_LENGTHS[f"SEQUENCE_LENGTH_{args.sequence_length}"])
    for attempt in range(args.retries):
        try:
            return model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
            )
        except Exception:
            if attempt == args.retries - 1:
                raise
            time.sleep(2**attempt)


def score_row(model, dna_client, genome, variant_scorers, row: dict, args, pd) -> list[dict]:
    row_summaries = []
    for alternate in (alt for alt in str(row["alternate_bases"]).split(",") if alt):
        scores = score_variant(model, dna_client, genome, variant_scorers, row, alternate, args)
        row_summaries.extend(
            summarise_scores(
                row["snp_id"],
                alternate,
                variant_scorers.tidy_scores(scores),
                pd,
                args.top_tracks_per_modality,
            )
        )
    return row_summaries


def main() -> int:
    args = parse_args()
    if args.top_tracks_per_modality < 0:
        raise SystemExit("--top-tracks-per-modality cannot be negative")
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    if args.checkpoint_every < 1:
        raise SystemExit("--checkpoint-every must be at least 1")
    api_key = os.environ.get(args.api_key_env)
    if not api_key:
        raise SystemExit(f"{args.api_key_env} is not set. AlphaGenome prediction was not started.")

    (
        pd,
        genome,
        dna_client,
        dna_model,
        variant_scorers,
        sdk_version,
    ) = required_packages()
    rows = load_rows(args.input, pd)
    input_sha256 = file_sha256(args.input)
    sequence_length_bases = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f"SEQUENCE_LENGTH_{args.sequence_length}"
    ]
    input_alternates = {
        row["snp_id"]: {
            alternate
            for alternate in str(row["alternate_bases"]).split(",")
            if alternate
        }
        for row in rows
    }
    cache_is_compatible = compatible_cache(
        read_json_object(args.status_out),
        args,
        input_sha256=input_sha256,
        sdk_version=sdk_version,
    )
    if not args.overwrite and args.out.exists() and not cache_is_compatible:
        print(
            "Existing AlphaGenome cache uses different alleles, scopes, model, "
            "SDK or sequence length; starting a fresh reproducible batch.",
            flush=True,
        )
    results = (
        []
        if args.overwrite or not cache_is_compatible
        else [
            result
            for result in read_json_list(args.out)
            if result.get("snpId") in input_alternates
            and result.get("alternate") in input_alternates[result["snpId"]]
        ]
    )
    scored_alternates: dict[str, set[str]] = {}
    for result in results:
        scored_alternates.setdefault(result.get("snpId"), set()).add(
            result.get("alternate")
        )
    complete_ids = {
        snp_id
        for snp_id, alternates in input_alternates.items()
        if alternates and alternates.issubset(scored_alternates.get(snp_id, set()))
    }
    requested_ids = set(args.snp_id)
    to_score = [
        row
        for row in rows
        if row["snp_id"] not in complete_ids
        and (not requested_ids or row["snp_id"] in requested_ids)
    ]
    if args.max_variants:
        to_score = to_score[: args.max_variants]

    model = dna_client.create(
        api_key,
        model_version=dna_model.ModelVersion[args.model_version],
    )
    scorer_names = [
        str(scorer)
        for scorer in variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()
    ]
    failed_ids: set[str] = set()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(score_row, model, dna_client, genome, variant_scorers, row, args, pd): row for row in to_score}
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            row = futures[future]
            try:
                row_summaries = future.result()
                if row_summaries:
                    results = [
                        result
                        for result in results
                        if result.get("snpId") != row["snp_id"]
                    ]
                    results.extend(row_summaries)
                print(f"Scored {index}/{len(to_score)}: {row['snp_id']} ({len(row_summaries)} track summaries)", flush=True)
            except Exception:
                failed_ids.add(row["snp_id"])
                print(f"Skipped {index}/{len(to_score)}: {row['snp_id']} after {args.retries} attempts", flush=True)
            if index % args.checkpoint_every == 0 or index == len(to_score):
                atomic_json(args.out, results)
                atomic_json(
                    args.status_out,
                    coverage_manifest(
                        rows,
                        results,
                        failed_ids,
                        args,
                        input_sha256=input_sha256,
                        sdk_version=sdk_version,
                        sequence_length_bases=sequence_length_bases,
                        scorer_names=scorer_names,
                    ),
                )

    manifest = coverage_manifest(
        rows,
        results,
        failed_ids,
        args,
        input_sha256=input_sha256,
        sdk_version=sdk_version,
        sequence_length_bases=sequence_length_bases,
        scorer_names=scorer_names,
    )
    atomic_json(args.out, results)
    atomic_json(args.status_out, manifest)
    print(f"Published {manifest['scoredSnps']}/{manifest['totalInputSnps']} SNPs and {manifest['summaryRows']} AlphaGenome summaries")
    return 0 if not failed_ids else 1


if __name__ == "__main__":
    raise SystemExit(main())
