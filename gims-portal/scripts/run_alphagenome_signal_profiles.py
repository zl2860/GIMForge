#!/usr/bin/env python3
"""Publish actual AlphaGenome REF/ALT local signal tracks for portal viewing.

This is deliberately separate from the official 1 MB scalar scoring job. It
calls ``predict_variant`` for three explicit target entities (stomach, primary
T cell, primary B cell), keeps ATAC and RNA-seq outputs, and mean-pools a 16-kb
central visualization profile to display-ready bins. It never creates a
fabricated curve from a scalar score.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import importlib.metadata
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from alphagenome_scopes import SCORING_MODEL_VERSION


PORTAL_ROOT = Path(__file__).resolve().parents[1]
SIGNAL_PROFILE_VERSION = 2
SCOPES = {
    "gastric_tissue": {"label": "Stomach tissue", "ontology": "UBERON:0000945"},
    "tcell": {"label": "Primary T cell", "ontology": "CL:0000084"},
    "bcell": {"label": "Primary B cell", "ontology": "CL:0000236"},
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_input.tsv")
    parser.add_argument("--out-dir", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_signals")
    parser.add_argument("--status-out", type=Path, default=PORTAL_ROOT / "public" / "data" / "alphagenome_signal_coverage.json")
    parser.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY")
    parser.add_argument(
        "--model-version",
        choices=["ALL_FOLDS", "FOLD_0", "FOLD_1", "FOLD_2", "FOLD_3"],
        default=SCORING_MODEL_VERSION,
    )
    parser.add_argument("--snp-id", action="append", default=[])
    parser.add_argument("--max-variants", type=int, default=0)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--bin-size", type=int, default=128, help="Mean-pooling bin size, in bp, for saved profiles.")
    parser.add_argument("--max-tracks-per-modality", type=int, default=2)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


def load_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def track_label(metadata: dict) -> str:
    return str(metadata.get("name") or metadata.get("Assay title") or "AlphaGenome track")


def select_tracks(reference, alternate, bin_size: int, limit: int) -> list[dict]:
    values_ref = np.asarray(reference.values, dtype=np.float32)
    values_alt = np.asarray(alternate.values, dtype=np.float32)
    if values_ref.ndim != 2 or values_alt.shape != values_ref.shape:
        return []
    centre = values_ref.shape[0] // 2
    radius = min(250, centre)
    local_delta = np.mean(np.abs(values_alt[centre - radius: centre + radius + 1] - values_ref[centre - radius: centre + radius + 1]), axis=0)
    selected = np.argsort(local_delta)[::-1][:limit]
    records = []
    for track_index in selected:
        metadata = reference.metadata.iloc[int(track_index)].dropna().to_dict()
        pooled_ref = values_ref[:, track_index].reshape(-1, bin_size).mean(axis=1)
        pooled_alt = values_alt[:, track_index].reshape(-1, bin_size).mean(axis=1)
        records.append(
            {
                "name": track_label(metadata),
                "biosample": metadata.get("biosample_name"),
                "biosampleType": metadata.get("biosample_type"),
                "assay": metadata.get("Assay title"),
                "ontology": metadata.get("ontology_curie"),
                "ref": np.round(pooled_ref, 6).tolist(),
                "alt": np.round(pooled_alt, 6).tolist(),
                "centralMeanAbsoluteDelta": float(local_delta[track_index]),
            }
        )
    return records


def make_profile(model, row: dict, args, genome, dna_client) -> dict:
    variant = genome.Variant(
        chromosome=row["chromosome"],
        position=int(row["position"]),
        reference_bases=row["reference_bases"],
        alternate_bases=row["alternate_bases"].split(",")[0],
    )
    interval = variant.reference_interval.resize(dna_client.SUPPORTED_SEQUENCE_LENGTHS["SEQUENCE_LENGTH_16KB"])
    modalities = (dna_client.OutputType.ATAC, dna_client.OutputType.RNA_SEQ)
    scopes = {}
    for scope_id, scope in SCOPES.items():
        output = model.predict_variant(
            interval=interval,
            variant=variant,
            requested_outputs=modalities,
            ontology_terms=[scope["ontology"]],
        )
        tracks = []
        for modality in modalities:
            reference = output.reference.get(modality)
            alternate = output.alternate.get(modality)
            if reference is None or alternate is None or reference.values.shape[0] % args.bin_size:
                continue
            for track in select_tracks(reference, alternate, args.bin_size, args.max_tracks_per_modality):
                track["modality"] = modality.name
                tracks.append(track)
        scopes[scope_id] = {"label": scope["label"], "ontology": scope["ontology"], "tracks": tracks}
    return {
        "schemaVersion": SIGNAL_PROFILE_VERSION,
        "modelVersion": args.model_version,
        "snpId": row["snp_id"],
        "coordinate": {"chromosome": row["chromosome"], "position": int(row["position"]), "reference": row["reference_bases"], "alternate": row["alternate_bases"].split(",")[0]},
        "interval": {"start": interval.start, "end": interval.end, "width": interval.width, "binSize": args.bin_size, "variantOffset": int(row["position"]) - interval.start},
        "scopes": scopes,
}


def profile_matches(path: Path, row: dict, model_version: str) -> bool:
    try:
        profile = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    coordinate = profile.get("coordinate") or {}
    return (
        profile.get("schemaVersion") == SIGNAL_PROFILE_VERSION
        and profile.get("modelVersion") == model_version
        and coordinate.get("chromosome") == row["chromosome"]
        and int(coordinate.get("position") or -1) == int(row["position"])
        and coordinate.get("reference") == row["reference_bases"]
        and coordinate.get("alternate")
        == row["alternate_bases"].split(",")[0]
    )


def compatible_profile_ids(args, rows: list[dict]) -> set[str]:
    return {
        row["snp_id"]
        for row in rows
        if profile_matches(
            args.out_dir / f"{row['snp_id']}.json",
            row,
            args.model_version,
        )
    }


def write_manifest(
    args,
    total: int,
    published_ids: set[str],
    failed: int,
    *,
    input_sha256: str,
    sdk_version: str,
) -> None:
    atomic_json(args.status_out, {
        "schemaVersion": SIGNAL_PROFILE_VERSION,
        "source": "AlphaGenome predict_variant REF/ALT profiles",
        "sdkVersion": sdk_version,
        "modelVersion": args.model_version,
        "inputSha256": input_sha256,
        "entities": SCOPES,
        "modalities": ["ATAC", "RNA_SEQ"],
        "interval": "16KB central visualization profile",
        "binSize": args.bin_size,
        "totalInputSnps": total,
        "publishedSnps": len(published_ids),
        "publishedSnpIds": sorted(published_ids),
        "remainingSnps": max(total - len(published_ids), 0),
        "failedSnps": failed,
        "updatedAt": datetime.now(UTC).isoformat(),
    })


def main() -> int:
    args = parse_args()
    if args.bin_size < 1 or 16_384 % args.bin_size:
        raise SystemExit("--bin-size must be a positive divisor of 16384.")
    api_key = os.environ.get(args.api_key_env)
    if not api_key:
        raise SystemExit(f"{args.api_key_env} is not set.")
    from alphagenome.data import genome
    from alphagenome.models import dna_client, dna_model

    rows = load_rows(args.input)
    input_sha256 = file_sha256(args.input)
    sdk_version = importlib.metadata.version("alphagenome")
    requested = set(args.snp_id)
    if requested:
        rows = [row for row in rows if row["snp_id"] in requested]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    input_ids = {row["snp_id"] for row in rows}
    rows_by_id = {row["snp_id"]: row for row in rows}
    existing = set()
    if not args.overwrite:
        existing = {
            snp_id
            for snp_id, row in rows_by_id.items()
            if profile_matches(
                args.out_dir / f"{snp_id}.json",
                row,
                args.model_version,
            )
        }
    todo = [row for row in rows if row["snp_id"] not in existing]
    if args.max_variants:
        todo = todo[: args.max_variants]
    model = dna_client.create(
        api_key,
        model_version=dna_model.ModelVersion[args.model_version],
    )
    failed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(make_profile, model, row, args, genome, dna_client): row for row in todo}
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            row = futures[future]
            try:
                atomic_json(args.out_dir / f"{row['snp_id']}.json", future.result())
                print(f"Published {index}/{len(todo)} signal profiles: {row['snp_id']}", flush=True)
            except Exception as error:
                failed += 1
                print(f"Skipped {index}/{len(todo)} signal profile: {row['snp_id']} ({type(error).__name__})", flush=True)
            if index % 5 == 0 or index == len(todo):
                published_ids = compatible_profile_ids(args, rows)
                write_manifest(
                    args,
                    len(rows),
                    published_ids,
                    failed,
                    input_sha256=input_sha256,
                    sdk_version=sdk_version,
                )
    published_ids = compatible_profile_ids(args, rows)
    write_manifest(
        args,
        len(rows),
        published_ids,
        failed,
        input_sha256=input_sha256,
        sdk_version=sdk_version,
    )
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
