#!/usr/bin/env python3
"""Portable GIMForge → annotation → GIMs Portal workflow."""

from __future__ import annotations

import argparse
import csv
import getpass
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

from alphagenome_scopes import (
    SCORING_MODEL_VERSION,
    SCORING_SCOPE_FILTER_VERSION,
    SCORING_SEQUENCE_LENGTH,
)
from portal_pipeline import (
    annotate_cytobands,
    annotate_ensembl,
    build_data,
    bundle_snp_ids,
    read_bundle_manifest,
    read_tsv,
    table_path,
    value,
    write_alphagenome_input,
)


PORTAL_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PORTAL_ROOT / "portal.config.json"
DEFAULT_BUNDLE = PORTAL_ROOT / "portal-bundle"
DEFAULT_DATA_DIR = PORTAL_ROOT / "public" / "data"
ANNOTATION_FILES = (
    "ensembl_annotation_cache.json",
    "cytobands_grch38.json",
    "alphagenome_input.tsv",
    "alphagenome_scores_summary.json",
    "alphagenome_coverage.json",
    "alphagenome_index.json",
    "alphagenome_signal_coverage.json",
    "annotation_manifest.json",
)
ANNOTATION_DIRS = ("alphagenome_scores", "alphagenome_signals")
EXTRA_FILES = (
    "gim_single_cell_module_matrix.json",
    "gim_single_cell_module_scores.json",
)


def load_config(path: Path) -> tuple[dict, Path]:
    if not path.is_file():
        return {}, path.parent.resolve()
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Portal config must be a JSON object: {path}")
    return payload, path.parent.resolve()


def configured_path(
    args,
    config: dict,
    config_base: Path,
    name: str,
    *,
    default: Path | None = None,
) -> Path | None:
    candidate = getattr(args, name, None)
    if candidate is not None:
        return Path(candidate).expanduser().resolve()
    configured = config.get(name)
    if configured:
        path = Path(configured).expanduser()
        return (path if path.is_absolute() else config_base / path).resolve()
    return default.resolve() if default is not None else None


def find_table(directory: Path, stem: str) -> Path | None:
    for suffix in (".tsv.gz", ".tsv"):
        path = directory / f"{stem}{suffix}"
        if path.is_file():
            return path
    return None


def selected_table(
    args,
    config: dict,
    config_base: Path,
    results: Path,
    stem: str,
    *,
    required: bool = False,
) -> Path | None:
    override = configured_path(args, config, config_base, stem)
    path = override or find_table(results, stem)
    if required and (path is None or not path.is_file()):
        raise FileNotFoundError(
            f"Required GIMForge table is missing: {stem}.tsv(.gz) in {results}"
        )
    return path


def canonical_copy(source: Path | None, source_dir: Path, stem: str) -> Path | None:
    for suffix in (".tsv", ".tsv.gz"):
        stale = source_dir / f"{stem}{suffix}"
        if stale.exists():
            stale.unlink()
    if source is None:
        return None
    suffix = ".tsv.gz" if source.suffix == ".gz" else ".tsv"
    destination = source_dir / f"{stem}{suffix}"
    shutil.copy2(source, destination)
    return destination


def open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8")
    )


def infer_bim(manifest_path: Path | None) -> Path | None:
    if manifest_path is None or not manifest_path.is_file():
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    inputs = manifest.get("inputs") if isinstance(manifest, dict) else None
    prefix = inputs.get("analysis_bfile") if isinstance(inputs, dict) else None
    if not prefix:
        return None
    prefix_path = Path(str(prefix)).expanduser()
    candidate = (
        prefix_path
        if prefix_path.suffix == ".bim"
        else Path(f"{prefix_path}.bim")
    )
    return candidate.resolve() if candidate.is_file() else None


def write_variant_subset(
    bim_path: Path | None, snp_ids: set[str], output: Path
) -> int:
    rows = []
    if bim_path is not None:
        if not bim_path.is_file():
            raise FileNotFoundError(f"PLINK BIM file does not exist: {bim_path}")
        with bim_path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 6:
                    continue
                chromosome, snp_id, _cm, position, allele1, allele2 = fields[:6]
                if snp_id in snp_ids:
                    rows.append(
                        {
                            "snp_id": snp_id,
                            "chromosome": chromosome,
                            "position": position,
                            "allele1": allele1,
                            "allele2": allele2,
                        }
                    )
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "snp_id",
                "chromosome",
                "position",
                "allele1",
                "allele2",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_variant_table_subset(
    variant_path: Path, snp_ids: set[str], output: Path
) -> int:
    """Normalize the self-contained variant catalogue from modern GIMForge."""

    rows = []
    for row in read_tsv(variant_path):
        snp_id = value(row, "snp_id", "ID", "node_id")
        if not snp_id or str(snp_id) not in snp_ids:
            continue
        rows.append(
            {
                "snp_id": snp_id,
                "chromosome": value(row, "chromosome", "CHR", "chrom"),
                "position": value(row, "position", "BP", "pos"),
                "allele1": value(row, "allele1", "A1", "alt"),
                "allele2": value(row, "allele2", "A2", "ref"),
            }
        )
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "snp_id",
                "chromosome",
                "position",
                "allele1",
                "allele2",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def sha256_files(paths: list[Path], source_build: str) -> str:
    digest = hashlib.sha256(source_build.encode("utf-8"))
    for path in sorted(paths, key=lambda item: item.name):
        digest.update(path.name.encode("utf-8"))
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    return digest.hexdigest()[:20]


def read_json_list(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload if isinstance(payload, list) else []


def copy_seed_annotations(seed_dir: Path, annotation_dir: Path, snp_ids: set[str]) -> None:
    annotation_dir.mkdir(parents=True, exist_ok=True)
    cache_path = seed_dir / "ensembl_annotation_cache.json"
    if cache_path.is_file():
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
        filtered = {snp_id: cache[snp_id] for snp_id in snp_ids if snp_id in cache}
        (annotation_dir / cache_path.name).write_text(
            json.dumps(filtered, ensure_ascii=False), encoding="utf-8"
        )

    cytoband_path = seed_dir / "cytobands_grch38.json"
    if cytoband_path.is_file():
        shutil.copy2(cytoband_path, annotation_dir / cytoband_path.name)

    summary_path = seed_dir / "alphagenome_scores_summary.json"
    if summary_path.is_file():
        rows = [
            row
            for row in read_json_list(summary_path)
            if str(row.get("snpId")) in snp_ids
        ]
        (annotation_dir / summary_path.name).write_text(
            json.dumps(rows, ensure_ascii=False), encoding="utf-8"
        )

    for name in (
        "alphagenome_coverage.json",
        "alphagenome_signal_coverage.json",
    ):
        source = seed_dir / name
        if source.is_file():
            shutil.copy2(source, annotation_dir / name)

    for directory_name in ANNOTATION_DIRS:
        source_dir = seed_dir / directory_name
        destination_dir = annotation_dir / directory_name
        destination_dir.mkdir(parents=True, exist_ok=True)
        if not source_dir.is_dir():
            continue
        for snp_id in sorted(snp_ids):
            source = source_dir / f"{snp_id}.json"
            if source.is_file():
                shutil.copy2(source, destination_dir / source.name)


def prune_bundle_annotations(annotation_dir: Path, snp_ids: set[str]) -> None:
    """Retain resumable overlap while removing annotations for an older bundle."""

    cache_path = annotation_dir / "ensembl_annotation_cache.json"
    if cache_path.is_file():
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
        cache = {snp_id: cache[snp_id] for snp_id in snp_ids if snp_id in cache}
        cache_path.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
    summary_path = annotation_dir / "alphagenome_scores_summary.json"
    if summary_path.is_file():
        rows = [
            row
            for row in read_json_list(summary_path)
            if str(row.get("snpId")) in snp_ids
        ]
        summary_path.write_text(json.dumps(rows, ensure_ascii=False), encoding="utf-8")
    for directory_name in ANNOTATION_DIRS:
        directory = annotation_dir / directory_name
        if not directory.is_dir():
            continue
        for path in directory.glob("*.json"):
            if path.stem not in snp_ids:
                path.unlink()


def prepare_bundle(args, config: dict, config_base: Path) -> Path:
    results = configured_path(args, config, config_base, "results")
    if results is None or not results.is_dir():
        raise FileNotFoundError(
            "A completed GIMForge result directory is required. "
            "Use --results or set results in portal.config.json."
        )
    # Supplying --results means "use this run as-is"; current-project table
    # overrides in portal.config.json must not leak into an unrelated run.
    source_config = config if getattr(args, "results", None) is None else {}
    bundle = configured_path(
        args, source_config, config_base, "bundle", default=DEFAULT_BUNDLE
    )
    source_dir = bundle / "source"
    source_dir.mkdir(parents=True, exist_ok=True)

    copied: list[Path] = []
    required_stems = ("matrix_out", "edges", "members", "gim_summary")
    optional_stems = (
        "independent_signals",
        "regions",
        "region_summary",
        "region_metabolites",
    )
    for stem in required_stems + optional_stems:
        selected = selected_table(
            args,
            source_config,
            config_base,
            results,
            stem,
            required=stem in required_stems,
        )
        copied_path = canonical_copy(selected, source_dir, stem)
        if copied_path is not None:
            copied.append(copied_path)

    result_manifest = configured_path(
        args, source_config, config_base, "run_manifest"
    ) or (results / "run_manifest.json")
    copied_manifest = source_dir / "run_manifest.json"
    if result_manifest.is_file():
        shutil.copy2(result_manifest, copied_manifest)
        copied.append(copied_manifest)
    elif copied_manifest.exists():
        copied_manifest.unlink()

    members_path = table_path(source_dir, "members", required=True)
    member_rows = read_tsv(members_path)
    snp_ids = {
        str(value(row, "node_id"))
        for row in member_rows
        if str(value(row, "node_type", default="")).upper() == "SNP"
    }
    signals_path = table_path(source_dir, "independent_signals")
    if signals_path:
        snp_ids.update(
            str(snp_id)
            for row in read_tsv(signals_path)
            if (snp_id := value(row, "snp_id", "ID", "markername"))
        )
    explicit_bim = configured_path(args, source_config, config_base, "bim")
    result_variants = find_table(results, "variants")
    if explicit_bim is not None:
        variant_count = write_variant_subset(
            explicit_bim, snp_ids, source_dir / "variants.tsv"
        )
        variant_source = "explicit_bim"
    elif result_variants is not None:
        variant_count = write_variant_table_subset(
            result_variants, snp_ids, source_dir / "variants.tsv"
        )
        variant_source = "gimforge_variants_table"
    else:
        inferred_bim = infer_bim(
            result_manifest if result_manifest.is_file() else None
        )
        variant_count = write_variant_subset(
            inferred_bim, snp_ids, source_dir / "variants.tsv"
        )
        variant_source = (
            "run_manifest_analysis_bim"
            if inferred_bim is not None
            else "unavailable"
        )
    copied.append(source_dir / "variants.tsv")

    differential = configured_path(
        args, source_config, config_base, "differential_metabolomics"
    )
    differential_destination = source_dir / "differential_metabolomics.csv"
    if differential and differential.is_file():
        shutil.copy2(differential, differential_destination)
        copied.append(differential_destination)
    elif differential_destination.exists():
        differential_destination.unlink()

    extras_dir = bundle / "extras"
    extras_dir.mkdir(parents=True, exist_ok=True)
    extra_config = {
        "gim_single_cell_module_matrix.json": "single_cell_module_matrix",
        "gim_single_cell_module_scores.json": "single_cell_module_scores",
    }
    for destination_name, config_name in extra_config.items():
        source = configured_path(args, source_config, config_base, config_name)
        destination = extras_dir / destination_name
        if source and source.is_file():
            shutil.copy2(source, destination)
        elif destination.exists():
            destination.unlink()

    source_build = str(
        getattr(args, "source_build", None)
        or source_config.get("source_build")
        or "GRCh37"
    )
    bundle_id = sha256_files(copied, source_build)
    bundle_manifest = {
        "schemaVersion": 1,
        "bundleId": bundle_id,
        "createdAt": datetime.now(UTC).isoformat(),
        "sourceGenomeBuild": source_build,
        "sourceResultName": results.name,
        "nGimSnps": len(snp_ids),
        "nSourceCoordinates": variant_count,
        "variantSource": variant_source,
        "containsIndividualLevelData": False,
        "sourceTables": sorted(path.name for path in copied),
    }
    (bundle / "bundle_manifest.json").write_text(
        json.dumps(bundle_manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    prune_bundle_annotations(bundle / "annotations", snp_ids)
    reuse = configured_path(
        args, source_config, config_base, "reuse_annotations"
    )
    if reuse and reuse.is_dir():
        copy_seed_annotations(reuse, bundle / "annotations", snp_ids)

    print(
        f"Prepared portable bundle {bundle_id} with {len(snp_ids):,} GIM SNPs "
        f"at {bundle}"
    )
    if variant_count < len(snp_ids):
        print(
            f"Warning: source coordinates were found for {variant_count:,}/"
            f"{len(snp_ids):,} SNPs; Ensembl rsID annotation can still continue."
        )
    return bundle


def alpha_input_ids(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["snp_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("snp_id")
        }


def sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json_object(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


def score_cache_compatible(annotation_dir: Path, input_sha256: str | None) -> bool:
    coverage = read_json_object(annotation_dir / "alphagenome_coverage.json")
    return (
        coverage.get("schemaVersion") == 2
        and coverage.get("sequenceLength") == SCORING_SEQUENCE_LENGTH
        and coverage.get("modelVersion") == SCORING_MODEL_VERSION
        and bool(coverage.get("sdkVersion"))
        and coverage.get("inputSha256") == input_sha256
        and coverage.get("scopeFilterVersion") == SCORING_SCOPE_FILTER_VERSION
        and coverage.get("retainsOnlyTargetScopes") is True
    )


def signal_cache_ids(
    annotation_dir: Path,
    input_ids: set[str],
    input_sha256: str | None,
) -> set[str]:
    coverage = read_json_object(
        annotation_dir / "alphagenome_signal_coverage.json"
    )
    if not signal_cache_compatible(coverage, input_sha256):
        return set()
    published = coverage.get("publishedSnpIds")
    if not isinstance(published, list):
        return set()
    return {
        str(snp_id)
        for snp_id in published
        if str(snp_id) in input_ids
        and (
            annotation_dir / "alphagenome_signals" / f"{snp_id}.json"
        ).is_file()
    }


def signal_cache_compatible(
    coverage: dict,
    input_sha256: str | None,
) -> bool:
    return (
        coverage.get("schemaVersion") == 2
        and coverage.get("modelVersion") == SCORING_MODEL_VERSION
        and bool(coverage.get("sdkVersion"))
        and coverage.get("inputSha256") == input_sha256
        and isinstance(coverage.get("publishedSnpIds"), list)
    )


def existing_score_ids(path: Path, input_ids: set[str]) -> set[str]:
    return {
        str(row.get("snpId"))
        for row in read_json_list(path)
        if str(row.get("snpId")) in input_ids
    }


def score_cache_ids(
    annotation_dir: Path,
    input_ids: set[str],
    input_sha256: str | None,
) -> set[str]:
    if not score_cache_compatible(annotation_dir, input_sha256):
        return set()
    coverage = read_json_object(
        annotation_dir / "alphagenome_coverage.json"
    )
    published = coverage.get("scoredSnpIds")
    if isinstance(published, list):
        return {
            str(snp_id)
            for snp_id in published
            if str(snp_id) in input_ids
            and (
                annotation_dir / "alphagenome_scores" / f"{snp_id}.json"
            ).is_file()
        }
    scores_dir = annotation_dir / "alphagenome_scores"
    split_ids = (
        {
            path.stem
            for path in scores_dir.glob("*.json")
            if path.stem in input_ids
        }
        if scores_dir.is_dir()
        else set()
    )
    if split_ids:
        return split_ids
    # Compatibility fallback for a partial legacy cache that predates both the
    # compact ID manifest and per-SNP splitting.
    return existing_score_ids(
        annotation_dir / "alphagenome_scores_summary.json", input_ids
    )


def annotation_state(bundle: Path) -> dict:
    annotation_dir = bundle / "annotations"
    input_path = annotation_dir / "alphagenome_input.tsv"
    input_ids = alpha_input_ids(input_path)
    input_sha256 = sha256_file(input_path)
    scored_ids = score_cache_ids(
        annotation_dir, input_ids, input_sha256
    )
    signal_ids = signal_cache_ids(annotation_dir, input_ids, input_sha256)
    return {
        "input": input_ids,
        "scored": scored_ids,
        "signals": signal_ids,
        "missingScores": input_ids - scored_ids,
        "missingSignals": input_ids - signal_ids,
    }


def normalize_coverage_manifests(bundle: Path) -> None:
    """Make coverage counts describe this bundle, not a previous resumable run."""

    annotation_dir = bundle / "annotations"
    state = annotation_state(bundle)
    input_sha256 = sha256_file(annotation_dir / "alphagenome_input.tsv")
    score_compatible = score_cache_compatible(annotation_dir, input_sha256)
    coverage_path = annotation_dir / "alphagenome_coverage.json"
    coverage = (
        json.loads(coverage_path.read_text(encoding="utf-8"))
        if coverage_path.is_file()
        else {
            "schemaVersion": 1,
            "source": "AlphaGenome SDK RECOMMENDED_VARIANT_SCORERS",
        }
    )
    if not score_compatible:
        coverage = {
            "schemaVersion": 2,
            "source": "AlphaGenome SDK RECOMMENDED_VARIANT_SCORERS",
            "sequenceLength": SCORING_SEQUENCE_LENGTH,
            "modelVersion": SCORING_MODEL_VERSION,
            "scopeFilterVersion": SCORING_SCOPE_FILTER_VERSION,
            "inputSha256": input_sha256,
            "retainsOnlyTargetScopes": True,
            "cacheStatus": "pending_required_rescore",
        }
    coverage.update(
        {
            "totalInputSnps": len(state["input"]),
            "scoredSnps": len(state["scored"]),
            "scoredSnpIds": sorted(state["scored"]),
            "remainingSnps": len(state["missingScores"]),
            "summaryRows": (
                int(coverage.get("summaryRows") or 0)
                if score_compatible
                else 0
            ),
            "cacheCompatible": score_compatible,
            "requiredSequenceLength": SCORING_SEQUENCE_LENGTH,
            "requiredModelVersion": SCORING_MODEL_VERSION,
            "requiredScopeFilterVersion": SCORING_SCOPE_FILTER_VERSION,
            "completed": not state["missingScores"],
            "updatedAt": datetime.now(UTC).isoformat(),
        }
    )
    coverage_path.write_text(
        json.dumps(coverage, ensure_ascii=False), encoding="utf-8"
    )

    signal_path = annotation_dir / "alphagenome_signal_coverage.json"
    signal_coverage = (
        json.loads(signal_path.read_text(encoding="utf-8"))
        if signal_path.is_file()
        else {
            "schemaVersion": 1,
            "source": "AlphaGenome predict_variant REF/ALT profiles",
        }
    )
    signal_compatible = signal_cache_compatible(
        signal_coverage, input_sha256
    )
    if not signal_compatible:
        signal_coverage = {
            "schemaVersion": 2,
            "source": "AlphaGenome predict_variant REF/ALT profiles",
            "modelVersion": SCORING_MODEL_VERSION,
            "inputSha256": input_sha256,
            "cacheStatus": "pending_required_rescore",
        }
    signal_coverage.update(
        {
            "totalInputSnps": len(state["input"]),
            "publishedSnps": len(state["signals"]),
            "remainingSnps": len(state["missingSignals"]),
            "cacheCompatible": signal_compatible,
            "requiredModelVersion": SCORING_MODEL_VERSION,
            "completed": not state["missingSignals"],
            "updatedAt": datetime.now(UTC).isoformat(),
        }
    )
    signal_path.write_text(
        json.dumps(signal_coverage, ensure_ascii=False), encoding="utf-8"
    )


def alpha_python_path(
    args, config: dict, config_base: Path
) -> Path:
    configured = configured_path(args, config, config_base, "alpha_python")
    if configured is not None:
        return configured
    return PORTAL_ROOT / ".venv-alphagenome" / "bin" / "python"


def check_alpha_dependencies(alpha_python: Path) -> dict:
    """Validate the isolated AlphaGenome runtime before requesting a credential."""
    check = (
        "import importlib.metadata as metadata, json; "
        "import numpy, pandas; "
        "from alphagenome.data import genome; "
        "from alphagenome.models import dna_client, dna_model, variant_scorers; "
        "print(json.dumps({"
        "'alphagenome': metadata.version('alphagenome'), "
        "'pandas': pandas.__version__, "
        "'numpy': numpy.__version__"
        "}))"
    )
    result = subprocess.run(
        [str(alpha_python), "-c", check],
        cwd=PORTAL_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        details = (result.stderr or result.stdout).strip().splitlines()
        detail = details[-1] if details else "runtime import check failed"
        install = (
            f'"{alpha_python}" -m pip install '
            '"git+https://github.com/google-deepmind/alphagenome.git" '
            "pandas numpy"
        )
        raise RuntimeError(
            f"AlphaGenome annotation dependencies are missing or broken in "
            f"{alpha_python} ({detail}). Install them before entering an API key:\n"
            f"{install}"
        )
    try:
        versions = json.loads(result.stdout.strip())
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"AlphaGenome dependency check returned an unreadable response "
            f"from {alpha_python}."
        ) from error
    if not isinstance(versions, dict):
        raise RuntimeError(
            f"AlphaGenome dependency check returned an invalid response "
            f"from {alpha_python}."
        )
    return versions


def request_api_key(environment_name: str) -> tuple[str, bool]:
    existing = os.environ.get(environment_name)
    if existing:
        print(f"Using AlphaGenome credential from {environment_name}.")
        return existing, False
    if not sys.stdin.isatty():
        raise RuntimeError(
            f"AlphaGenome work remains, but no interactive terminal is available. "
            f"Set {environment_name} in the process environment."
        )
    secret = getpass.getpass(
        "AlphaGenome API key (hidden; used only by this process): "
    ).strip()
    if not secret:
        raise RuntimeError("No AlphaGenome API key was entered.")
    return secret, True


def run_process(command: list[str], *, env: dict | None = None) -> int:
    result = subprocess.run(command, cwd=PORTAL_ROOT, env=env, check=False)
    return result.returncode


def split_alpha_scores(alpha_python: Path, annotation_dir: Path) -> int:
    annotation_dir.mkdir(parents=True, exist_ok=True)
    summary = annotation_dir / "alphagenome_scores_summary.json"
    input_sha256 = sha256_file(annotation_dir / "alphagenome_input.tsv")
    if not summary.is_file() or not score_cache_compatible(
        annotation_dir, input_sha256
    ):
        (annotation_dir / "alphagenome_index.json").write_text(
            "[]", encoding="utf-8"
        )
        return 0
    return run_process(
        [
            str(alpha_python),
            str(PORTAL_ROOT / "scripts" / "split_alphagenome_portal_data.py"),
            "--input",
            str(summary),
            "--index-out",
            str(annotation_dir / "alphagenome_index.json"),
            "--scores-dir",
            str(annotation_dir / "alphagenome_scores"),
        ]
    )


def write_annotation_manifest(bundle: Path) -> dict:
    source_manifest = read_bundle_manifest(bundle)
    state = annotation_state(bundle)
    cache_path = bundle / "annotations" / "ensembl_annotation_cache.json"
    ensembl = (
        json.loads(cache_path.read_text(encoding="utf-8"))
        if cache_path.is_file()
        else {}
    )
    payload = {
        "schemaVersion": 2,
        "bundleId": source_manifest["bundleId"],
        "updatedAt": datetime.now(UTC).isoformat(),
        "nBundleSnps": len(bundle_snp_ids(bundle)),
        "nEnsemblResolvedRequests": len(ensembl),
        "cytobandsResolved": (
            bundle / "annotations" / "cytobands_grch38.json"
        ).is_file(),
        "nAlphaGenomeReady": len(state["input"]),
        "nAlphaGenomeScored": len(state["scored"]),
        "nSignalProfiles": len(state["signals"]),
        "alphaGenomeSequenceLength": SCORING_SEQUENCE_LENGTH,
        "alphaGenomeModelVersion": SCORING_MODEL_VERSION,
        "alphaGenomeScopeFilterVersion": SCORING_SCOPE_FILTER_VERSION,
        "alphaGenomeComplete": not state["missingScores"],
        "signalProfilesComplete": not state["missingSignals"],
        "credentialStored": False,
    }
    (bundle / "annotations" / "annotation_manifest.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return payload


def annotate_bundle(args, config: dict, config_base: Path, bundle: Path | None = None) -> int:
    bundle = bundle or configured_path(
        args, config, config_base, "bundle", default=DEFAULT_BUNDLE
    )
    read_bundle_manifest(bundle)
    annotation_dir = bundle / "annotations"
    annotation_dir.mkdir(parents=True, exist_ok=True)

    print("Resolving portable Ensembl/GRCh38 annotations…")
    annotate_ensembl(bundle)
    print("Resolving GRCh38 cytoband annotations…")
    annotate_cytobands(bundle)
    alpha_input = write_alphagenome_input(bundle)
    state = annotation_state(bundle)
    skip_alpha = bool(getattr(args, "skip_alphagenome", False))
    skip_signals = bool(getattr(args, "skip_signals", False))
    max_variants = int(getattr(args, "max_variants", 0) or 0)
    workers = int(getattr(args, "workers", 5) or 5)
    overwrite = bool(getattr(args, "overwrite", False))

    if skip_alpha:
        normalize_coverage_manifests(bundle)
        manifest = write_annotation_manifest(bundle)
        print(
            f"Ensembl stage complete; AlphaGenome was skipped. "
            f"{manifest['nAlphaGenomeReady']:,} variants are ready."
        )
        return 0

    needs_scores = bool(state["missingScores"]) or overwrite
    needs_signals = (bool(state["missingSignals"]) or overwrite) and not skip_signals
    if not needs_scores and not needs_signals:
        alpha_python = alpha_python_path(args, config, config_base)
        split_alpha_scores(alpha_python if alpha_python.is_file() else Path(sys.executable), annotation_dir)
        normalize_coverage_manifests(bundle)
        manifest = write_annotation_manifest(bundle)
        print(
            f"Annotations are already complete for "
            f"{manifest['nAlphaGenomeScored']:,} variants; no API key was requested."
        )
        return 0

    alpha_python = alpha_python_path(args, config, config_base)
    if not alpha_python.is_file():
        raise FileNotFoundError(
            f"AlphaGenome Python environment is missing: {alpha_python}. "
            "Create .venv-alphagenome and install the official SDK plus pandas/numpy."
        )
    versions = check_alpha_dependencies(alpha_python)
    print(
        "AlphaGenome runtime ready: "
        f"SDK {versions.get('alphagenome', 'unknown')}, "
        f"pandas {versions.get('pandas', 'unknown')}, "
        f"numpy {versions.get('numpy', 'unknown')}."
    )
    api_env_name = str(
        getattr(args, "api_key_env", None)
        or config.get("api_key_env")
        or "ALPHAGENOME_API_KEY"
    )
    api_key, prompted = request_api_key(api_env_name)
    child_env = os.environ.copy()
    child_env[api_env_name] = api_key
    status = 0
    common_limit = ["--max-variants", str(max_variants)] if max_variants else []

    if needs_scores:
        score_command = [
            str(alpha_python),
            str(PORTAL_ROOT / "scripts" / "run_alphagenome_scoring.py"),
            "--input",
            str(alpha_input),
            "--out",
            str(annotation_dir / "alphagenome_scores_summary.json"),
            "--status-out",
            str(annotation_dir / "alphagenome_coverage.json"),
            "--api-key-env",
            api_env_name,
            "--sequence-length",
            SCORING_SEQUENCE_LENGTH,
            "--model-version",
            SCORING_MODEL_VERSION,
            "--workers",
            str(workers),
            *common_limit,
        ]
        if overwrite:
            score_command.append("--overwrite")
        status = max(status, run_process(score_command, env=child_env))
        split_status = split_alpha_scores(alpha_python, annotation_dir)
        status = max(status, split_status)

    if needs_signals:
        signal_command = [
            str(alpha_python),
            str(PORTAL_ROOT / "scripts" / "run_alphagenome_signal_profiles.py"),
            "--input",
            str(alpha_input),
            "--out-dir",
            str(annotation_dir / "alphagenome_signals"),
            "--status-out",
            str(annotation_dir / "alphagenome_signal_coverage.json"),
            "--api-key-env",
            api_env_name,
            "--model-version",
            SCORING_MODEL_VERSION,
            "--workers",
            str(workers),
            *common_limit,
        ]
        if overwrite:
            signal_command.append("--overwrite")
        status = max(status, run_process(signal_command, env=child_env))

    child_env.pop(api_env_name, None)
    if prompted:
        api_key = ""
    normalize_coverage_manifests(bundle)
    manifest = write_annotation_manifest(bundle)
    print(
        f"Annotation bundle now contains {manifest['nAlphaGenomeScored']:,}/"
        f"{manifest['nAlphaGenomeReady']:,} scalar predictions and "
        f"{manifest['nSignalProfiles']:,} signal profiles."
    )
    return status


def sync_json_directory(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    source_names = {path.name for path in source.glob("*.json")} if source.is_dir() else set()
    for stale in destination.glob("*.json"):
        if stale.name not in source_names:
            stale.unlink()
    if source.is_dir():
        for path in source.glob("*.json"):
            shutil.copy2(path, destination / path.name)


def publish_bundle_assets(bundle: Path, out_dir: Path) -> None:
    annotation_dir = bundle / "annotations"
    input_sha256 = sha256_file(annotation_dir / "alphagenome_input.tsv")
    score_compatible = score_cache_compatible(annotation_dir, input_sha256)
    signal_coverage = read_json_object(
        annotation_dir / "alphagenome_signal_coverage.json"
    )
    signals_compatible = signal_cache_compatible(
        signal_coverage, input_sha256
    )
    for name in ANNOTATION_FILES:
        if name == "annotation_manifest.json":
            continue
        source = annotation_dir / name
        destination = out_dir / name
        if name == "alphagenome_scores_summary.json":
            # Keep the multi-gigabyte reproducibility table in the portable
            # annotation bundle. The browser uses the scientifically
            # equivalent per-SNP lazy files and compact index, so copying this
            # file into both public/ and dist/ would only duplicate storage.
            if destination.exists():
                destination.unlink()
            continue
        if source.is_file():
            shutil.copy2(source, destination)
        elif destination.exists() and name.startswith("alphagenome_"):
            destination.unlink()
    for directory_name in ANNOTATION_DIRS:
        compatible = (
            score_compatible
            if directory_name == "alphagenome_scores"
            else signals_compatible
        )
        sync_json_directory(
            annotation_dir / directory_name
            if compatible
            else annotation_dir / f".stale-{directory_name}",
            out_dir / directory_name,
        )
    extras_dir = bundle / "extras"
    for name in EXTRA_FILES:
        source = extras_dir / name
        destination = out_dir / name
        if source.is_file():
            shutil.copy2(source, destination)
        elif destination.exists():
            destination.unlink()
    if not (out_dir / "alphagenome_index.json").is_file():
        (out_dir / "alphagenome_index.json").write_text("[]", encoding="utf-8")


def build_bundle(args, config: dict, config_base: Path, bundle: Path | None = None) -> int:
    bundle = bundle or configured_path(
        args, config, config_base, "bundle", default=DEFAULT_BUNDLE
    )
    manifest = read_bundle_manifest(bundle)
    out_dir = configured_path(
        args, config, config_base, "out_dir", default=DEFAULT_DATA_DIR
    )
    annotation_dir = bundle / "annotations"
    alpha_python = alpha_python_path(args, config, config_base)
    split_alpha_scores(
        alpha_python if alpha_python.is_file() else Path(sys.executable),
        annotation_dir,
    )
    normalize_coverage_manifests(bundle)
    stats = build_data(bundle, out_dir)
    publish_bundle_assets(bundle, out_dir)
    annotation_manifest = write_annotation_manifest(bundle)
    complete = (
        stats["annotationStatus"] == "resolved"
        and annotation_manifest["alphaGenomeComplete"]
        and (
            annotation_manifest["signalProfilesComplete"]
            or bool(getattr(args, "allow_missing_signals", False))
        )
    )
    build_manifest = {
        "schemaVersion": 1,
        "bundleId": manifest["bundleId"],
        "builtAt": datetime.now(UTC).isoformat(),
        "complete": complete,
        "core": stats,
        "annotations": annotation_manifest,
    }
    (out_dir / "portal_build_manifest.json").write_text(
        json.dumps(build_manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    if getattr(args, "require_complete", False) and not complete:
        raise RuntimeError(
            "Portal data is partial. Finish the annotate stage or build without "
            "--require-complete."
        )
    if getattr(args, "skip_vite", False):
        print(f"Portal data prepared at {out_dir}; Vite build was skipped.")
        return 0
    status = run_process(["npm", "run", "build"])
    if status == 0:
        print(
            f"Built GIMs Portal from bundle {manifest['bundleId']} at "
            f"{PORTAL_ROOT / 'dist'}"
        )
        if not complete:
            print(
                "Note: the portal was built with partial optional annotations; "
                "pending variants remain visible as pending."
            )
    return status


def add_config_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="JSON config; defaults to portal.config.json beside package.json.",
    )


def add_bundle_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--bundle", type=Path, help="Portable bundle directory.")


def add_prepare_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--results", type=Path, help="Completed GIMForge result directory.")
    parser.add_argument("--bim", type=Path, help="Optional PLINK BIM for source coordinates.")
    parser.add_argument("--source-build", default=None, help="Source coordinate build, e.g. GRCh37.")
    for stem in (
        "matrix_out",
        "edges",
        "members",
        "gim_summary",
        "independent_signals",
        "regions",
        "region_summary",
        "region_metabolites",
    ):
        parser.add_argument(f"--{stem.replace('_', '-')}", dest=stem, type=Path)
    parser.add_argument("--run-manifest", type=Path)
    parser.add_argument("--differential-metabolomics", type=Path)
    parser.add_argument("--single-cell-module-matrix", type=Path)
    parser.add_argument("--single-cell-module-scores", type=Path)
    parser.add_argument("--reuse-annotations", type=Path)


def add_annotation_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--alpha-python", type=Path)
    parser.add_argument("--api-key-env", default=None)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--max-variants", type=int, default=0)
    parser.add_argument("--skip-alphagenome", action="store_true")
    parser.add_argument("--skip-signals", action="store_true")
    parser.add_argument("--overwrite", action="store_true")


def add_build_options(
    parser: argparse.ArgumentParser, *, include_alpha_python: bool = True
) -> None:
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--skip-vite", action="store_true")
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument("--allow-missing-signals", action="store_true")
    if include_alpha_python:
        parser.add_argument("--alpha-python", type=Path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gims-portal",
        description=(
            "Turn standard GIMForge output into a portable, annotated static portal."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare",
        help="export a portable bundle from GIMForge results (no internet required)",
    )
    add_config_option(prepare)
    add_bundle_option(prepare)
    add_prepare_options(prepare)

    annotate = subparsers.add_parser(
        "annotate",
        help="add Ensembl and AlphaGenome annotations on an internet-connected machine",
    )
    add_config_option(annotate)
    add_bundle_option(annotate)
    add_annotation_options(annotate)

    build = subparsers.add_parser(
        "build", help="build portal data and static web assets from a bundle"
    )
    add_config_option(build)
    add_bundle_option(build)
    add_build_options(build)

    all_parser = subparsers.add_parser(
        "all", help="prepare, annotate, and build on one internet-connected machine"
    )
    add_config_option(all_parser)
    add_bundle_option(all_parser)
    add_prepare_options(all_parser)
    add_annotation_options(all_parser)
    add_build_options(all_parser, include_alpha_python=False)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config_path = args.config.expanduser().resolve()
    try:
        config, config_base = load_config(config_path)
        if args.command == "prepare":
            prepare_bundle(args, config, config_base)
            return 0
        if args.command == "annotate":
            return annotate_bundle(args, config, config_base)
        if args.command == "build":
            return build_bundle(args, config, config_base)
        bundle = prepare_bundle(args, config, config_base)
        annotate_status = annotate_bundle(args, config, config_base, bundle)
        build_status = build_bundle(args, config, config_base, bundle)
        return max(annotate_status, build_status)
    except (FileNotFoundError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
