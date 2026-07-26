"""Reference locus construction from existing mGWAS summary statistics."""

from __future__ import annotations

import csv
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Mapping, Sequence

from .io import as_float, as_int, iter_table, table_header
from .parameters import GIMParameters
from .plink import common_filters, require_bfile, require_executable, run
from .progress import progress


class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, index: int) -> int:
        while self.parent[index] != index:
            self.parent[index] = self.parent[self.parent[index]]
            index = self.parent[index]
        return index

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _normalise_sumstats_file(
    path: str | Path,
    *,
    columns: Mapping[str, str] | None = None,
    region_p: float,
    fixed_trait: str | None = None,
) -> list[dict[str, object]]:
    """Read one long-format or single-trait summary-statistic file."""

    columns = {
        "metabolite": "metabolite",
        "chromosome": "chromosome",
        "position": "position",
        "snp_id": "snp_id",
        "p": "p",
        **dict(columns or {}),
    }
    header = table_header(path)
    if not header:
        raise ValueError("mGWAS summary statistics are empty.")
    aliases = {
        "metabolite": ("metabolite", "Metabolite", "trait", "phenotype"),
        "chromosome": ("chromosome", "CHROM", "#CHROM", "chr"),
        "position": ("position", "POS", "BP", "pos"),
        "snp_id": ("snp_id", "ID", "SNP", "rsid", "MarkerName"),
        "p": ("p", "P", "p_value", "Pvalue"),
    }
    for target, source in list(columns.items()):
        if target == "metabolite" and fixed_trait is not None:
            continue
        if source not in header:
            inferred = next((candidate for candidate in aliases[target] if candidate in header), None)
            if inferred is None:
                raise ValueError(f"mGWAS summary statistics are missing a {target} column; expected {source}.")
            columns[target] = inferred

    output: list[dict[str, object]] = []
    for row in iter_table(path):
        chromosome = str(row[columns["chromosome"]]).removeprefix("chr")
        position, p_value = as_int(row[columns["position"]]), as_float(row[columns["p"]])
        if chromosome not in {str(value) for value in range(1, 23)}:
            continue
        if position is None or p_value is None or not 0 < p_value <= region_p:
            continue
        trait = fixed_trait if fixed_trait is not None else str(row[columns["metabolite"]])
        if not trait:
            continue
        output.append(
            {
                "metabolite": trait,
                "chromosome": chromosome,
                "position": position,
                "snp_id": str(row[columns["snp_id"]]),
                "p": p_value,
            }
        )
    return output


def _deduplicate_sumstats(rows: Sequence[Mapping[str, object]], *, region_p: float) -> list[dict[str, object]]:
    if not rows:
        raise ValueError(f"No autosomal mGWAS associations have p <= {region_p:.3g}.")
    unique: dict[tuple[object, ...], dict[str, object]] = {}
    for item in rows:
        row = dict(item)
        key = tuple(row[name] for name in ("metabolite", "chromosome", "position", "snp_id"))
        if key not in unique or row["p"] < unique[key]["p"]:
            unique[key] = row
    return sorted(unique.values(), key=lambda row: (int(row["chromosome"]), row["position"], row["p"], row["snp_id"]))


def normalise_sumstats(
    path: str | Path,
    *,
    columns: Mapping[str, str] | None = None,
    region_p: float,
) -> list[dict[str, object]]:
    """Read a stacked, long-format multi-trait summary-statistic table."""

    return _deduplicate_sumstats(
        _normalise_sumstats_file(path, columns=columns, region_p=region_p),
        region_p=region_p,
    )


def normalise_sumstats_manifest(
    manifest_path: str | Path,
    *,
    columns: Mapping[str, str] | None = None,
    region_p: float,
) -> list[dict[str, object]]:
    """Read one summary-statistic file per trait without concatenating them."""

    manifest_path = Path(manifest_path)
    header = table_header(manifest_path)
    if not {"trait_id", "path"}.issubset(header):
        raise ValueError("Summary-statistic manifest requires trait_id and path columns.")
    combined: list[dict[str, object]] = []
    seen_traits: set[str] = set()
    entries = list(iter_table(manifest_path))
    for entry_index, entry in enumerate(entries, start=1):
        row_number = entry_index + 1
        trait_id = str(entry.get("trait_id", "")).strip()
        supplied_path = str(entry.get("path", "")).strip()
        if not trait_id or not supplied_path:
            raise ValueError(f"Summary-statistic manifest has an empty trait_id/path at row {row_number}.")
        if trait_id in seen_traits:
            raise ValueError(f"Summary-statistic manifest contains duplicated trait_id {trait_id!r}.")
        seen_traits.add(trait_id)
        sumstats_path = Path(supplied_path)
        if not sumstats_path.is_absolute():
            sumstats_path = manifest_path.parent / sumstats_path
        if not sumstats_path.is_file():
            raise FileNotFoundError(f"Summary-statistic file for {trait_id!r} does not exist: {sumstats_path}")
        trait_rows = _normalise_sumstats_file(
            sumstats_path,
            columns=columns,
            region_p=region_p,
            fixed_trait=trait_id,
        )
        combined.extend(trait_rows)
        progress(
            f"GWAS input {entry_index}/{len(entries)}: {trait_id} "
            f"({len(trait_rows):,} eligible associations)"
        )
    return _deduplicate_sumstats(combined, region_p=region_p)


def _read_sumstats_input(
    *,
    sumstats_path: str | Path | None,
    sumstats_manifest: str | Path | None,
    columns: Mapping[str, str] | None,
    region_p: float,
) -> list[dict[str, object]]:
    if (sumstats_path is None) == (sumstats_manifest is None):
        raise ValueError("Provide exactly one of sumstats_path or sumstats_manifest.")
    if sumstats_manifest is not None:
        return normalise_sumstats_manifest(
            sumstats_manifest,
            columns=columns,
            region_p=region_p,
        )
    return normalise_sumstats(sumstats_path, columns=columns, region_p=region_p)  # type: ignore[arg-type]


def _write_assoc(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("ID", "P"), delimiter="\t")
        writer.writeheader()
        writer.writerows({"ID": row["snp_id"], "P": row["p"]} for row in rows)


def _read_whitespace_table(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        rows = [line.split() for line in handle if line.strip()]
    if len(rows) < 2:
        return []
    return [dict(zip(rows[0], values)) for values in rows[1:]]


def _select_sentinels(
    sumstats: list[dict[str, object]], *, bfile: Path, plink2: Path, parameters: GIMParameters, scratch: Path, verbose: bool
) -> list[dict[str, object]]:
    """One r² >= 0.1 clump leader per metabolite-specific signal."""

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in sumstats:
        grouped[str(row["metabolite"])].append(row)
    selected: list[dict[str, object]] = []
    trait_groups = sorted(grouped.items())
    for index, (metabolite, rows) in enumerate(trait_groups, start=1):
        progress(
            f"Sentinel clumping {index}/{len(trait_groups)}: {metabolite} "
            f"({len(rows):,} eligible associations)"
        )
        assoc_path, prefix = scratch / f"clump_{index:05d}.tsv", scratch / f"clump_{index:05d}"
        _write_assoc(assoc_path, rows)
        arguments = [
            "--bfile", bfile,
            *common_filters(mac_min=parameters.mac_min, geno_missing_max=parameters.geno_missing_max),
            "--clump", assoc_path,
            "--clump-p1", str(parameters.sentinel_p),
            "--clump-p2", str(parameters.sentinel_p),
            "--clump-r2", str(parameters.sentinel_clump_r2),
            "--clump-kb", str(parameters.sentinel_clump_window_kb),
            "--clump-unphased",
            "--threads", str(parameters.threads),
            "--out", prefix,
        ]
        run(plink2, arguments, context=f"sentinel selection for {metabolite}", quiet=not verbose)
        clump_file = prefix.with_suffix(".clumps")
        leaders = _read_whitespace_table(clump_file)
        lookup = {str(row["snp_id"]): row for row in rows}
        for leader in leaders:
            snp_id = leader.get("ID")
            if snp_id in lookup:
                selected.append(dict(lookup[snp_id]))
        for path in (assoc_path, clump_file, prefix.with_suffix(".log")):
            path.unlink(missing_ok=True)
    if not selected:
        raise ValueError("PLINK2 did not retain any sentinel variant. Check variant IDs and genotype QC.")
    return selected


def clump_sentinels(
    sumstats_path: str | Path | None = None,
    *,
    sumstats_manifest: str | Path | None = None,
    ld_bfile: str | Path,
    plink2: str | Path,
    parameters: GIMParameters,
    columns: Mapping[str, str] | None = None,
    verbose: bool = False,
) -> list[dict[str, object]]:
    """Select metabolite-specific sentinel variants with the configured LD clump."""

    parameters.validate()
    bfile, executable = require_bfile(ld_bfile), require_executable(plink2)
    progress("Reading GWAS summary statistics and filtering sentinel candidates")
    sumstats = _read_sumstats_input(
        sumstats_path=sumstats_path,
        sumstats_manifest=sumstats_manifest,
        columns=columns,
        region_p=parameters.sentinel_p,
    )
    n_traits = len({str(row["metabolite"]) for row in sumstats})
    progress(f"Eligible input loaded: {len(sumstats):,} associations across {n_traits:,} traits")
    with tempfile.TemporaryDirectory(prefix="gimforge_clump_") as directory:
        selected = _select_sentinels(
            sumstats,
            bfile=bfile,
            plink2=executable,
            parameters=parameters,
            scratch=Path(directory),
            verbose=verbose,
        )
    progress(f"Sentinel clumping complete: {len(selected):,} trait-specific leaders")
    return selected


def _ld_summaries(
    sentinels: list[dict[str, object]], *, bfile: Path, plink2: Path, parameters: GIMParameters, scratch: Path, verbose: bool
) -> tuple[dict[str, tuple[int, int]], list[tuple[str, str]]]:
    """Calculate only sentinel LD spans and sentinel-to-sentinel merge links."""

    sentinel_ids = sorted({str(row["snp_id"]) for row in sentinels})
    spans: dict[str, list[int]] = defaultdict(list)
    merge_pairs: list[tuple[str, str]] = []
    batch_size = 100
    known_sentinels = set(sentinel_ids)
    batches = list(range(0, len(sentinel_ids), batch_size))
    for batch_index, start in enumerate(batches, start=1):
        identifiers = sentinel_ids[start : start + batch_size]
        progress(f"LD span calculation {batch_index}/{len(batches)} ({len(identifiers):,} sentinels)")
        list_path, prefix = scratch / f"ld_{batch_index:05d}.ids", scratch / f"ld_{batch_index:05d}"
        list_path.write_text("\n".join(identifiers) + "\n")
        arguments = [
            "--bfile", bfile,
            *common_filters(mac_min=parameters.mac_min, geno_missing_max=parameters.geno_missing_max),
            "--ld-snp-list", list_path,
            "--r2-unphased",
            "--ld-window", "999999999",
            "--ld-window-kb", str(parameters.ld_window_kb),
            "--ld-window-r2", str(parameters.ld_span_r2),
            "--threads", str(parameters.threads),
            "--out", prefix,
        ]
        run(plink2, arguments, context=f"LD span batch {batch_index}", quiet=not verbose)
        vcor_path = prefix.with_suffix(".vcor")
        for row in _read_whitespace_table(vcor_path):
            left, right = row.get("ID_A"), row.get("ID_B")
            r2, position = as_float(row.get("UNPHASED_R2")), as_int(row.get("POS_B"))
            if not left or not right or r2 is None or position is None or r2 < parameters.ld_span_r2:
                continue
            if left != right:
                spans[left].append(position)
            if r2 > parameters.cross_metabolite_merge_r2 and right in known_sentinels and left != right:
                merge_pairs.append((left, right))
        for path in (list_path, vcor_path, prefix.with_suffix(".log")):
            path.unlink(missing_ok=True)
    return ({snp: (min(values), max(values)) for snp, values in spans.items() if values}, merge_pairs)


def define_regions(
    sumstats_path: str | Path | None = None,
    *,
    sumstats_manifest: str | Path | None = None,
    sentinels_path: str | Path | None = None,
    ld_bfile: str | Path,
    plink2: str | Path,
    parameters: GIMParameters,
    columns: Mapping[str, str] | None = None,
    verbose: bool = False,
) -> dict[str, list[dict[str, object]]]:
    """Define candidate GIM loci from mGWAS summaries or precomputed sentinels."""

    parameters.validate()
    bfile, plink2 = require_bfile(ld_bfile), require_executable(plink2)
    supplied_inputs = (sumstats_path, sumstats_manifest, sentinels_path)
    if sum(value is not None for value in supplied_inputs) != 1:
        raise ValueError("Provide exactly one of sumstats_path, sumstats_manifest, or sentinels_path.")
    if sentinels_path is not None:
        progress("Reading precomputed trait-specific sentinels; sentinel clumping will be skipped")
        sentinels = normalise_sumstats(
            sentinels_path,
            columns=columns,
            region_p=1.0,
        )
        progress(f"Precomputed sentinel input loaded: {len(sentinels):,} trait-specific leaders")
    else:
        progress("Reading GWAS summary statistics and filtering sentinel candidates")
        sumstats = _read_sumstats_input(
            sumstats_path=sumstats_path,
            sumstats_manifest=sumstats_manifest,
            columns=columns,
            region_p=parameters.sentinel_p,
        )
        n_traits = len({str(row["metabolite"]) for row in sumstats})
        progress(f"Eligible input loaded: {len(sumstats):,} associations across {n_traits:,} traits")
        sentinels = []
    with tempfile.TemporaryDirectory(prefix="gimforge_regions_") as directory:
        scratch = Path(directory)
        if sentinels_path is None:
            sentinels = _select_sentinels(
                sumstats,
                bfile=bfile,
                plink2=plink2,
                parameters=parameters,
                scratch=scratch,
                verbose=verbose,
            )
        progress("Constructing sentinel LD spans and cross-trait LD links")
        spans, merge_pairs = _ld_summaries(sentinels, bfile=bfile, plink2=plink2, parameters=parameters, scratch=scratch, verbose=verbose)

    initial: list[dict[str, object]] = []
    for index, sentinel in enumerate(sentinels, start=1):
        snp_id, position = str(sentinel["snp_id"]), int(sentinel["position"])
        if snp_id in spans:
            start, end = min(position, spans[snp_id][0]), max(position, spans[snp_id][1])
            definition = f"ld_span_r2_ge_{parameters.ld_span_r2:g}"
        else:
            half_width_bp = parameters.no_ld_half_width_kb * 1000
            start, end, definition = max(1, position - half_width_bp), position + half_width_bp, "sentinel_fallback_window_no_ld_neighbour"
        initial.append({**sentinel, "initial_region_id": f"initial_{index:05d}", "start": start, "end": end, "initial_definition": definition})

    union_find = _UnionFind(len(initial))
    by_snp: dict[str, list[int]] = defaultdict(list)
    for index, item in enumerate(initial):
        by_snp[str(item["snp_id"])].append(index)
    for indices in by_snp.values():
        for index in indices[1:]:
            union_find.union(indices[0], index)
    for left, right in merge_pairs:
        for left_index in by_snp.get(left, []):
            for right_index in by_snp.get(right, []):
                if initial[left_index]["chromosome"] == initial[right_index]["chromosome"]:
                    union_find.union(left_index, right_index)

    grouped: dict[int, list[dict[str, object]]] = defaultdict(list)
    for index, item in enumerate(initial):
        grouped[union_find.find(index)].append(item)
    pre_regions: list[dict[str, object]] = []
    for group in grouped.values():
        chromosome = str(group[0]["chromosome"])
        start, end = min(int(item["start"]) for item in group), max(int(item["end"]) for item in group)
        pre_regions.append(
            {
                "chromosome": chromosome,
                "start_before_margin": start,
                "end_before_margin": end,
                "start": max(1, start - parameters.region_padding_kb * 1000),
                "end": end + parameters.region_padding_kb * 1000,
                "initial_region_ids": [str(item["initial_region_id"]) for item in group],
            }
        )

    merged: list[dict[str, object]] = []
    for item in sorted(pre_regions, key=lambda row: (int(row["chromosome"]), int(row["start"]))):
        if merged and item["chromosome"] == merged[-1]["chromosome"] and int(item["start"]) <= int(merged[-1]["end"]):
            previous = merged[-1]
            previous["start_before_margin"] = min(int(previous["start_before_margin"]), int(item["start_before_margin"]))
            previous["end_before_margin"] = max(int(previous["end_before_margin"]), int(item["end_before_margin"]))
            previous["start"] = min(int(previous["start"]), int(item["start"]))
            previous["end"] = max(int(previous["end"]), int(item["end"]))
            previous["initial_region_ids"].extend(item["initial_region_ids"])
        else:
            merged.append(dict(item))

    initial_to_final: dict[str, str] = {}
    regions: list[dict[str, object]] = []
    for index, item in enumerate(merged, start=1):
        region_id = f"GIMForge_region_{index:03d}"
        for initial_id in item["initial_region_ids"]:
            initial_to_final[initial_id] = region_id
        regions.append(
            {
                "region_id": region_id,
                "chromosome": item["chromosome"],
                "start": item["start"],
                "end": item["end"],
                "start_before_margin": item["start_before_margin"],
                "end_before_margin": item["end_before_margin"],
                "n_initial_regions": len(item["initial_region_ids"]),
            }
        )
    memberships = []
    for item in initial:
        memberships.append(
            {
                "region_id": initial_to_final[str(item["initial_region_id"])],
                "metabolite": item["metabolite"],
                "sentinel_id": item["snp_id"],
                "sentinel_position": item["position"],
                "sentinel_p": item["p"],
                "initial_region_id": item["initial_region_id"],
                "initial_definition": item["initial_definition"],
            }
        )
    progress(
        f"Region construction complete: {len(sentinels):,} sentinels merged into "
        f"{len(regions):,} candidate regions"
    )
    return {"sentinels": sentinels, "regions": regions, "region_metabolites": memberships}
