"""Individual-level conditional analyses for GIM locus definition."""

from __future__ import annotations

import math
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .io import as_float, iter_table, table_header
from .parameters import GIMParameters
from .plink import common_filters, read_fam_ids, run, write_ids
from .progress import progress


def _nonempty(value: object) -> bool:
    return value is not None and str(value).strip().upper() not in {"", "NA", "NAN", "-9"}


def _validated_keep(
    *, phenotypes: str | Path, covariates: str | Path, bfile: str | Path, scratch: Path
) -> tuple[Path, list[str]]:
    """Build the only persistent analysis temporary: a compact sample keep file."""

    for label, path in (("Phenotype", Path(phenotypes)), ("Covariate", Path(covariates))):
        suffixes = [suffix.lower() for suffix in path.suffixes]
        if ".csv" in suffixes or suffixes[-1:] == [".gz"]:
            raise ValueError(
                f"{label} input must be an uncompressed, headered tab-separated file for PLINK2; got {path.name}."
            )
    phenotype_columns, covariate_columns = table_header(phenotypes), table_header(covariates)
    for label, columns in (("phenotype", phenotype_columns), ("covariate", covariate_columns)):
        if not {"FID", "IID"}.issubset(columns):
            raise ValueError(f"{label} input requires FID and IID columns.")
    if len(phenotype_columns) <= 2:
        raise ValueError("Phenotype input must contain at least one metabolite column.")
    if len(covariate_columns) <= 2:
        raise ValueError("Covariate input must contain at least one covariate column.")

    phenotype_ids: list[tuple[str, str]] = []
    for row_number, row in enumerate(iter_table(phenotypes), start=2):
        identifier = (row["FID"], row["IID"])
        if not all(_nonempty(value) for value in identifier):
            raise ValueError(f"Phenotype input has an empty FID/IID at row {row_number}.")
        for column in phenotype_columns[2:]:
            value = row.get(column)
            parsed = as_float(value) if _nonempty(value) else None
            if _nonempty(value) and (parsed is None or not math.isfinite(parsed)):
                raise ValueError(f"Phenotype column {column!r} has a nonnumeric value at row {row_number}.")
        phenotype_ids.append(identifier)
    if len(phenotype_ids) != len(set(phenotype_ids)):
        raise ValueError("Phenotype input has duplicated FID/IID rows.")
    seen_covariate_ids: set[tuple[str, str]] = set()
    covariate_ids: set[tuple[str, str]] = set()
    for row_number, row in enumerate(iter_table(covariates), start=2):
        identifier = (row["FID"], row["IID"])
        if not all(_nonempty(value) for value in identifier):
            raise ValueError(f"Covariate input has an empty FID/IID at row {row_number}.")
        if identifier in seen_covariate_ids:
            raise ValueError("Covariate input has duplicated FID/IID rows.")
        seen_covariate_ids.add(identifier)
        if all(_nonempty(row.get(column)) for column in covariate_columns[2:]):
            for column in covariate_columns[2:]:
                parsed = as_float(row.get(column))
                if parsed is None or not math.isfinite(parsed):
                    raise ValueError(f"Covariate column {column!r} has a nonnumeric value at row {row_number}.")
            covariate_ids.add(identifier)
    common = sorted(set(phenotype_ids).intersection(covariate_ids).intersection(read_fam_ids(bfile)))
    if len(common) < 3:
        raise ValueError("Fewer than three genotype-matched samples with complete covariates remain.")
    keep_file = scratch / "analysis_samples.keep"
    write_ids(keep_file, common)
    return keep_file, phenotype_columns[2:]


def _read_glm(path: Path) -> list[dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        header = handle.readline().strip().split()
        output: list[dict[str, object]] = []
        for line in handle:
            values = line.strip().split()
            if len(values) != len(header):
                continue
            row = dict(zip(header, values))
            if row.get("TEST") != "ADD" or row.get("ERRCODE", ".") not in {".", ""}:
                continue
            p_value = as_float(row.get("P"))
            if p_value is None or not 0 <= p_value <= 1:
                continue
            output.append(
                {
                    "snp_id": row.get("ID", ""),
                    "p": p_value,
                    "beta": as_float(row.get("BETA")),
                    "se": as_float(row.get("SE")),
                    "n": as_float(row.get("OBS_CT")),
                }
            )
    return output


class _PLINKConditionalRunner:
    """Run small locus-level conditional models and delete each report at once."""

    def __init__(
        self,
        *,
        bfile: Path,
        plink2: Path,
        phenotypes: Path,
        covariates: Path,
        keep_file: Path,
        parameters: GIMParameters,
        scratch: Path,
        verbose: bool,
    ) -> None:
        self.bfile = bfile
        self.plink2 = plink2
        self.phenotypes = phenotypes
        self.covariates = covariates
        self.keep_file = keep_file
        self.parameters = parameters
        self.scratch = scratch
        self.verbose = verbose
        self.index = 0

    def scan(
        self,
        *,
        chromosome: str,
        start: int,
        end: int,
        metabolites: Sequence[str],
        candidates: Sequence[str] | None,
        conditioned_on: Sequence[str],
        label: str,
    ) -> dict[str, list[dict[str, object]]]:
        output: dict[str, list[dict[str, object]]] = {}
        # One phenotype per default batch bounds temporary PLINK reports even
        # for large multi-metabolite loci. Changing the batch size only affects
        # scheduling and peak disk use, not fitted models or GIM membership.
        batch_size = self.parameters.metabolite_batch_size
        for batch_start in range(0, len(metabolites), batch_size):
            batch = list(metabolites[batch_start : batch_start + batch_size])
            self.index += 1
            prefix = self.scratch / f"glm_{self.index:07d}"
            temporary: list[Path] = []
            arguments: list[str | Path] = [
                "--bfile", self.bfile,
                "--keep", self.keep_file,
                "--pheno", self.phenotypes,
                "--pheno-name", *batch,
                "--covar", self.covariates,
                "--covar-variance-standardize",
                "--chr", chromosome,
                "--from-bp", str(start),
                "--to-bp", str(end),
                *common_filters(mac_min=self.parameters.mac_min, geno_missing_max=self.parameters.geno_missing_max, chromosome=None),
            ]
            if candidates is not None:
                variant_ids = sorted(set(candidates).union(conditioned_on))
                extract = self.scratch / f"glm_{self.index:07d}.extract"
                extract.write_text("\n".join(variant_ids) + "\n")
                temporary.append(extract)
                # PLINK2 requires this explicit acknowledgement when a small
                # candidate set is intersected with the locus coordinate range.
                arguments.extend(["--extract", extract, "--force-intersect"])
            if conditioned_on:
                condition_file = self.scratch / f"glm_{self.index:07d}.condition"
                condition_file.write_text("\n".join(conditioned_on) + "\n")
                temporary.append(condition_file)
                arguments.extend(["--condition-list", condition_file])
            arguments.extend(["--glm", "hide-covar", "omit-ref", "skip-invalid-pheno", "--threads", str(self.parameters.threads), "--out", prefix])
            run(self.plink2, arguments, context=label, quiet=not self.verbose)
            for metabolite in batch:
                matches = list(self.scratch.glob(f"{prefix.name}.{metabolite}.glm.linear"))
                output[metabolite] = _read_glm(matches[0]) if len(matches) == 1 else []
                for report in matches:
                    report.unlink(missing_ok=True)
            for artifact in [*temporary, prefix.with_suffix(".log")]:
                artifact.unlink(missing_ok=True)
        return output


def _top_association(rows: Iterable[Mapping[str, object]], excluded: set[str]) -> dict[str, object] | None:
    viable = [dict(row) for row in rows if row["snp_id"] not in excluded and row.get("p") is not None]
    return min(viable, key=lambda row: float(row["p"])) if viable else None


def _forward_select(
    runner: _PLINKConditionalRunner,
    *,
    region: Mapping[str, object],
    metabolite: str,
) -> list[dict[str, object]]:
    selected: list[str] = []
    trace: list[dict[str, object]] = []
    while True:
        results = runner.scan(
            chromosome=str(region["chromosome"]), start=int(region["start"]), end=int(region["end"]),
            metabolites=[metabolite], candidates=None, conditioned_on=selected,
            label=f"forward conditional analysis for {region['region_id']}/{metabolite}",
        )[metabolite]
        best = _top_association(results, set(selected))
        if best is None or float(best["p"]) > runner.parameters.conditional_p:
            return trace
        selected.append(str(best["snp_id"]))
        trace.append(
            {
                "region_id": region["region_id"], "metabolite": metabolite,
                "forward_order": len(selected), "snp_id": best["snp_id"],
                "beta": best["beta"], "se": best["se"], "p": best["p"], "n": best["n"],
            }
        )


def _full_model_prune(
    runner: _PLINKConditionalRunner,
    *,
    region: Mapping[str, object],
    metabolite: str,
    forward: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    selected = [str(row["snp_id"]) for row in forward]
    retained: list[dict[str, object]] = []
    for snp_id in selected:
        association = runner.scan(
            chromosome=str(region["chromosome"]), start=int(region["start"]), end=int(region["end"]),
            metabolites=[metabolite], candidates=[snp_id], conditioned_on=[item for item in selected if item != snp_id],
            label=f"full conditional model for {region['region_id']}/{metabolite}/{snp_id}",
        )[metabolite]
        match = next((row for row in association if row["snp_id"] == snp_id), None)
        if match is not None and float(match["p"]) <= runner.parameters.conditional_p:
            retained.append({
                "region_id": region["region_id"], "metabolite": metabolite, "snp_id": snp_id,
                "selection_source": "forward_then_full_model", **match,
            })
    if not retained and len(selected) == 1 and runner.parameters.force_single_forward_lead:
        retained.append({
            "region_id": region["region_id"], "metabolite": metabolite, "snp_id": selected[0],
            "selection_source": "forced_single_forward_lead", "beta": "", "se": "", "p": "", "n": "",
        })
    return retained


def _build_matrix(
    runner: _PLINKConditionalRunner,
    *,
    region: Mapping[str, object],
    metabolites: Sequence[str],
    independent: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Iteratively construct the ordered V_R × M_R condition matrix."""

    remaining = sorted({str(row["snp_id"]) for row in independent})
    selected: list[str] = []
    matrix, markers = [], []
    while remaining:
        scan = runner.scan(
            chromosome=str(region["chromosome"]), start=int(region["start"]), end=int(region["end"]),
            metabolites=metabolites, candidates=remaining, conditioned_on=selected,
            label=f"V_R x M_R conditional matrix for {region['region_id']}, order {len(selected) + 1}",
        )
        candidates = []
        for metabolite, rows in scan.items():
            for row in rows:
                if row["snp_id"] in remaining:
                    candidates.append({**row, "metabolite": metabolite})
        best = _top_association(candidates, set())
        if best is None or float(best["p"]) > runner.parameters.conditional_p:
            break
        snp_id = str(best["snp_id"])
        order = len(selected) + 1
        markers.append(
            {
                "region_id": region["region_id"], "marker_order": order, "snp_id": snp_id,
                "trigger_metabolite": best["metabolite"], "trigger_beta": best["beta"],
                "trigger_se": best["se"], "trigger_p": best["p"], "conditioned_on_n": len(selected),
            }
        )
        for metabolite in metabolites:
            match = next((row for row in scan.get(metabolite, []) if row["snp_id"] == snp_id), None)
            matrix.append(
                {
                    "region_id": region["region_id"], "marker_order": order, "snp_id": snp_id,
                    "metabolite": metabolite, "beta": match["beta"] if match else "", "se": match["se"] if match else "",
                    "p": match["p"] if match else "", "n": match["n"] if match else "",
                    "conditioned_on_n": len(selected), "testable": bool(match),
                }
            )
        selected.append(snp_id)
        remaining.remove(snp_id)
    return markers, matrix


def run_conditional_analysis(
    *,
    regions: Sequence[Mapping[str, object]],
    region_metabolites: Sequence[Mapping[str, object]],
    bfile: str | Path,
    plink2: str | Path,
    phenotypes: str | Path,
    covariates: str | Path,
    parameters: GIMParameters,
    only_regions: set[str] | None = None,
    verbose: bool = False,
) -> dict[str, list[dict[str, object]]]:
    """Run all conditional stages after regions have been defined."""

    bfile, plink2, phenotypes, covariates = Path(bfile), Path(plink2), Path(phenotypes), Path(covariates)
    if not phenotypes.is_file() or not covariates.is_file():
        raise FileNotFoundError("Phenotype or covariate input does not exist.")
    region_members: dict[str, list[str]] = defaultdict(list)
    for item in region_metabolites:
        region_id, metabolite = str(item["region_id"]), str(item["metabolite"])
        if metabolite not in region_members[region_id]:
            region_members[region_id].append(metabolite)
    all_forward, all_independent, all_markers, all_matrix, summaries = [], [], [], [], []
    with tempfile.TemporaryDirectory(prefix="gimforge_conditional_") as directory:
        scratch = Path(directory)
        progress("Validating phenotype, covariate, and genotype sample identifiers")
        keep_file, available_metabolites = _validated_keep(phenotypes=phenotypes, covariates=covariates, bfile=bfile, scratch=scratch)
        runner = _PLINKConditionalRunner(
            bfile=bfile, plink2=plink2, phenotypes=phenotypes, covariates=covariates, keep_file=keep_file,
            parameters=parameters, scratch=scratch, verbose=verbose,
        )
        selected_regions = [
            region
            for region in regions
            if only_regions is None or str(region["region_id"]) in only_regions
        ]
        for region_index, region in enumerate(selected_regions, start=1):
            region_id = str(region["region_id"])
            metabolites = sorted(set(region_members.get(region_id, [])))
            progress(
                f"Conditional region {region_index}/{len(selected_regions)}: {region_id} "
                f"({len(metabolites):,} traits)"
            )
            missing = sorted(set(metabolites).difference(available_metabolites))
            if missing:
                raise ValueError(f"{region_id}: mGWAS metabolites absent from phenotype table: {', '.join(missing)}")
            forward, independent = [], []
            for metabolite_index, metabolite in enumerate(metabolites, start=1):
                progress(
                    f"{region_id}: forward/full model trait "
                    f"{metabolite_index}/{len(metabolites)} ({metabolite})"
                )
                trace = _forward_select(runner, region=region, metabolite=metabolite)
                forward.extend(trace)
                independent.extend(_full_model_prune(runner, region=region, metabolite=metabolite, forward=trace))
            progress(f"{region_id}: constructing ordered V_R x M_R conditional matrix")
            markers, matrix = _build_matrix(runner, region=region, metabolites=metabolites, independent=independent)
            all_forward.extend(forward)
            all_independent.extend(independent)
            all_markers.extend(markers)
            all_matrix.extend(matrix)
            summaries.append(
                {
                    "region_id": region_id, "chromosome": region["chromosome"], "start": region["start"], "end": region["end"],
                    "n_metabolites": len(metabolites), "n_forward_signals": len(forward),
                    "n_independent_signals": len(independent), "n_matrix_markers": len(markers),
                }
            )
    return {
        "forward_trace": all_forward, "independent_signals": all_independent, "matrix_markers": all_markers,
        "matrix_out": all_matrix, "region_summary": summaries,
    }
