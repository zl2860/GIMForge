"""End-to-end orchestration: supplied summaries → loci → conditional matrix → GIMs."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

from .components import components_from_matrix
from .conditional import run_conditional_analysis
from .io import write_json, write_table
from .parameters import GIMParameters
from .plink import require_bfile, require_executable, version
from .progress import progress
from .regions import define_regions
from .report import write_report


_FIELDS = {
    "sentinels": ["metabolite", "chromosome", "position", "snp_id", "p"],
    "regions": ["region_id", "chromosome", "start", "end", "start_before_margin", "end_before_margin", "n_initial_regions"],
    "region_metabolites": ["region_id", "metabolite", "sentinel_id", "sentinel_position", "sentinel_p", "initial_region_id", "initial_definition"],
    "forward_trace": ["region_id", "metabolite", "forward_order", "snp_id", "beta", "se", "p", "n"],
    "independent_signals": ["region_id", "metabolite", "snp_id", "selection_source", "beta", "se", "p", "n"],
    "matrix_markers": ["region_id", "marker_order", "snp_id", "trigger_metabolite", "trigger_beta", "trigger_se", "trigger_p", "conditioned_on_n"],
    "matrix_out": ["region_id", "marker_order", "snp_id", "metabolite", "beta", "se", "p", "n", "conditioned_on_n", "testable"],
    "region_summary": ["region_id", "chromosome", "start", "end", "n_metabolites", "n_forward_signals", "n_independent_signals", "n_matrix_markers"],
    "edges": ["region_id", "marker_order", "snp_id", "metabolite", "beta", "se", "p", "n", "conditioned_on_n", "testable", "gim_id"],
    "members": ["gim_id", "region_id", "node_type", "node_id", "marker_order"],
    "gim_summary": ["gim_id", "region_id", "n_snps", "n_metabolites", "snps", "metabolites"],
}


def _write_result(output: Path, result: Mapping[str, Sequence[Mapping[str, object]]]) -> None:
    for name, fields in _FIELDS.items():
        if name in result:
            suffix = ".tsv" if name in {"regions", "region_metabolites", "sentinels", "region_summary", "gim_summary"} else ".tsv.gz"
            write_table(output / f"{name}{suffix}", result[name], fieldnames=fields)


def run_gim(
    *,
    sumstats: str | Path | None = None,
    sumstats_manifest: str | Path | None = None,
    sentinels: str | Path | None = None,
    bfile: str | Path,
    ld_bfile: str | Path,
    ancestry: str,
    ld_panel_name: str | None = None,
    phenotypes: str | Path,
    covariates: str | Path,
    output: str | Path,
    plink2: str | Path,
    parameters: GIMParameters,
    sumstats_columns: Mapping[str, str] | None = None,
    only_regions: set[str] | None = None,
    max_regions: int | None = None,
    verbose: bool = False,
) -> dict[str, object]:
    """Execute the GIM method; mGWAS is supplied input, never re-run."""

    progress("Validating parameters, genotype prefixes, and PLINK2")
    parameters.validate()
    if ancestry not in {"AFR", "AMR", "EAS", "EUR", "SAS", "CUSTOM"}:
        raise ValueError("ancestry must be one of AFR, AMR, EAS, EUR, SAS, or CUSTOM.")
    output, bfile, plink2 = Path(output), require_bfile(bfile), require_executable(plink2)
    ld_bfile = require_bfile(ld_bfile)
    panel_label = ld_panel_name or ld_bfile.name
    if output.exists() and any(path.name != "run.log" for path in output.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output}. Choose a new output directory.")
    output.mkdir(parents=True, exist_ok=True)
    progress("Defining candidate regions")
    region_result = define_regions(
        sumstats,
        sumstats_manifest=sumstats_manifest,
        sentinels_path=sentinels,
        ld_bfile=ld_bfile,
        plink2=plink2,
        parameters=parameters,
        columns=sumstats_columns,
        verbose=verbose,
    )
    regions = region_result["regions"]
    if only_regions is not None:
        regions = [row for row in regions if row["region_id"] in only_regions]
    if max_regions is not None:
        if max_regions < 1:
            raise ValueError("max_regions must be at least one.")
        regions = regions[:max_regions]
    selected_ids = {str(row["region_id"]) for row in regions}
    region_result["regions"] = regions
    region_result["region_metabolites"] = [row for row in region_result["region_metabolites"] if str(row["region_id"]) in selected_ids]
    progress(f"Starting individual-level conditional analysis for {len(regions):,} regions")
    conditional = run_conditional_analysis(
        regions=regions, region_metabolites=region_result["region_metabolites"], bfile=bfile, plink2=plink2,
        phenotypes=phenotypes, covariates=covariates, parameters=parameters, only_regions=selected_ids, verbose=verbose,
    )
    progress("Assembling significant SNP–trait edges into GIM connected components")
    components = components_from_matrix(conditional["matrix_out"], conditional_p=parameters.conditional_p) if conditional["matrix_out"] else {"edges": [], "members": [], "gim_summary": []}
    result: dict[str, object] = {**region_result, **conditional, **components}
    progress("Writing result tables")
    _write_result(output, result)  # type: ignore[arg-type]
    manifest = {
        "method": "GIMForge — Genetically Influenced Metabotypes",
        "reference": "Surendran et al., Nature Medicine 2022, DOI: 10.1038/s41591-022-02046-0",
        "mGWAS": (
            "Precomputed trait-specific clump leaders were used for region definition; no GWAS or sentinel clumping was run."
            if sentinels is not None
            else "Supplied summary statistics were used only for trait-specific sentinel clumping and region definition; no GWAS was run."
        ),
        "conditional_model": "PLINK2 additive linear regression with identical supplied covariates (variance-standardized for numerical stability); forward selection, full-model pruning, and ordered V_R x M_R matrix.",
        "genotype_backend": "PLINK 1 binary hard-call input. This backend cannot apply imputation INFO filtering; MAC filtering is applied.",
        "parameters": parameters.to_dict(),
        "inputs": {
            "sumstats": str(Path(sumstats).resolve()) if sumstats is not None else None,
            "sumstats_manifest": str(Path(sumstats_manifest).resolve()) if sumstats_manifest is not None else None,
            "sentinels": str(Path(sentinels).resolve()) if sentinels is not None else None,
            "analysis_bfile": str(bfile.resolve()),
            "ld_reference_bfile": str(ld_bfile.resolve()),
            "ld_reference_panel_name": panel_label,
            "ld_reference_ancestry": ancestry,
            "phenotypes": str(Path(phenotypes).resolve()),
            "covariates": str(Path(covariates).resolve()),
        },
        "ld_reference": "explicitly_supplied_reference_panel",
        "region_seed_input": (
            "precomputed_trait_specific_sentinels"
            if sentinels is not None
            else (
                "per_trait_summary_statistics_with_internal_trait_specific_clumping"
                if sumstats_manifest is not None
                else "stacked_summary_statistics_with_internal_trait_specific_clumping"
            )
        ),
        "tools": {"plink2": str(plink2.resolve()), "plink2_version": version(plink2)},
    }
    write_json(output / "run_manifest.json", manifest)
    progress("Writing interactive HTML report")
    write_report(
        output / "report.html", matrix_out=result["matrix_out"], members=result["members"],
        gim_summary=result["gim_summary"], edges=result["edges"], regions=result["regions"], conditional_p=parameters.conditional_p,
        metadata={
            "LD reference panel": panel_label,
            "LD reference ancestry": ancestry,
            "LD reference prefix": str(ld_bfile.resolve()),
            "Individual-level genotype": str(bfile.resolve()),
            "PLINK2": version(plink2),
            "Sentinel P threshold": parameters.sentinel_p,
            "Sentinel clumping": (
                "precomputed trait-specific leaders"
                if sentinels is not None
                else f"r² {parameters.sentinel_clump_r2} / {parameters.sentinel_clump_window_kb:,} kb"
            ),
            "LD span r² / window": f"{parameters.ld_span_r2} / {parameters.ld_window_kb:,} kb",
            "Cross-metabolite merge r²": parameters.cross_metabolite_merge_r2,
            "No-LD half-window / padding": f"{parameters.no_ld_half_width_kb:,} / {parameters.region_padding_kb:,} kb",
            "Conditional P threshold": parameters.conditional_p,
        },
    )
    progress(
        f"Run complete: {len(result['gim_summary']):,} GIMs; "
        f"report written to {(output / 'report.html').resolve()}"
    )
    return result
