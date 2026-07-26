"""Command-line interface for Linux execution."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .components import components_from_matrix
from .dependencies import doctor, resolve_bolt, resolve_plink2
from .io import read_json, read_table, write_json, write_table
from .parameters import parameters_from_args
from .pipeline import _FIELDS, run_gim
from .plink import version
from .progress import configure_progress_log, progress, start_progress_session
from .regions import clump_sentinels
from .report import write_report

_ANCESTRIES = ("AFR", "AMR", "EAS", "EUR", "SAS", "CUSTOM")


def _sumstats_column_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sumstats-metabolite-column", default="metabolite")
    parser.add_argument("--sumstats-chromosome-column", default="chromosome")
    parser.add_argument("--sumstats-position-column", default="position")
    parser.add_argument("--sumstats-snp-column", default="snp_id")
    parser.add_argument("--sumstats-p-column", default="p")


def _sumstats_columns(args: argparse.Namespace) -> dict[str, str]:
    return {
        "metabolite": args.sumstats_metabolite_column,
        "chromosome": args.sumstats_chromosome_column,
        "position": args.sumstats_position_column,
        "snp_id": args.sumstats_snp_column,
        "p": args.sumstats_p_column,
    }


def _parameter_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("GIM parameters")
    group.add_argument("--sentinel-p", type=float, default=1.25e-11, help="mGWAS P threshold for metabolite-specific sentinel candidates")
    group.add_argument("--conditional-p", type=float, default=1.24741348813236e-8, help="threshold for conditional selection and GIM edges")
    group.add_argument("--sentinel-clump-r2", type=float, default=0.1, help="r² for same-metabolite sentinel clumping")
    group.add_argument("--sentinel-clump-window-kb", type=int, default=1_000_000, help="physical window for sentinel clumping")
    group.add_argument("--ld-span-r2", type=float, default=0.1, help="r² threshold for sentinel LD spans")
    group.add_argument("--cross-metabolite-merge-r2", type=float, default=0.6, help="r² for cross-metabolite sentinel merging")
    group.add_argument("--ld-window-kb", type=int, default=1_000_000, help="LD search window used to construct sentinel spans")
    group.add_argument("--no-ld-half-width-kb", type=int, default=500, help="sentinel half-width when no LD neighbour exists")
    group.add_argument("--region-padding-kb", type=int, default=250, help="padding added after merging LD spans")
    group.add_argument("--mac-min", type=int, default=10)
    group.add_argument("--geno-missing-max", type=float, default=None)
    group.add_argument("--threads", type=int, default=1)
    group.add_argument("--metabolite-batch-size", type=int, default=1, help="PLINK phenotypes per temporary result batch; larger values trade disk for speed")
    group.add_argument(
        "--genetic-model",
        choices=("additive", "dominant", "recessive"),
        default="additive",
        help="genotype coding used for both tested and conditioning variants",
    )
    group.add_argument(
        "--regression-model",
        choices=("linear", "mixed"),
        default="linear",
        help="linear uses PLINK2; mixed uses additive BOLT-LMM-inf",
    )
    group.add_argument(
        "--mixed-backend",
        choices=("bolt-lmm",),
        default="bolt-lmm",
        help="mixed-model engine",
    )
    group.add_argument("--no-force-single-forward-lead", dest="force_single_forward_lead", action="store_false")
    group.set_defaults(force_single_forward_lead=True)


def _parameters(args: argparse.Namespace):
    return parameters_from_args({
        "sentinel_p": args.sentinel_p, "conditional_p": args.conditional_p,
        "sentinel_clump_r2": args.sentinel_clump_r2,
        "sentinel_clump_window_kb": args.sentinel_clump_window_kb,
        "ld_span_r2": args.ld_span_r2,
        "cross_metabolite_merge_r2": args.cross_metabolite_merge_r2,
        "ld_window_kb": args.ld_window_kb,
        "no_ld_half_width_kb": args.no_ld_half_width_kb,
        "region_padding_kb": args.region_padding_kb,
        "mac_min": args.mac_min, "geno_missing_max": args.geno_missing_max, "threads": args.threads,
        "metabolite_batch_size": args.metabolite_batch_size,
        "force_single_forward_lead": args.force_single_forward_lead,
        "genetic_model": args.genetic_model,
        "regression_model": args.regression_model,
        "mixed_backend": args.mixed_backend,
    })


def _write_component_result(
    output: Path,
    matrix: list[dict[str, str]],
    conditional_p: float,
    reference_panel: str | None = None,
    ancestry: str | None = None,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    components = components_from_matrix(matrix, conditional_p=conditional_p)
    for name in ("edges", "members", "gim_summary"):
        suffix = ".tsv" if name == "gim_summary" else ".tsv.gz"
        write_table(output / f"{name}{suffix}", components[name], fieldnames=_FIELDS[name])
    write_table(output / "matrix_out.tsv.gz", components["matrix_out"], fieldnames=_FIELDS["matrix_out"])
    write_json(
        output / "run_manifest.json",
        {
            "method": "GIMForge components from an existing ordered conditional matrix",
            "conditional_p": conditional_p,
            "inputs": {
                "ld_reference_bfile": "not available from matrix-only input",
                "ld_reference_panel_name": reference_panel or "not recorded in matrix input",
                "ld_reference_ancestry": ancestry or "not recorded in matrix input",
            },
        },
    )
    write_report(
        output / "report.html", matrix_out=components["matrix_out"], members=components["members"],
        gim_summary=components["gim_summary"], edges=components["edges"], conditional_p=conditional_p,
        metadata={
            "Input": "existing ordered matrix_out",
            "LD reference panel": reference_panel or "not recorded in matrix input",
            "LD reference ancestry": ancestry or "not recorded in matrix input",
            "Conditional P threshold": conditional_p,
        },
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="gimforge", description="Define genetically influenced metabotypes from supplied mGWAS summaries and individual-level data.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    clump_parser = subparsers.add_parser("clump", help="select trait-specific sentinel variants with an explicit LD reference panel")
    clump_sumstats = clump_parser.add_mutually_exclusive_group(required=True)
    clump_sumstats.add_argument("--sumstats", help="stacked long-format multi-trait summary-statistic table")
    clump_sumstats.add_argument("--sumstats-manifest", help="TSV with trait_id and path, one summary-statistic file per trait")
    clump_parser.add_argument("--ld-bfile", required=True, help="PLINK 1 binary population-matched LD reference-panel prefix")
    clump_parser.add_argument("--ancestry", required=True, choices=_ANCESTRIES, help="ancestry represented by --ld-bfile")
    clump_parser.add_argument("--ld-panel-name", help="human-readable panel release/name recorded in the manifest")
    clump_parser.add_argument("--out", required=True, help="new output directory for sentinels.tsv, run_manifest.json, and run.log")
    clump_parser.add_argument("--plink2", help="PLINK2 executable; defaults to plink2 in PATH")
    clump_parser.add_argument("--sentinel-p", type=float, default=1.25e-11)
    clump_parser.add_argument("--sentinel-clump-r2", type=float, default=0.1)
    clump_parser.add_argument("--sentinel-clump-window-kb", type=int, default=1_000_000)
    clump_parser.add_argument("--mac-min", type=int, default=10)
    clump_parser.add_argument("--geno-missing-max", type=float, default=None)
    clump_parser.add_argument("--threads", type=int, default=1)
    clump_parser.add_argument("--verbose", action="store_true")
    _sumstats_column_arguments(clump_parser)
    run_parser = subparsers.add_parser("run", help="define regions, run conditional analysis, and write GIMs plus report.html")
    run_sumstats = run_parser.add_mutually_exclusive_group(required=True)
    run_sumstats.add_argument("--sumstats", help="stacked long-format multi-trait summary-statistic table")
    run_sumstats.add_argument("--sumstats-manifest", help="TSV with trait_id and path, one summary-statistic file per trait")
    run_sumstats.add_argument("--sentinels", help="precomputed trait-specific clump leaders from gimforge clump")
    run_parser.add_argument("--bfile", required=True, help="PLINK 1 binary prefix")
    run_parser.add_argument("--ld-bfile", required=True, help="PLINK 1 binary population-matched LD reference-panel prefix")
    run_parser.add_argument("--ancestry", required=True, choices=_ANCESTRIES, help="ancestry represented by --ld-bfile")
    run_parser.add_argument("--ld-panel-name", help="human-readable panel label recorded in the manifest/report, e.g. '1000 Genomes Phase 3 EUR'")
    run_parser.add_argument("--phenotypes", required=True, help="uncompressed headered TSV: FID IID plus numeric trait columns")
    run_parser.add_argument("--covariates", required=True, help="uncompressed headered TSV: FID IID plus numeric covariates")
    run_parser.add_argument("--out", required=True)
    run_parser.add_argument("--plink2", help="PLINK2 executable; defaults to plink2 in PATH")
    run_parser.add_argument("--bolt", help="BOLT-LMM executable; required with --regression-model mixed")
    run_parser.add_argument(
        "--bolt-model-bfile",
        help="genome-wide PLINK BED/BIM/FAM prefix for the BOLT random effect; defaults to --bfile",
    )
    run_parser.add_argument(
        "--bolt-genetic-map",
        help="BOLT-LMM genetic map matching the analysis genome build; required for mixed",
    )
    _sumstats_column_arguments(run_parser)
    run_parser.add_argument("--only-region", action="append", help="run one or more already-numbered regions")
    run_parser.add_argument("--max-regions", type=int, help="run the first N regions; useful for a smoke test")
    run_parser.add_argument("--verbose", action="store_true")
    _parameter_arguments(run_parser)
    component_parser = subparsers.add_parser("components", help="create GIMs and report.html from an existing ordered matrix_out")
    component_parser.add_argument("--matrix-out", required=True)
    component_parser.add_argument("--out", required=True)
    component_parser.add_argument("--conditional-p", type=float, default=1.24741348813236e-8)
    component_parser.add_argument("--reference-panel", help="reference-panel provenance label for a precomputed matrix")
    component_parser.add_argument("--ancestry", choices=_ANCESTRIES, help="reference-panel ancestry for a precomputed matrix")
    report_parser = subparsers.add_parser("report", help="(re)write report.html from a completed result directory")
    report_parser.add_argument("--results", required=True)
    report_parser.add_argument("--out", help="HTML path; defaults to RESULTS/report.html")
    report_parser.add_argument("--conditional-p", type=float, help="override the threshold recorded in run_manifest.json")
    doctor_parser = subparsers.add_parser(
        "doctor", help="check Python, PLINK2, and optional mixed-model dependencies"
    )
    doctor_parser.add_argument("--plink2", help="PLINK2 executable; defaults to plink2 in PATH")
    doctor_parser.add_argument("--regression-model", choices=("linear", "mixed"), default="linear")
    doctor_parser.add_argument("--bolt", help="BOLT-LMM executable; checked when regression model is mixed")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = ["gimforge", *(argv if argv is not None else sys.argv[1:])]
    try:
        if args.command == "clump":
            output = Path(args.out)
            if output.exists() and any(path.name != "run.log" for path in output.iterdir()):
                raise FileExistsError(f"Output directory is not empty: {output}. Choose a new output directory.")
            output.mkdir(parents=True, exist_ok=True)
            configure_progress_log(output / "run.log")
            start_progress_session(
                version=__version__, command=command, log_path=output / "run.log"
            )
            progress("Starting trait-specific sentinel clumping")
            progress("Checking PLINK2 dependency")
            executable = resolve_plink2(args.plink2)
            parameters = parameters_from_args(
                {
                    "sentinel_p": args.sentinel_p,
                    "sentinel_clump_r2": args.sentinel_clump_r2,
                    "sentinel_clump_window_kb": args.sentinel_clump_window_kb,
                    "mac_min": args.mac_min,
                    "geno_missing_max": args.geno_missing_max,
                    "threads": args.threads,
                }
            )
            sentinels = clump_sentinels(
                args.sumstats,
                sumstats_manifest=args.sumstats_manifest,
                ld_bfile=args.ld_bfile,
                plink2=executable,
                parameters=parameters,
                columns=_sumstats_columns(args),
                verbose=args.verbose,
            )
            write_table(output / "sentinels.tsv", sentinels, fieldnames=_FIELDS["sentinels"])
            write_json(
                output / "run_manifest.json",
                {
                    "method": "GIMForge trait-specific sentinel LD clumping",
                    "parameters": parameters.to_dict(),
                    "inputs": {
                        "sumstats": str(Path(args.sumstats).resolve()) if args.sumstats else None,
                        "sumstats_manifest": str(Path(args.sumstats_manifest).resolve()) if args.sumstats_manifest else None,
                        "ld_reference_bfile": str(Path(args.ld_bfile).resolve()),
                        "ld_reference_panel_name": args.ld_panel_name or Path(args.ld_bfile).name,
                        "ld_reference_ancestry": args.ancestry,
                    },
                    "tools": {
                        "plink2": str(executable.resolve()),
                        "plink2_version": version(executable),
                    },
                    "outputs": {"n_sentinels": len(sentinels)},
                },
            )
            progress(f"Clumping outputs written to {output.resolve()}")
        elif args.command == "run":
            output = Path(args.out)
            if output.exists() and any(path.name != "run.log" for path in output.iterdir()):
                raise FileExistsError(f"Output directory is not empty: {output}. Choose a new output directory.")
            output.mkdir(parents=True, exist_ok=True)
            configure_progress_log(output / "run.log")
            start_progress_session(
                version=__version__, command=command, log_path=output / "run.log"
            )
            progress("Starting GIMForge run")
            parameters = _parameters(args)
            progress(
                "Conditional configuration: "
                f"regression={parameters.regression_model}; "
                f"genetic_model={parameters.genetic_model}; "
                f"P<={parameters.conditional_p:.6g}"
            )
            progress("Checking PLINK2 dependency")
            if args.regression_model == "mixed":
                progress("Checking BOLT-LMM dependency")
            run_gim(
                sumstats=args.sumstats, sumstats_manifest=args.sumstats_manifest, sentinels=args.sentinels,
                bfile=args.bfile, ld_bfile=args.ld_bfile,
                ancestry=args.ancestry, ld_panel_name=args.ld_panel_name,
                phenotypes=args.phenotypes, covariates=args.covariates,
                output=args.out, plink2=resolve_plink2(args.plink2), parameters=parameters,
                bolt=resolve_bolt(args.bolt) if args.regression_model == "mixed" else None,
                bolt_model_bfile=args.bolt_model_bfile,
                bolt_genetic_map=args.bolt_genetic_map,
                sumstats_columns=_sumstats_columns(args),
                only_regions=set(args.only_region) if args.only_region else None, max_regions=args.max_regions, verbose=args.verbose,
            )
        elif args.command == "components":
            output = Path(args.out)
            if output.exists() and any(path.name != "run.log" for path in output.iterdir()):
                raise FileExistsError(f"Output directory is not empty: {output}. Choose a new output directory.")
            output.mkdir(parents=True, exist_ok=True)
            configure_progress_log(output / "run.log")
            start_progress_session(
                version=__version__, command=command, log_path=output / "run.log"
            )
            progress("Reading and normalising existing ordered conditional matrix")
            _write_component_result(output, read_table(args.matrix_out), args.conditional_p, args.reference_panel, args.ancestry)
            progress(f"Component outputs and report written to {output.resolve()}")
        elif args.command == "report":
            results = Path(args.results)
            configure_progress_log(results / "run.log")
            start_progress_session(
                version=__version__, command=command, log_path=results / "run.log"
            )
            progress("Regenerating interactive HTML report")
            matrix = read_table(results / "matrix_out.tsv.gz")
            edges_path = results / "edges.tsv.gz"
            edges = read_table(edges_path) if edges_path.is_file() else []
            members = read_table(results / "members.tsv.gz")
            summary = read_table(results / "gim_summary.tsv")
            regions_path = results / "regions.tsv"
            manifest_path = results / "run_manifest.json"
            manifest = read_json(manifest_path) if manifest_path.is_file() else {}
            inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}
            tools = manifest.get("tools", {}) if isinstance(manifest, dict) else {}
            recorded_parameters = manifest.get("parameters", {}) if isinstance(manifest, dict) else {}
            recorded_p = recorded_parameters.get("conditional_p") if isinstance(recorded_parameters, dict) else None
            if recorded_p is None and isinstance(manifest, dict):
                recorded_p = manifest.get("conditional_p")
            conditional_p = args.conditional_p if args.conditional_p is not None else float(recorded_p or 1.24741348813236e-8)
            metadata = {
                "LD reference panel": inputs.get("ld_reference_panel_name", "not recorded"),
                "LD reference ancestry": inputs.get("ld_reference_ancestry", "not recorded"),
                "LD reference prefix": inputs.get("ld_reference_bfile", "not recorded"),
                "Individual-level genotype": inputs.get("analysis_bfile", "not recorded"),
                "PLINK2": tools.get("plink2_version", "not recorded"),
                "BOLT-LMM": tools.get("bolt_version") or "not used",
                "Genetic model": recorded_parameters.get("genetic_model", "additive"),
                "Regression model": recorded_parameters.get("regression_model", "linear"),
                "Conditional P threshold": conditional_p,
            }
            write_report(
                args.out or results / "report.html", matrix_out=matrix, members=members,
                gim_summary=summary, edges=edges, regions=read_table(regions_path) if regions_path.is_file() else [],
                conditional_p=conditional_p, metadata=metadata,
            )
            progress(f"Report written to {Path(args.out or results / 'report.html').resolve()}")
        else:
            healthy, rows = doctor(
                args.plink2,
                regression_model=args.regression_model,
                bolt=args.bolt,
            )
            for dependency, state, detail in rows:
                print(f"{dependency}\t{state}\t{detail}")
            return 0 if healthy else 2
    except (FileNotFoundError, FileExistsError, ValueError, RuntimeError) as error:
        progress(f"ERROR: {error}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
