from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gims_portal import (  # noqa: E402
    check_alpha_dependencies,
    prepare_bundle,
    request_api_key,
    score_cache_ids,
)
from portal_pipeline import build_data, write_alphagenome_input  # noqa: E402


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class PortalWorkflowTest(unittest.TestCase):
    def test_completed_score_ids_do_not_require_loading_full_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            annotation_dir = Path(temporary)
            scores_dir = annotation_dir / "alphagenome_scores"
            scores_dir.mkdir()
            (scores_dir / "rs1.json").write_text("[]", encoding="utf-8")
            (annotation_dir / "alphagenome_coverage.json").write_text(
                json.dumps(
                    {
                        "schemaVersion": 2,
                        "sequenceLength": "1MB",
                        "modelVersion": "ALL_FOLDS",
                        "sdkVersion": "0.7.0",
                        "inputSha256": "test-sha",
                        "scopeFilterVersion": 2,
                        "retainsOnlyTargetScopes": True,
                        "scoredSnpIds": ["rs1"],
                    }
                ),
                encoding="utf-8",
            )
            with patch(
                "gims_portal.existing_score_ids",
                side_effect=AssertionError("full summary should not be loaded"),
            ):
                ids = score_cache_ids(
                    annotation_dir,
                    {"rs1"},
                    "test-sha",
                )
            self.assertEqual(ids, {"rs1"})

    def test_annotation_dependencies_are_checked_before_api_prompt(self):
        failed = SimpleNamespace(
            returncode=1,
            stdout="",
            stderr="ModuleNotFoundError: No module named 'alphagenome'",
        )
        with patch("gims_portal.subprocess.run", return_value=failed):
            with self.assertRaisesRegex(
                RuntimeError,
                r"Install them before entering an API key",
            ) as raised:
                check_alpha_dependencies(Path("/tmp/alpha-python"))
        self.assertIn("-m pip install", str(raised.exception))
        self.assertIn("pandas numpy", str(raised.exception))

    def test_annotation_dependency_versions_are_reported(self):
        passed = SimpleNamespace(
            returncode=0,
            stdout=json.dumps(
                {
                    "alphagenome": "0.7.0",
                    "pandas": "3.0.0",
                    "numpy": "2.0.0",
                }
            ),
            stderr="",
        )
        with patch("gims_portal.subprocess.run", return_value=passed):
            versions = check_alpha_dependencies(Path("/tmp/alpha-python"))
        self.assertEqual(versions["alphagenome"], "0.7.0")

    def test_api_key_can_be_supplied_without_command_arguments(self):
        original = os.environ.get("TEST_ALPHAGENOME_KEY")
        os.environ["TEST_ALPHAGENOME_KEY"] = "memory-only-secret"
        try:
            secret, prompted = request_api_key("TEST_ALPHAGENOME_KEY")
        finally:
            if original is None:
                os.environ.pop("TEST_ALPHAGENOME_KEY", None)
            else:
                os.environ["TEST_ALPHAGENOME_KEY"] = original
        self.assertEqual(secret, "memory-only-secret")
        self.assertFalse(prompted)

    def test_standard_gimforge_bundle_builds_without_original_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / "gimforge"
            bundle = root / "bundle"
            output = root / "web-data"
            results.mkdir()

            matrix = {
                "region_id": "GIMForge_region_001",
                "marker_order": "1",
                "snp_id": "rs1",
                "metabolite": "metabolite_a",
                "a1_freq": "0.2",
                "maf": "0.2",
                "maf_class": "common",
                "beta": "0.5",
                "se": "0.1",
                "p": "1e-8",
                "n": "100",
                "conditioned_on_n": "0",
                "testable": "TRUE",
            }
            write_tsv(results / "matrix_out.tsv", list(matrix), [matrix])
            write_tsv(
                results / "edges.tsv",
                [*matrix, "gim_id"],
                [{**matrix, "gim_id": "GIMForge_region_001_GIM_001"}],
            )
            write_tsv(
                results / "members.tsv",
                [
                    "gim_id",
                    "region_id",
                    "node_type",
                    "node_id",
                    "marker_order",
                    "maf",
                    "maf_class",
                ],
                [
                    {
                        "gim_id": "GIMForge_region_001_GIM_001",
                        "region_id": "GIMForge_region_001",
                        "node_type": "SNP",
                        "node_id": "rs1",
                        "marker_order": "1",
                        "maf": "0.2",
                        "maf_class": "common",
                    },
                    {
                        "gim_id": "GIMForge_region_001_GIM_001",
                        "region_id": "GIMForge_region_001",
                        "node_type": "metabolite",
                        "node_id": "metabolite_a",
                        "marker_order": "",
                        "maf": "",
                        "maf_class": "",
                    },
                ],
            )
            write_tsv(
                results / "gim_summary.tsv",
                [
                    "gim_id",
                    "region_id",
                    "n_snps",
                    "n_metabolites",
                    "snps",
                    "metabolites",
                ],
                [
                    {
                        "gim_id": "GIMForge_region_001_GIM_001",
                        "region_id": "GIMForge_region_001",
                        "n_snps": "1",
                        "n_metabolites": "1",
                        "snps": "rs1",
                        "metabolites": "metabolite_a",
                    }
                ],
            )
            write_tsv(
                results / "independent_signals.tsv",
                [
                    "region_id",
                    "metabolite",
                    "snp_id",
                    "selection_source",
                    "maf",
                    "maf_class",
                    "beta",
                    "se",
                    "p",
                    "n",
                ],
                [
                    {
                        "region_id": "GIMForge_region_001",
                        "metabolite": "metabolite_a",
                        "snp_id": "rs1",
                        "selection_source": "forward_then_full_model",
                        "maf": "0.2",
                        "maf_class": "common",
                        "beta": "0.5",
                        "se": "0.1",
                        "p": "1e-8",
                        "n": "100",
                    }
                ],
            )
            write_tsv(
                results / "regions.tsv",
                [
                    "region_id",
                    "chromosome",
                    "start",
                    "end",
                    "start_before_margin",
                    "end_before_margin",
                    "n_initial_regions",
                ],
                [
                    {
                        "region_id": "GIMForge_region_001",
                        "chromosome": "1",
                        "start": "50",
                        "end": "150",
                        "start_before_margin": "75",
                        "end_before_margin": "125",
                        "n_initial_regions": "1",
                    }
                ],
            )
            write_tsv(
                results / "variants.tsv",
                [
                    "snp_id",
                    "chromosome",
                    "position",
                    "allele1",
                    "allele2",
                ],
                [
                    {
                        "snp_id": "rs1",
                        "chromosome": "1",
                        "position": "100",
                        "allele1": "A",
                        "allele2": "G",
                    }
                ],
            )
            (results / "run_manifest.json").write_text(
                json.dumps(
                    {
                        "method": "GIMForge",
                        "inputs": {},
                    }
                ),
                encoding="utf-8",
            )

            args = SimpleNamespace(
                results=results,
                bundle=bundle,
                bim=None,
                source_build="GRCh37",
                matrix_out=None,
                edges=None,
                members=None,
                gim_summary=None,
                independent_signals=None,
                regions=None,
                region_summary=None,
                region_metabolites=None,
                run_manifest=None,
                differential_metabolomics=None,
                single_cell_module_matrix=None,
                single_cell_module_scores=None,
                reuse_annotations=None,
            )
            prepare_bundle(args, {}, root)
            manifest = json.loads(
                (bundle / "bundle_manifest.json").read_text(encoding="utf-8")
            )
            self.assertFalse(manifest["containsIndividualLevelData"])
            self.assertEqual(manifest["nGimSnps"], 1)
            self.assertEqual(manifest["nSourceCoordinates"], 1)
            self.assertEqual(
                manifest["variantSource"], "gimforge_variants_table"
            )

            annotation_dir = bundle / "annotations"
            annotation_dir.mkdir()
            (annotation_dir / "ensembl_annotation_cache.json").write_text(
                json.dumps(
                    {
                        "rs1": {
                            "variation": {
                                "var_class": "SNP",
                                "MAF": 0.18,
                                "mappings": [
                                    {
                                        "assembly_name": "GRCh38",
                                        "seq_region_name": "1",
                                        "start": 110,
                                        "allele_string": "A/G",
                                    }
                                ],
                            },
                            "vep": {
                                "most_severe_consequence": "intron_variant",
                                "transcript_consequences": [
                                    {
                                        "gene_id": "ENSG000001",
                                        "gene_symbol": "GENE1",
                                        "transcript_id": "ENST000001",
                                        "consequence_terms": ["intron_variant"],
                                    }
                                ],
                            },
                        }
                    }
                ),
                encoding="utf-8",
            )
            alpha_input = write_alphagenome_input(bundle)
            self.assertIn("chr1", alpha_input.read_text(encoding="utf-8"))

            stats = build_data(bundle, output)
            self.assertEqual(stats["nGimEntities"], 1)
            self.assertEqual(stats["nIndependentSnps"], 1)
            self.assertEqual(stats["annotationStatus"], "resolved")
            entity = json.loads(
                (output / "gim_entities.json").read_text(encoding="utf-8")
            )[0]
            self.assertEqual(entity["heatmap"][0]["p"], 1e-8)
            snp = json.loads(
                (output / "gim_snps.json").read_text(encoding="utf-8")
            )[0]
            self.assertEqual(snp["studyMaf"], 0.2)
            self.assertEqual(snp["studyMafClass"], "common")


if __name__ == "__main__":
    unittest.main()
