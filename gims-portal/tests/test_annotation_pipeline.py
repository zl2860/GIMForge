from __future__ import annotations

import sys
import json
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from alphagenome_scopes import target_scopes  # noqa: E402
from gims_portal import score_cache_compatible  # noqa: E402
from portal_pipeline import select_study_alternates  # noqa: E402
from split_alphagenome_portal_data import summarise_scope  # noqa: E402


class StudyAlleleTests(unittest.TestCase):
    def test_only_observed_alt_is_retained(self):
        alternates, orientation = select_study_alternates(
            "T", ["A", "G"], ["T", "G"]
        )
        self.assertEqual(alternates, ["G"])
        self.assertEqual(orientation, "forward")

    def test_reverse_complement_is_resolved(self):
        alternates, orientation = select_study_alternates(
            "T", ["C"], ["A", "G"]
        )
        self.assertEqual(alternates, ["C"])
        self.assertEqual(orientation, "reverse_complement")

    def test_unobserved_ensembl_alts_are_not_used_as_fallback(self):
        alternates, orientation = select_study_alternates(
            "T", ["A", "C"], ["T", "G"]
        )
        self.assertEqual(alternates, [])
        self.assertIsNone(orientation)


class BiologicalScopeTests(unittest.TestCase):
    def test_target_scope_classification(self):
        self.assertEqual(
            target_scopes(
                {"outputType": "ATAC", "biosample": "mucosa of stomach"}
            ),
            ["gastric_tissue"],
        )
        self.assertEqual(
            target_scopes({"outputType": "ATAC", "biosample": "AGS"}),
            ["gastric_cancer"],
        )
        self.assertEqual(
            target_scopes(
                {"outputType": "RNA_SEQ", "biosample": "CD14-positive monocyte"}
            ),
            ["immune"],
        )
        self.assertEqual(
            target_scopes({"outputType": "ATAC", "biosample": "neuron"}),
            [],
        )
        self.assertEqual(
            target_scopes({"outputType": "SPLICE_SITES", "track": "acceptor"}),
            ["shared"],
        )


class PrioritisationTests(unittest.TestCase):
    def test_absolute_quantile_preserves_strong_negative_effect(self):
        summary = summarise_scope(
            [
                {
                    "outputType": "ATAC",
                    "rankingScore": -0.99,
                    "rawScore": -0.1,
                    "track": "negative",
                },
                {
                    "outputType": "ATAC",
                    "rankingScore": 0.80,
                    "rawScore": 0.2,
                    "track": "positive",
                },
            ]
        )
        self.assertEqual(summary["maxAbsRankingScore"], 0.99)
        self.assertEqual(summary["rankingScoreAtMaximum"], -0.99)
        self.assertEqual(summary["topTrack"], "negative")

    def test_official_merged_splicing_formula(self):
        summary = summarise_scope(
            [
                {
                    "alternate": "G",
                    "geneId": "ENSG1",
                    "outputType": "SPLICE_SITES",
                    "rawScore": 0.2,
                },
                {
                    "alternate": "G",
                    "geneId": "ENSG1",
                    "outputType": "SPLICE_SITE_USAGE",
                    "rawScore": 0.3,
                },
                {
                    "alternate": "G",
                    "geneId": "ENSG1",
                    "outputType": "SPLICE_JUNCTIONS",
                    "rawScore": 0.5,
                },
            ]
        )
        self.assertAlmostEqual(summary["splicingCombined"]["score"], 0.6)


class ReproducibilityTests(unittest.TestCase):
    def test_legacy_100kb_cache_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            annotation_dir = Path(temporary)
            (annotation_dir / "alphagenome_coverage.json").write_text(
                json.dumps(
                    {
                        "schemaVersion": 1,
                        "sequenceLength": "100KB",
                        "completed": True,
                    }
                ),
                encoding="utf-8",
            )
            self.assertFalse(
                score_cache_compatible(annotation_dir, "new-input-digest")
            )


if __name__ == "__main__":
    unittest.main()
