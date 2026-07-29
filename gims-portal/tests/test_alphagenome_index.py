import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "split_alphagenome_portal_data.py"
SPEC = importlib.util.spec_from_file_location("split_alphagenome_portal_data", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AlphaGenomeIndexTest(unittest.TestCase):
    def test_scope_summary_preserves_multimodal_ranking_evidence(self):
        rows = [
            {
                "outputType": "ATAC",
                "rankingScore": 0.91,
                "rawScore": -2.5,
                "track": "UBERON:0000945 ATAC-seq",
                "biosample": "stomach",
            },
            {
                "outputType": "ATAC",
                "rankingScore": 0.82,
                "rawScore": 4.0,
                "track": "other ATAC-seq",
                "biosample": "other",
            },
            {
                "outputType": "RNA_SEQ",
                "rankingScore": 0.97,
                "rawScore": 0.4,
                "track": "UBERON:0000945 RNA-seq",
                "biosample": "stomach",
            },
        ]

        summary = MODULE.summarise_scope(rows)

        self.assertEqual(summary["nTracks"], 3)
        self.assertEqual(summary["nModalities"], 2)
        self.assertEqual(summary["modalities"], ["ATAC", "RNA_SEQ"])
        self.assertEqual(summary["maxRankingScore"], 0.97)
        self.assertEqual(summary["modalityStats"][0]["outputType"], "RNA_SEQ")
        atac = next(
            item
            for item in summary["modalityStats"]
            if item["outputType"] == "ATAC"
        )
        self.assertEqual(atac["nTracks"], 2)
        self.assertEqual(atac["maxAbsRawScore"], 4.0)


if __name__ == "__main__":
    unittest.main()
