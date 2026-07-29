import tempfile
import unittest
from pathlib import Path

from gimforge.pipeline import portal_variant_rows


class PortalVariantExportTests(unittest.TestCase):
    def test_exports_only_requested_bim_rows_without_sample_data(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "cohort"
            prefix.with_suffix(".bim").write_text(
                "1 rs1 0 101 A G\n"
                "2 rs2 0 202 C T\n"
                "3 unused 0 303 G A\n",
                encoding="utf-8",
            )

            rows = portal_variant_rows(prefix, {"rs2", "rs1", "missing"})

        self.assertEqual(
            rows,
            [
                {
                    "snp_id": "rs1",
                    "chromosome": "1",
                    "position": "101",
                    "allele1": "A",
                    "allele2": "G",
                },
                {
                    "snp_id": "rs2",
                    "chromosome": "2",
                    "position": "202",
                    "allele1": "C",
                    "allele2": "T",
                },
            ],
        )


if __name__ == "__main__":
    unittest.main()
