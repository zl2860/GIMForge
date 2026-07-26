import tempfile
import unittest
from pathlib import Path

from gimforge.regions import normalise_sumstats_manifest


class SummaryStatisticInputTests(unittest.TestCase):
    def test_manifest_reads_one_file_per_trait_without_concatenation(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_sumstats_test_") as directory:
            root = Path(directory)
            (root / "M001.tsv").write_text(
                "#CHROM\tPOS\tID\tP\n1\t100\trs1\t1e-12\n1\t200\trs2\t0.5\n"
            )
            (root / "M002.tsv").write_text(
                "#CHROM\tPOS\tID\tP\n2\t300\trs3\t2e-13\n"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "trait_id\tpath\nM001\tM001.tsv\nM002\tM002.tsv\n"
            )
            rows = normalise_sumstats_manifest(manifest, region_p=1e-10)
        self.assertEqual(
            [(row["metabolite"], row["snp_id"]) for row in rows],
            [("M001", "rs1"), ("M002", "rs3")],
        )


if __name__ == "__main__":
    unittest.main()
