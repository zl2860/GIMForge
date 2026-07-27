import unittest

from gimforge.plink import common_filters


class PlinkFilterTests(unittest.TestCase):
    def test_optional_genotype_qc_filters_are_forwarded(self):
        arguments = common_filters(
            mac_min=11,
            maf_min=0.0001,
            hwe_p_min=1e-6,
            geno_missing_max=0.05,
        )

        self.assertEqual(
            arguments,
            [
                "--chr",
                "1-22",
                "--snps-only",
                "just-acgt",
                "--max-alleles",
                "2",
                "--mac",
                "11",
                "--maf",
                "0.0001",
                "--hwe",
                "1e-06",
                "--geno",
                "0.05",
            ],
        )

    def test_optional_genotype_qc_filters_can_be_disabled(self):
        arguments = common_filters(
            mac_min=10,
            maf_min=None,
            hwe_p_min=None,
            geno_missing_max=None,
            chromosome=None,
        )

        self.assertNotIn("--maf", arguments)
        self.assertNotIn("--hwe", arguments)
        self.assertNotIn("--geno", arguments)


if __name__ == "__main__":
    unittest.main()
