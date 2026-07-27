import unittest

from gimforge.frequency import classify_maf, maf_from_a1_frequency


class FrequencyTests(unittest.TestCase):
    def test_reference_maf_class_boundaries(self):
        self.assertEqual(classify_maf(0.01), "rare")
        self.assertEqual(classify_maf(0.0100001), "low_frequency")
        self.assertEqual(classify_maf(0.05), "low_frequency")
        self.assertEqual(classify_maf(0.0500001), "common")
        self.assertEqual(classify_maf(""), "unknown")

    def test_maf_is_derived_from_either_allele_frequency(self):
        self.assertAlmostEqual(maf_from_a1_frequency(0.008), 0.008)
        self.assertAlmostEqual(maf_from_a1_frequency(0.992), 0.008)
        self.assertIsNone(maf_from_a1_frequency(1.2))


if __name__ == "__main__":
    unittest.main()
