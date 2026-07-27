import unittest

from gimforge.components import components_from_matrix


class ComponentDefinitionTests(unittest.TestCase):
    def test_significant_bipartite_components_are_gims(self):
        matrix = [
            {"region_id": "R", "marker_order": 1, "snp_id": "v1", "metabolite": "a", "maf": 0.005, "p": 1e-10},
            {"region_id": "R", "marker_order": 2, "snp_id": "v2", "metabolite": "a", "maf": 0.03, "p": 1e-9},
            {"region_id": "R", "marker_order": 3, "snp_id": "v3", "metabolite": "b", "maf": 0.2, "p": 1e-12},
            {"region_id": "R", "marker_order": 1, "snp_id": "v1", "metabolite": "b", "maf": 0.005, "p": 0.4},
        ]
        result = components_from_matrix(matrix, conditional_p=1e-8)
        self.assertEqual(len(result["gim_summary"]), 2)
        first, second = result["gim_summary"]
        self.assertEqual((first["snps"], first["metabolites"]), ("v1;v2", "a"))
        self.assertEqual((second["snps"], second["metabolites"]), ("v3", "b"))
        self.assertEqual(first["n_rare_snps"], 1)
        self.assertEqual(first["n_low_frequency_snps"], 1)
        self.assertEqual(first["n_common_snps"], 0)
        self.assertEqual(second["n_common_snps"], 1)
        self.assertEqual(first["snp_maf_classes"], "v1=rare;v2=low_frequency")
        snp_members = [
            row for row in result["members"] if row["node_type"] == "SNP"
        ]
        self.assertEqual(
            {row["node_id"]: row["maf_class"] for row in snp_members},
            {"v1": "rare", "v2": "low_frequency", "v3": "common"},
        )
        self.assertEqual(len(result["edges"]), 3)

    def test_import_normalises_uppercase_p_for_report_and_output(self):
        matrix = [
            {
                "region": "R",
                "order": "1",
                "ID": "v1",
                "trait": "a",
                "BETA": "0.25",
                "SE": "0.05",
                "P": "1e-10",
                "N": "100",
                "A1_FREQ": "0.9",
            }
        ]
        result = components_from_matrix(matrix, conditional_p=1e-8)
        row = result["matrix_out"][0]
        self.assertEqual(row["p"], 1e-10)
        self.assertEqual(row["snp_id"], "v1")
        self.assertEqual(row["metabolite"], "a")
        self.assertEqual(row["region_id"], "R")
        self.assertEqual(row["beta"], "0.25")
        self.assertAlmostEqual(row["maf"], 0.1)
        self.assertEqual(row["maf_class"], "common")
        self.assertEqual(len(result["edges"]), 1)


if __name__ == "__main__":
    unittest.main()
