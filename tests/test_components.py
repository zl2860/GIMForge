import unittest

from gimforge.components import components_from_matrix


class ComponentDefinitionTests(unittest.TestCase):
    def test_significant_bipartite_components_are_gims(self):
        matrix = [
            {"region_id": "R", "marker_order": 1, "snp_id": "v1", "metabolite": "a", "p": 1e-10},
            {"region_id": "R", "marker_order": 2, "snp_id": "v2", "metabolite": "a", "p": 1e-9},
            {"region_id": "R", "marker_order": 3, "snp_id": "v3", "metabolite": "b", "p": 1e-12},
            {"region_id": "R", "marker_order": 1, "snp_id": "v1", "metabolite": "b", "p": 0.4},
        ]
        result = components_from_matrix(matrix, conditional_p=1e-8)
        self.assertEqual(len(result["gim_summary"]), 2)
        first, second = result["gim_summary"]
        self.assertEqual((first["snps"], first["metabolites"]), ("v1;v2", "a"))
        self.assertEqual((second["snps"], second["metabolites"]), ("v3", "b"))
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
            }
        ]
        result = components_from_matrix(matrix, conditional_p=1e-8)
        row = result["matrix_out"][0]
        self.assertEqual(row["p"], 1e-10)
        self.assertEqual(row["snp_id"], "v1")
        self.assertEqual(row["metabolite"], "a")
        self.assertEqual(row["region_id"], "R")
        self.assertEqual(row["beta"], "0.25")
        self.assertEqual(len(result["edges"]), 1)


if __name__ == "__main__":
    unittest.main()
