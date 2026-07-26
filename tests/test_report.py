import tempfile
import unittest
from pathlib import Path

from gimforge.report import write_report


class ReportTests(unittest.TestCase):
    def test_report_contains_interactive_heatmap_without_catalogue_table(self):
        matrix = [
            {
                "region_id": "region_1",
                "marker_order": 1,
                "snp_id": snp,
                "metabolite": metabolite,
                "beta": beta,
                "se": 0.1,
                "p": p,
                "n": 100,
                "conditioned_on_n": 0,
                "testable": "TRUE",
            }
            for snp, metabolite, beta, p in (
                ("rs1", "M1", 0.4, 1e-9),
                ("rs1", "M2", -0.1, 0.2),
                ("rs2", "M1", 0.2, 0.1),
                ("rs2", "M2", -0.5, 2e-10),
            )
        ]
        members = [
            {"gim_id": "region_1_GIM_001", "region_id": "region_1", "node_type": "SNP", "node_id": "rs1", "marker_order": 1},
            {"gim_id": "region_1_GIM_001", "region_id": "region_1", "node_type": "SNP", "node_id": "rs2", "marker_order": 2},
            {"gim_id": "region_1_GIM_001", "region_id": "region_1", "node_type": "metabolite", "node_id": "M1", "marker_order": ""},
            {"gim_id": "region_1_GIM_001", "region_id": "region_1", "node_type": "metabolite", "node_id": "M2", "marker_order": ""},
        ]
        summary = [
            {
                "gim_id": "region_1_GIM_001",
                "region_id": "region_1",
                "n_snps": 2,
                "n_metabolites": 2,
                "snps": "rs1;rs2",
                "metabolites": "M1;M2",
            }
        ]
        with tempfile.TemporaryDirectory(prefix="gimforge_report_test_") as directory:
            path = write_report(
                Path(directory) / "report.html",
                matrix_out=matrix,
                members=members,
                gim_summary=summary,
                conditional_p=5e-8,
                metadata={"LD reference ancestry": "EAS"},
            )
            report = path.read_text()
        self.assertIn("Retained GIM associations only", report)
        self.assertIn("Colour encodes conditional β", report)
        self.assertIn("LD reference ancestry</th><td>EAS", report)
        self.assertNotIn("<h2>All GIMs</h2>", report)

    def test_report_restores_retained_p_from_edges_for_legacy_matrix(self):
        matrix = [
            {
                "region_id": "R",
                "marker_order": 1,
                "snp_id": "rs1",
                "metabolite": "M1",
                "beta": 0.4,
                "p": "",
                "testable": "TRUE",
            }
        ]
        edges = [{**matrix[0], "p": 1e-9, "gim_id": "R_GIM_001"}]
        members = [
            {"gim_id": "R_GIM_001", "region_id": "R", "node_type": "SNP", "node_id": "rs1", "marker_order": 1},
            {"gim_id": "R_GIM_001", "region_id": "R", "node_type": "metabolite", "node_id": "M1", "marker_order": ""},
        ]
        summary = [{"gim_id": "R_GIM_001", "region_id": "R", "n_snps": 1, "n_metabolites": 1}]
        with tempfile.TemporaryDirectory(prefix="gimforge_report_test_") as directory:
            path = write_report(
                Path(directory) / "report.html",
                matrix_out=matrix,
                edges=edges,
                members=members,
                gim_summary=summary,
                conditional_p=5e-8,
            )
            report = path.read_text()
        self.assertIn('"p": 1e-09', report)


if __name__ == "__main__":
    unittest.main()
