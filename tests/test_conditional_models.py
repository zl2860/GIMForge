import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from gimforge.conditional import (
    _BOLTConditionalRunner,
    _PLINKConditionalRunner,
    _read_bolt_stats,
    _read_glm,
)
from gimforge.parameters import parameters_from_args


class ConditionalModelTests(unittest.TestCase):
    def test_plink_primary_row_is_read_under_all_supported_encodings(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_model_test_") as directory:
            report = Path(directory) / "test.glm.linear"
            for genetic_model, test_label in (
                ("additive", "ADD"),
                ("dominant", "DOM"),
                ("recessive", "REC"),
            ):
                report.write_text(
                    "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE P ERRCODE\n"
                    f"1 100 rs1 A G G {test_label} 100 0.4 0.1 1e-9 .\n"
                )
                rows = _read_glm(report, genetic_model)
                self.assertEqual(rows[0]["snp_id"], "rs1")
                self.assertEqual(rows[0]["p"], 1e-9)

    def test_dominant_coding_is_applied_to_test_and_condition_variants(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_model_test_") as directory:
            scratch = Path(directory)
            runner = _PLINKConditionalRunner(
                bfile=scratch / "cohort",
                plink2=scratch / "plink2",
                phenotypes=scratch / "phenotypes.tsv",
                covariates=scratch / "covariates.tsv",
                keep_file=scratch / "keep.tsv",
                parameters=parameters_from_args({"genetic_model": "dominant"}),
                scratch=scratch,
                verbose=False,
            )
            with patch("gimforge.conditional.run_plink") as mocked:
                runner.scan(
                    chromosome="1",
                    start=1,
                    end=200,
                    metabolites=["M1"],
                    candidates=["rs1"],
                    conditioned_on=["rs2"],
                    label="test",
                )
            arguments = list(mocked.call_args.args[1])
            condition_index = arguments.index("--condition-list")
            self.assertEqual(arguments[condition_index + 2], "dominant")
            glm_index = arguments.index("--glm")
            self.assertIn("dominant", arguments[glm_index + 1 :])

    def test_bolt_infinitesimal_output_is_normalised(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_model_test_") as directory:
            report = Path(directory) / "bolt.stats"
            report.write_text(
                "SNP CHR BP ALLELE1 ALLELE0 BETA SE P_BOLT_LMM_INF\n"
                "rs1 1 100 G A 0.25 0.05 2e-10\n"
            )
            rows = _read_bolt_stats(report)
        self.assertEqual(rows[0]["snp_id"], "rs1")
        self.assertEqual(rows[0]["p"], 2e-10)
        self.assertEqual(rows[0]["beta"], 0.25)

    def test_bolt_runner_requests_infinitesimal_additive_statistic(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_model_test_") as directory:
            scratch = Path(directory)
            model = scratch / "model"
            model.with_suffix(".fam").write_text(
                "F1 I1 0 0 0 -9\nF2 I2 0 0 0 -9\nF3 I3 0 0 0 -9\n"
            )
            keep = scratch / "keep.tsv"
            keep.write_text("F1\tI1\nF2\tI2\nF3\tI3\n")
            covariates = scratch / "covariates.tsv"
            covariates.write_text("FID\tIID\tPC1\nF1\tI1\t0\nF2\tI2\t1\nF3\tI3\t-1\n")
            runner = _BOLTConditionalRunner(
                bfile=scratch / "cohort",
                model_bfile=model,
                plink2=scratch / "plink2",
                bolt=scratch / "bolt",
                genetic_map=scratch / "map.gz",
                phenotypes=scratch / "phenotypes.tsv",
                covariates=covariates,
                keep_file=keep,
                parameters=parameters_from_args({"regression_model": "mixed"}),
                scratch=scratch,
                verbose=False,
            )
            bgen, sample = scratch / "test.bgen", scratch / "test.sample"
            captured: list[str] = []

            def fake_bolt(_executable, arguments, *, context, quiet):
                del context, quiet
                captured.extend(map(str, arguments))
                stats = Path(
                    next(
                        argument.split("=", 1)[1]
                        for argument in captured
                        if argument.startswith("--statsFileBgenSnps=")
                    )
                )
                stats.write_text(
                    "SNP CHR BP ALLELE1 ALLELE0 BETA SE P_BOLT_LMM_INF\n"
                    "rs1 1 100 G A 0.25 0.05 2e-10\n"
                )

            with (
                patch.object(
                    runner,
                    "_candidate_bgen",
                    return_value=(bgen, sample, [bgen, sample]),
                ),
                patch.object(
                    runner,
                    "_condition_covariates",
                    return_value=(covariates, ["PC1"], []),
                ),
                patch("gimforge.conditional.run_bolt", side_effect=fake_bolt),
            ):
                rows = runner.scan(
                    chromosome="1",
                    start=1,
                    end=200,
                    metabolites=["M1"],
                    candidates=["rs1"],
                    conditioned_on=[],
                    label="test",
                )

        self.assertEqual(rows["M1"][0]["snp_id"], "rs1")
        self.assertIn("--lmmInfOnly", captured)
        self.assertIn("--qCovarCol=PC1", captured)
        self.assertIn("--bgenRefFirst", captured)


if __name__ == "__main__":
    unittest.main()
