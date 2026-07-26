import unittest

from gimforge.cli import build_parser


class CliTests(unittest.TestCase):
    def test_gim_expansion_is_metabotype(self):
        self.assertIn("genetically influenced metabotypes", build_parser().description.lower())

    def test_clump_subcommand_exposes_reference_and_clump_controls(self):
        args = build_parser().parse_args(
            [
                "clump",
                "--sumstats",
                "mgwas.tsv.gz",
                "--ld-bfile",
                "reference/EAS",
                "--ancestry",
                "EAS",
                "--out",
                "clump-result",
                "--sentinel-p",
                "1e-10",
                "--sentinel-clump-r2",
                "0.2",
                "--sentinel-clump-window-kb",
                "5000",
            ]
        )
        self.assertEqual(args.command, "clump")
        self.assertEqual(args.ancestry, "EAS")
        self.assertEqual(args.sentinel_p, 1e-10)
        self.assertEqual(args.sentinel_clump_r2, 0.2)
        self.assertEqual(args.sentinel_clump_window_kb, 5000)

    def test_manifest_can_replace_stacked_sumstats(self):
        args = build_parser().parse_args(
            [
                "clump",
                "--sumstats-manifest",
                "mgwas_files.tsv",
                "--ld-bfile",
                "reference/EAS",
                "--ancestry",
                "EAS",
                "--out",
                "clump-result",
            ]
        )
        self.assertIsNone(args.sumstats)
        self.assertEqual(args.sumstats_manifest, "mgwas_files.tsv")

    def test_run_can_resume_from_trait_specific_sentinels(self):
        args = build_parser().parse_args(
            [
                "run",
                "--sentinels",
                "clumping/sentinels.tsv",
                "--bfile",
                "data/cohort",
                "--ld-bfile",
                "reference/EAS",
                "--ancestry",
                "EAS",
                "--phenotypes",
                "data/phenotypes.tsv",
                "--covariates",
                "data/covariates.tsv",
                "--out",
                "result",
            ]
        )
        self.assertIsNone(args.sumstats)
        self.assertIsNone(args.sumstats_manifest)
        self.assertEqual(args.sentinels, "clumping/sentinels.tsv")

    def test_run_exposes_genetic_and_regression_models(self):
        args = build_parser().parse_args(
            [
                "run",
                "--sentinels",
                "clumping/sentinels.tsv",
                "--bfile",
                "data/cohort",
                "--ld-bfile",
                "reference/EAS",
                "--ancestry",
                "EAS",
                "--phenotypes",
                "data/phenotypes.tsv",
                "--covariates",
                "data/covariates.tsv",
                "--out",
                "result",
                "--genetic-model",
                "dominant",
                "--regression-model",
                "linear",
            ]
        )
        self.assertEqual(args.genetic_model, "dominant")
        self.assertEqual(args.regression_model, "linear")

    def test_mixed_run_accepts_bolt_inputs(self):
        args = build_parser().parse_args(
            [
                "run",
                "--sentinels",
                "clumping/sentinels.tsv",
                "--bfile",
                "data/cohort",
                "--ld-bfile",
                "reference/EAS",
                "--ancestry",
                "EAS",
                "--phenotypes",
                "data/phenotypes.tsv",
                "--covariates",
                "data/covariates.tsv",
                "--out",
                "result",
                "--regression-model",
                "mixed",
                "--bolt",
                "/opt/bolt/bolt",
                "--bolt-model-bfile",
                "data/cohort",
                "--bolt-genetic-map",
                "reference/genetic_map_hg38.txt.gz",
            ]
        )
        self.assertEqual(args.regression_model, "mixed")
        self.assertEqual(args.bolt, "/opt/bolt/bolt")
        self.assertEqual(args.bolt_model_bfile, "data/cohort")


if __name__ == "__main__":
    unittest.main()
