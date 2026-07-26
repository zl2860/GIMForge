import unittest

from gimforge.cli import build_parser


class CliTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
