import tempfile
import unittest
from pathlib import Path

from gimforge.conditional import _validated_keep


class IndividualTableTests(unittest.TestCase):
    def test_tabular_individual_inputs_are_matched_by_fid_iid(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_input_test_") as directory:
            root = Path(directory)
            phenotype = root / "phenotypes.tsv"
            covariate = root / "covariates.tsv"
            phenotype.write_text("FID\tIID\tM001\nF1\tI1\t0.4\nF1\tI2\tNA\nF1\tI3\t-0.2\n")
            covariate.write_text("FID\tIID\tage\tPC1\nF1\tI1\t50\t0.1\nF1\tI2\t60\t-0.1\nF1\tI3\t55\t0.2\n")
            (root / "cohort.fam").write_text(
                "F1 I1 0 0 0 -9\nF1 I2 0 0 0 -9\nF1 I3 0 0 0 -9\n"
            )
            keep, traits = _validated_keep(
                phenotypes=phenotype,
                covariates=covariate,
                bfile=root / "cohort",
                scratch=root,
            )
            self.assertEqual(traits, ["M001"])
            self.assertEqual(keep.read_text(), "F1\tI1\nF1\tI2\nF1\tI3\n")

    def test_csv_phenotype_is_rejected_before_plink(self):
        with tempfile.TemporaryDirectory(prefix="gimforge_input_test_") as directory:
            root = Path(directory)
            phenotype = root / "phenotypes.csv"
            covariate = root / "covariates.tsv"
            phenotype.write_text("FID,IID,M001\nF1,I1,0.4\n")
            covariate.write_text("FID\tIID\tage\nF1\tI1\t50\n")
            (root / "cohort.fam").write_text("F1 I1 0 0 0 -9\n")
            with self.assertRaisesRegex(ValueError, "tab-separated"):
                _validated_keep(
                    phenotypes=phenotype,
                    covariates=covariate,
                    bfile=root / "cohort",
                    scratch=root,
                )


if __name__ == "__main__":
    unittest.main()
