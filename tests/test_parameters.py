import unittest

from gimforge.parameters import parameters_from_args


class ParameterTests(unittest.TestCase):
    def test_region_controls_are_independently_configurable(self):
        parameters = parameters_from_args(
            {
                "sentinel_p": 1e-10,
                "sentinel_clump_r2": 0.08,
                "sentinel_clump_window_kb": 750_000,
                "ld_span_r2": 0.12,
                "cross_metabolite_merge_r2": 0.7,
                "ld_window_kb": 900_000,
                "no_ld_half_width_kb": 400,
                "region_padding_kb": 200,
            }
        )
        self.assertEqual(parameters.sentinel_p, 1e-10)
        self.assertEqual(parameters.sentinel_clump_window_kb, 750_000)
        self.assertEqual(parameters.ld_span_r2, 0.12)
        self.assertEqual(parameters.region_padding_kb, 200)

    def test_cross_metabolite_merge_r2_must_be_stricter(self):
        with self.assertRaisesRegex(ValueError, "stricter"):
            parameters_from_args(
                {
                    "ld_span_r2": 0.2,
                    "cross_metabolite_merge_r2": 0.1,
                }
            )

    def test_linear_models_support_three_genotype_encodings(self):
        for genetic_model in ("additive", "dominant", "recessive"):
            parameters = parameters_from_args(
                {
                    "regression_model": "linear",
                    "genetic_model": genetic_model,
                }
            )
            self.assertEqual(parameters.genetic_model, genetic_model)

    def test_bolt_mixed_model_is_explicitly_additive_only(self):
        parameters = parameters_from_args(
            {
                "regression_model": "mixed",
                "genetic_model": "additive",
            }
        )
        self.assertEqual(parameters.mixed_backend, "bolt-lmm")
        with self.assertRaisesRegex(ValueError, "additive allele dosage only"):
            parameters_from_args(
                {
                    "regression_model": "mixed",
                    "genetic_model": "recessive",
                }
            )


if __name__ == "__main__":
    unittest.main()
