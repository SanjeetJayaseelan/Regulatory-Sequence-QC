import unittest

from regulatory_qc.models import QCConfig


class ConfigValidationTests(unittest.TestCase):
    def test_complexity_k_must_be_positive(self):
        with self.assertRaisesRegex(ValueError, "complexity_k"):
            QCConfig(complexity_k=0)

    def test_motif_scan_step_must_be_positive(self):
        with self.assertRaisesRegex(ValueError, "motif_scan_step"):
            QCConfig(motif_scan_step=0)

    def test_score_weights_must_be_finite(self):
        with self.assertRaisesRegex(ValueError, "finite"):
            QCConfig(score_weights={"target_motif_strength": float("nan")})

    def test_score_weights_reject_unknown_component_names(self):
        with self.assertRaisesRegex(ValueError, "Unknown score weight"):
            QCConfig(score_weights={"target_motif_strenght": 1.0})


if __name__ == "__main__":
    unittest.main()
