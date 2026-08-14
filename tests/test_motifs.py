import json
import tempfile
import unittest
from pathlib import Path

from regulatory_qc.models import MotifSpec
from regulatory_qc.motifs import load_motifs, reverse_complement, scan_motifs


class MotifTests(unittest.TestCase):
    def setUp(self):
        self.motif = MotifSpec(
            id="M1",
            matrix=(
                (0.95, 0.02, 0.02, 0.01),
                (0.02, 0.95, 0.02, 0.01),
                (0.02, 0.02, 0.95, 0.01),
            ),
            threshold=0.80,
        )

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ACGT"), "ACGT")
        self.assertEqual(reverse_complement("AAA"), "TTT")

    def test_scans_forward_and_reverse_strands(self):
        forward = scan_motifs("TTTACGTTT", [self.motif])
        self.assertTrue(any(hit.strand == "+" and hit.start == 3 for hit in forward))
        reverse = scan_motifs("TTTCGTAAAT", [self.motif])
        self.assertTrue(any(hit.strand == "-" and hit.start == 3 for hit in reverse))

    def test_unwanted_role_is_preserved(self):
        unwanted = MotifSpec(id="bad", matrix=self.motif.matrix, threshold=0.8, role="unwanted")
        hits = scan_motifs("ACG", [unwanted])
        self.assertEqual(hits[0].role, "unwanted")

    def test_condition_changes_effective_motif_role(self):
        conditional = MotifSpec(
            id="PTF1A",
            matrix=self.motif.matrix,
            threshold=0.8,
            conditions=("PTF1A_NEGATIVE",),
            unwanted_conditions=("KRAS_MAPK_ERK", "HNF4G_FOXA1", "GATA6"),
        )

        negative_control_hits = scan_motifs("ACG", [conditional], condition="PTF1A_NEGATIVE")
        pdac_hits = scan_motifs("ACG", [conditional], condition="GATA6")
        unrelated_hits = scan_motifs("ACG", [conditional], condition="UNKNOWN")

        self.assertEqual(negative_control_hits[0].role, "target")
        self.assertEqual(pdac_hits[0].role, "unwanted")
        self.assertEqual(unrelated_hits[0].role, "neutral")

    def test_loads_condition_and_source_metadata(self):
        payload = [{
            "id": "MA0001.1",
            "matrix": [[3, 1, 0, 0]],
            "matrix_type": "counts",
            "conditions": ["GATA6"],
            "unwanted_conditions": ["PTF1A_NEGATIVE"],
            "metadata": {"database": "JASPAR", "source_url": "https://example.test/motif"},
        }]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "motifs.json"
            path.write_text(json.dumps(payload))
            motif = load_motifs(path)[0]

        self.assertEqual(motif.conditions, ("GATA6",))
        self.assertEqual(motif.unwanted_conditions, ("PTF1A_NEGATIVE",))
        self.assertEqual(motif.metadata["database"], "JASPAR")

    def test_rejects_non_finite_matrix_values(self):
        with self.assertRaisesRegex(ValueError, "finite"):
            MotifSpec(id="bad", matrix=((float("nan"), 0.0, 0.0, 1.0),))

    def test_rejects_negative_matrix_values(self):
        with self.assertRaisesRegex(ValueError, "non-negative"):
            MotifSpec(id="bad", matrix=((-1.0, 1.0, 1.0, 1.0),))

    def test_rejects_non_finite_pseudocount(self):
        with self.assertRaisesRegex(ValueError, "pseudocount must be finite"):
            MotifSpec(
                id="bad",
                matrix=((1.0, 1.0, 1.0, 1.0),),
                matrix_type="counts",
                pseudocount=float("nan"),
            )

    def test_rejects_negative_pseudocount(self):
        with self.assertRaisesRegex(ValueError, "pseudocount must be non-negative"):
            MotifSpec(
                id="bad",
                matrix=((1.0, 1.0, 1.0, 1.0),),
                matrix_type="counts",
                pseudocount=-0.1,
            )

    def test_rejects_probability_rows_without_weight(self):
        with self.assertRaisesRegex(ValueError, "positive total"):
            MotifSpec(id="bad", matrix=((0.0, 0.0, 0.0, 0.0),))

    def test_rejects_count_rows_without_weight_or_pseudocount(self):
        with self.assertRaisesRegex(ValueError, "positive total"):
            MotifSpec(
                id="bad",
                matrix=((0.0, 0.0, 0.0, 0.0),),
                matrix_type="counts",
                pseudocount=0.0,
            )


if __name__ == "__main__":
    unittest.main()
