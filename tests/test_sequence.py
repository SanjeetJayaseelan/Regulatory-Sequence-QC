import unittest

from regulatory_qc.models import QCConfig
from regulatory_qc.sequence import longest_homopolymer, normalized_kmer_entropy, sequence_checks


class SequenceChecksTests(unittest.TestCase):
    def setUp(self):
        self.config = QCConfig(
            min_length=10,
            max_length=20,
            min_gc=0.25,
            max_gc=0.75,
            max_homopolymer=4,
            repeat_k=3,
            max_repeat_count=2,
            max_repeat_fraction=0.8,
            min_complexity=0.2,
        )

    def test_length_and_alphabet_are_reported(self):
        _, issues = sequence_checks("ACGTN", self.config)
        codes = {issue.code for issue in issues}
        self.assertIn("invalid_alphabet", codes)
        self.assertIn("invalid_length", codes)

    def test_gc_and_homopolymer(self):
        metrics, issues = sequence_checks("AAAAACCCCC", self.config)
        self.assertEqual(metrics["length"], 10)
        self.assertAlmostEqual(metrics["gc_content"], 0.5)
        self.assertEqual(longest_homopolymer("AAAACG"), (4, "A"))
        self.assertIn("long_homopolymer", {issue.code for issue in issues})

    def test_complexity_is_bounded(self):
        self.assertEqual(normalized_kmer_entropy("AAAAAAAAAA", 2), 0.0)
        self.assertGreater(normalized_kmer_entropy("ACGTACGTAC", 2), 0.0)


if __name__ == "__main__":
    unittest.main()
