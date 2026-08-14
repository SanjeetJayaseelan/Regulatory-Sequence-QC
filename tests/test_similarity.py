import unittest

from regulatory_qc.similarity import find_batch_duplicates, sequence_similarity


class SimilarityTests(unittest.TestCase):
    def test_equal_length_similarity_uses_hamming_identity(self):
        similarity, method = sequence_similarity("ACGTACGT", "ACGTACGA")
        self.assertEqual(method, "hamming_identity")
        self.assertAlmostEqual(similarity, 0.875)

    def test_different_lengths_use_kmer_jaccard(self):
        similarity, method = sequence_similarity("ACGTACGT", "ACGTACGTAA")
        self.assertEqual(method, "kmer_jaccard")
        self.assertGreater(similarity, 0.0)

    def test_batch_matches_are_symmetric(self):
        matches = find_batch_duplicates([("a", "ACGT"), ("b", "ACGA")], 0.7, 2)
        self.assertEqual(matches["a"][0].other_id, "b")
        self.assertEqual(matches["b"][0].other_id, "a")


if __name__ == "__main__":
    unittest.main()
