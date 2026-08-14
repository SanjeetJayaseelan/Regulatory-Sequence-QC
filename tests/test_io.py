import json
import tempfile
import unittest
from pathlib import Path

from regulatory_qc.io import load_sequences


class InputValidationTests(unittest.TestCase):
    def test_fasta_header_requires_an_id(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.fasta"
            path.write_text(">\nACGT\n")

            with self.assertRaisesRegex(ValueError, "has no ID"):
                load_sequences(path)

    def test_csv_rejects_a_missing_sequence_value(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.csv"
            path.write_text("id,sequence\ncandidate_1\n")

            with self.assertRaisesRegex(ValueError, "missing a sequence"):
                load_sequences(path)

    def test_csv_rejects_extra_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.csv"
            path.write_text("id,sequence\ncandidate_1,ACGT,unexpected\n")

            with self.assertRaisesRegex(ValueError, "extra value"):
                load_sequences(path)

    def test_csv_rejects_a_missing_id_value(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.csv"
            path.write_text("id,sequence\n,ACGT\n")

            with self.assertRaisesRegex(ValueError, "missing an ID"):
                load_sequences(path)

    def test_json_rejects_a_null_sequence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.json"
            path.write_text(json.dumps([{"id": "candidate_1", "sequence": None}]))

            with self.assertRaisesRegex(ValueError, "missing a sequence"):
                load_sequences(path)


if __name__ == "__main__":
    unittest.main()
