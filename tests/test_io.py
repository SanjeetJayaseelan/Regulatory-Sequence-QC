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

    def test_jsonl_assigns_ids_and_preserves_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "generated_batch.jsonl"
            records = [
                {"sequence": "ACGT", "condition": "control"},
                {"id": "", "sequence": "TGCA", "model_probability": 0.75},
                {"id": "provided_id", "sequence": "GATTACA"},
            ]
            path.write_text("\n".join(json.dumps(record) for record in records) + "\n")

            loaded = load_sequences(path)

        self.assertEqual(
            [record.id for record in loaded],
            ["generated_batch_000001", "generated_batch_000002", "provided_id"],
        )
        self.assertEqual(loaded[0].metadata["condition"], "control")
        self.assertEqual(loaded[1].metadata["model_probability"], 0.75)

    def test_jsonl_reports_the_line_number_for_invalid_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "generated_batch.jsonl"
            path.write_text('{"sequence": "ACGT"}\nnot-json\n')

            with self.assertRaisesRegex(ValueError, "line 2 contains invalid JSON"):
                load_sequences(path)

    def test_jsonl_requires_a_sequence_on_each_line(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "generated_batch.jsonl"
            path.write_text(json.dumps({"condition": "control"}) + "\n")

            with self.assertRaisesRegex(ValueError, "line 1 is missing a sequence"):
                load_sequences(path)

    def test_jsonl_requires_an_object_on_each_line(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "generated_batch.jsonl"
            path.write_text(json.dumps(["ACGT"]) + "\n")

            with self.assertRaisesRegex(ValueError, "line 1 must contain an object"):
                load_sequences(path)


if __name__ == "__main__":
    unittest.main()
