import json
import tempfile
import unittest
from pathlib import Path

from regulatory_qc.io import load_sequences, result_to_report
from regulatory_qc.models import MotifSpec, QCConfig, SequenceRecord
from regulatory_qc.pipeline import run_qc


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.config = QCConfig(
            min_length=6,
            max_length=30,
            min_gc=0.2,
            max_gc=0.8,
            max_homopolymer=5,
            repeat_k=3,
            max_repeat_count=5,
            max_repeat_fraction=1.0,
            min_complexity=0.0,
            duplicate_similarity=0.9,
        )
        self.motif = MotifSpec(
            id="target",
            matrix=((0.95, 0.02, 0.02, 0.01), (0.02, 0.95, 0.02, 0.01), (0.02, 0.02, 0.95, 0.01)),
            threshold=0.8,
        )

    def test_pipeline_returns_explanations_and_duplicates(self):
        records = [SequenceRecord("a", "ACGTACG"), SequenceRecord("b", "ACGTACG")]
        results, summary = run_qc(records, [self.motif], config=self.config)
        self.assertEqual(summary["candidate_count"], 2)
        self.assertEqual(summary["failed_count"], 2)
        self.assertTrue(any(issue.code == "near_duplicate" for issue in results[0].issues))
        self.assertTrue(results[0].score_explanation)
        self.assertIn("target_motif_strength", results[0].score_components)

    def test_report_is_json_serializable(self):
        results, _ = run_qc([SequenceRecord("a", "ACGTACG")], [self.motif], config=self.config)
        report = result_to_report(results, self.config, 1)
        json.dumps(report)

    def test_only_matching_condition_receives_target_motif_credit(self):
        conditional_motif = MotifSpec(
            id="conditional",
            matrix=self.motif.matrix,
            threshold=0.8,
            conditions=("KRAS_MAPK_ERK",),
        )
        records = [
            SequenceRecord("matching", "TTTACGTTT", {"condition": "KRAS_MAPK_ERK"}),
            SequenceRecord("other", "TTTACGAAA", {"condition": "GATA6"}),
        ]

        results, _ = run_qc(records, [conditional_motif], config=self.config)

        self.assertGreater(results[0].score_components["target_motif_strength"], 0.0)
        self.assertEqual(results[1].score_components["target_motif_strength"], 0.0)
        self.assertEqual(results[1].metrics["motifs"]["neutral_hit_count"], 1)

    def test_missing_condition_is_reported_for_condition_specific_motifs(self):
        conditional_motif = MotifSpec(
            id="conditional",
            matrix=self.motif.matrix,
            threshold=0.8,
            conditions=("GATA6",),
        )

        results, _ = run_qc(
            [SequenceRecord("missing-label", "TTTACGTTT")],
            [conditional_motif],
            config=self.config,
        )

        self.assertTrue(any(issue.code == "missing_condition" for issue in results[0].issues))

    def test_fasta_loader(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.fasta"
            path.write_text(">one|condition=GATA6|model=conditional description\nACGT\n>two\nTGCA\n")
            records = load_sequences(path)
        self.assertEqual([record.id for record in records], ["one", "two"])
        self.assertEqual(records[0].metadata["condition"], "GATA6")
        self.assertEqual(records[0].metadata["model"], "conditional")


if __name__ == "__main__":
    unittest.main()
