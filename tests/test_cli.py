import io
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

from regulatory_qc.cli import main


class CommandLineTests(unittest.TestCase):
    def test_rejects_an_empty_candidate_file(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "empty.fasta"
            output_path = Path(directory) / "report.json"
            input_path.write_text("")

            with redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as raised:
                    main(["--input", str(input_path), "--output", str(output_path)])

            self.assertEqual(raised.exception.code, 2)
            self.assertIn("contains no sequences", stderr.getvalue())
            self.assertFalse(output_path.exists())

    def test_reports_invalid_input_without_a_traceback(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "invalid.fasta"
            output_path = Path(directory) / "report.json"
            input_path.write_text(">\nACGT\n")

            with redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as raised:
                    main(["--input", str(input_path), "--output", str(output_path)])

            self.assertEqual(raised.exception.code, 2)
            self.assertIn("has no ID", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
