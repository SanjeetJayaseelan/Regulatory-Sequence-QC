# Regulatory Sequence QC

This project checks generated regulatory DNA sequences before they move into later analysis. It applies the same quality-control rules to every candidate and produces a JSON report with measurements, warnings, blocking errors, motif matches, similarity results, and a score from 0 to 100.

## What the pipeline checks

- Sequence length and invalid DNA letters
- GC content
- Long homopolymers
- Repeated DNA patterns and low complexity
- Exact and near-duplicate sequences
- Similarity to a supplied training set
- Transcription-factor motif matches on both DNA strands

## Requirements and installation

Use Python 3.10 or newer.

```bash
git clone https://github.com/SanjeetJayaseelan/Regulatory-Sequence-QC.git
cd Regulatory-Sequence-QC
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

## Try the included example

The example sequences are intentionally shorter than real candidates, so this command changes the length limits:

```bash
regulatory-qc \
  --input examples/candidates.fasta \
  --motifs examples/motifs.json \
  --output results.json \
  --min-length 4 \
  --max-length 20
```

For real project files, the command will look like this:

```bash
regulatory-qc \
  --input generated_sequences.fasta \
  --motifs validated_motifs.json \
  --training training_sequences.fasta \
  --output results.json
```

The training file and motif file are optional. Without them, the pipeline still performs the sequence-level checks.

## Sequence input

The pipeline accepts FASTA, CSV, and JSON files.

FASTA headers can include condition and model information:

```text
>candidate_001|condition=GATA6|model=conditional
ACGTACGTACGT
```

A CSV file must contain `id` and `sequence` columns. The `condition` column is optional:

```csv
id,sequence,condition
candidate_001,ACGTACGTACGT,GATA6
```

A JSON file can contain a list of records:

```json
[
  {
    "id": "candidate_001",
    "sequence": "ACGTACGTACGT",
    "condition": "GATA6"
  }
]
```

The current condition labels are:

- `KRAS_MAPK_ERK`
- `HNF4G_FOXA1`
- `GATA6`
- `PTF1A_NEGATIVE`

## Motif input

The motif file is JSON. Each position-weight matrix uses one row per DNA position and four columns in A, C, G, T order.

```json
[
  {
    "id": "GATA6_example",
    "name": "Example GATA6 target motif",
    "matrix_type": "probabilities",
    "matrix": [
      [0.90, 0.03, 0.03, 0.04],
      [0.05, 0.80, 0.10, 0.05]
    ],
    "threshold": 0.80,
    "role": "target",
    "conditions": ["GATA6"]
  }
]
```

A motif can have one of three roles:

- `target`: expected for the candidate's condition
- `unwanted`: reported with a score penalty
- `neutral`: reported without changing the score

The `conditions` and `unwanted_conditions` fields control how a motif is interpreted for each candidate label. See [examples/motifs.json](examples/motifs.json) for a complete example.

## Current default thresholds

The software currently uses:

- Length: 200–500 base pairs
- GC content: 30%–70%
- Maximum homopolymer: 8 bases
- Repeat size: 8 bases
- Maximum repeat count: 3
- Maximum repeated fraction: 35%
- Near-duplicate similarity: 95%

These are software defaults, not final biological claims. Replace them with the team's approved thresholds when those values are finalized.

## Output

The pipeline writes one JSON report containing:

- Pass or fail status
- Sequence measurements
- Problems and warnings
- Motif locations and scores
- Similar sequences
- Score components
- A plain-language explanation of the final score

## Tests

Run the complete test suite with:

```bash
python -m unittest discover -s tests -v
```

GitHub also runs the tests automatically on supported Python versions whenever code is pushed or a pull request is opened.

## Important limitation

A motif match only means that part of a sequence resembles a supplied binding pattern. It does not prove transcription-factor binding or biological activity.

This pipeline is a quality-control and ranking tool. Final conclusions still require validated motif data, approved thresholds, appropriate reference sequences, and laboratory testing.
