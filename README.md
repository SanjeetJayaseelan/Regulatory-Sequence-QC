# Regulatory Sequence QC

I built this project to check generated regulatory DNA sequences before they move into later analysis.


**What the pipeline checks
**
For each sequence, the pipeline checks:

- The sequence length
- Invalid DNA letters
- The percentage of G and C bases
- Long runs of the same base
- Repeated DNA patterns
- Exact and near-duplicate sequences
- Similarity to the training sequences
- Transcription-factor motif matches on both DNA strands

The pipeline also gives each candidate a score from 0 to 100. The report explains each part of that score.

** Running the pipeline
**
Use Python 3.10 or a newer version.

```bash
python -m regulatory_qc.cli \
  --input candidates.fasta \
  --motifs examples/motifs.json \
  --training training.fasta \
  --output results.json
```

The default candidate length is 200 to 500 base pairs. You can change the limits with command-line options.

The example sequences are shorter than real candidates. Use this command for the example files:

```bash
python -m regulatory_qc.cli \
  --input examples/candidates.fasta \
  --motifs examples/motifs.json \
  --output results.json \
  --min-length 4 \
  --max-length 20
```

**Sequence input
**
The pipeline accepts FASTA, CSV, and JSON files.

A CSV file must contain `id` and `sequence` columns. It can also contain a `condition` column.

```json
[
  {
    "id": "candidate_001",
    "sequence": "ACGT...",
    "condition": "GATA6"
  }
]
```

A FASTA header can include the condition and model information:

```text
>candidate_001|condition=GATA6|model=conditional
ACGT...
```

The current condition labels are:

- `KRAS_MAPK_ERK`
- `HNF4G_FOXA1`
- `GATA6`
- `PTF1A_NEGATIVE`

**Motif input
**
The motif file tells the pipeline which DNA patterns to find. Each matrix uses position rows with A, C, G, and T columns.

A motif can have one of three roles:

- `target`: The condition is expected to contain this motif.
- `unwanted`: The motif receives a score penalty.
- `neutral`: The pipeline reports the motif but does not change the score.

The `conditions` and `unwanted_conditions` fields control how the pipeline treats a motif for each candidate label.

## Output

The pipeline writes one JSON report. The report includes:

- Pass or fail status
- Sequence measurements
- Problems and warnings
- Motif locations and scores
- Similar sequences
- Score components
- A plain explanation of the final score

## Tests

Run the full test suite with:

```bash
python -m unittest discover -s tests -v
```

## Important limitation

A motif match only means that a DNA section resembles a known binding pattern. It does not prove binding or biological activity.

This pipeline is a quality-control and ranking tool. Laboratory experiments and validated activity models are still necessary for biological conclusions.
