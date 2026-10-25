# ProbeAudit

ProbeAudit connects legacy microarray probe labels to constituent sequence evidence. It distinguishes absent intended support, shared intended and alternative support, and orientation or genomic-placement questions, then ranks probe sets for targeted investigation.

This repository contains the retained code, documentation, data and analyses for **Sequence auditing reveals alternative-gene support in apparently straightforward lncRNA microarray assignments**, by Yuanning Zhang and Xinqiu Li.

The fixed version 1.2.0 archive is available at [Zenodo](https://doi.org/10.5281/zenodo.23473178). The original 46 archive files are preserved here; this README and the ignore file provide repository navigation.

## Findings and uses

- In the additional sample of 100 sets, 15 have exact alternative-gene matches: five full-support and ten partial-support cases.
- CEP170P1 and CYTOR illustrate different attribution problems: absent intended support versus complete shared support.
- PVT1 demonstrates why orientation should be inspected separately from alternative support.
- Per-gene records distinguish support concentrated on one alternative from support dispersed across several identifiers.

## Reproduce the retained results

Python 3.10 or later is required. Aggregation, annotation regeneration and statistics use the Python standard library.

```sh
cd legacy_annotation
python regenerate_annotation.py
python verify_annotation_r4.py
cd ../evidence
python prioritize.py --out results
python -m unittest -v test_prioritize
cd ..
python statistics/estimate_frequencies.py
```

Open `evidence/results/report.html` in a browser to inspect the offline report. The report links rankings to oligonucleotide sequences, gene and transcript identifiers, and available locations.

## Choose an investigation

Run these commands in `evidence/`:

```sh
python prioritize.py --mode consensus --min-consensus 7 --out consensus7
python prioritize.py --mode ambiguity --reference ensembl --out ambiguity
python prioritize.py --mode orientation --out orientation
python prioritize.py --probe 225799_at --out CYTOR
```

Alternative support highlights competing assignments; deficit highlights incomplete intended support; ambiguity highlights shared compatibility; consensus highlights named alternatives supported across references; orientation highlights reverse-only intended matches. Genomic components separate multiple locations from missing contiguous sites. The manuscript defines the scores and compares their purposes; the report preserves component values, filters and weights.

## Audit new sets

New sets require complete eleven-oligonucleotide inventories, documented reference releases and gene mappings, transcript-hit tables and genomic counts matching the supplied schema. The retained sequence engine is in `legacy_annotation/retained_sequence_code/`; its dependencies are separate from the standard-library aggregation. Unmodified RefSeq FASTA requires accession-to-gene mapping. See [the complete workflow](README.txt) and [input schemas and ranking options](evidence/README.txt).

## Repository contents

- `legacy_annotation/`: sequence engine, reference manifests, eligibility frame, selection records and annotation-regeneration programs.
- `evidence/inputs/`: retained constituent matches and inventories for 114 sets.
- `evidence/results/`: gene support, set metrics, full rankings and the interactive report.
- `statistics/`: selection verification and frequency estimates.
- `SHA256SUMS.txt`: checksums for the original archived files.

## Scope and rights

Frequency estimates concern the defined 2,209-set eligibility frame under the outcome-unrelated hash-ordering approximation. Rankings were developed from the retained cases and guide exploration. Sequence compatibility supports attribution review; experimental specificity requires assay evidence.

Original code is MIT licensed; original documentation and analytical outputs are CC BY 4.0; third-party rights apply. Reference FASTA files remain with their providers and are identified by the archived manifests.
