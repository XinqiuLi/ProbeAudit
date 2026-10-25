ProbeAudit repository release 1.2.0
Reference-sequence auditing of lncRNA-associated microarray probes
Yuanning Zhang and Xinqiu Li
19 October 2026

The repository connects original probe labels to constituent sequence support,
question-specific rankings and frequency estimates. It contains all retained
inputs needed to reproduce the reported aggregation and statistical results.

REPRODUCE THE REPORTED RESULTS
Requires Python 3.10 or later. The aggregation, statistical and annotation
regeneration programs below use the Python standard library.
1. In legacy_annotation/: python regenerate_annotation.py
   Then: python verify_annotation_r4.py
2. In evidence/: python prioritize.py --out results
   Then: python -m unittest -v test_prioritize
3. From the repository root: python statistics/estimate_frequencies.py
The statistics program verifies all 2,209 selection hashes, reproduces the
100 selected probe IDs, counts unique sets and writes estimates and Wilson
intervals. The 15-set union includes 13 Ensembl sets and 11 RefSeq sets with
nine overlapping. Five have full alternative support and ten partial support.

AUDIT NEW PROBE SETS
1. Retrieve each original probe ID, complete annotation and eleven sequences.
2. Obtain references identified by releases and checksums. The retained engine
   and reference manifest are in legacy_annotation/retained_sequence_code/.
   Its engine imports and command-line help document alignment dependencies.
3. Generate original-orientation and reverse-complement transcript matches and
   contiguous genomic counts, documenting reference gene mappings. Unmodified
   RefSeq FASTA requires the appropriate accession-to-gene mapping.
4. Prepare the same input columns as evidence/inputs/, including a complete
   inventory, candidate labels, reference coverage and independently obtained
   summary counts. Empty group files retain their headers.
5. In evidence/ run python prioritize.py --inputs NEW_INPUTS --out NEW_RESULTS.
   Inspect report.html and the gene_support.csv/probe_metrics.csv tables.
Detailed input schemas, filter behavior and examples are in evidence/README.txt.

FILE GUIDE
legacy_annotation/: published selection and annotation-regeneration materials.
evidence/inputs/: all 114 sets and retained oligonucleotide-level matches.
evidence/results/: all gene-support combinations, metrics, complete rankings
and an offline report with expandable sequences, transcripts and locations.
statistics/: executable selection verification and frequency estimation.
CHANGELOG.txt: release scope and provenance.
SHA256SUMS.txt: checksums of every distributed file except this checksum list.

INTERPRETATION
The sample estimates sequence-pattern frequencies in the 2,209-set eligible
population under an outcome-unrelated hash-ordering approximation. Wilson
intervals use a binomial approximation without finite-population correction.
The historical 14 sets illustrate patterns. Rankings guide specific inspection
questions and were developed from the retained cases. Sequence support guides
attribution review; experimental specificity requires assay evidence.

SOURCE AND RIGHTS
Version 1.0.0 is archived at https://doi.org/10.5281/zenodo.23273760.
The release retains its licensing boundary: original code MIT; original
documentation and analytical outputs CC BY 4.0; third-party rights apply.
The large reference FASTA files remain with their original providers and are
identified through the manifests. AI assistance with code, computational
checks and documentation drafting is disclosed in the associated manuscript.
