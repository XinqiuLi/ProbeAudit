ProbeAudit evidence explorer 1.1.0
Reference-sequence auditing of lncRNA-associated microarray probes

This directory extends the version 1.0.0 annotation archive (DOI 10.5281/zenodo.23273760) with per-gene evidence aggregation and user-selectable prioritization. It uses retained Ensembl 116 and RefSeq GRCh38.p14 RS_2025_08 matches. It aggregates retained matches into inspectable gene-support counts. Python 3.10 or later; standard library only for this extension. The separate archived alignment engine has additional dependencies.

REPRODUCTION
python prioritize.py --out results
python -m unittest -v test_prioritize

EXAMPLES
python prioritize.py --mode consensus --min-consensus 7 --out consensus7
python prioritize.py --mode ambiguity --reference ensembl --out ambiguity
python prioritize.py --mode orientation --out orientation
python prioritize.py --mode alternative --group extension --min-alternative 9 --out alternative9
python prioritize.py --mode deficit --reference refseq --out refseq_deficit
python prioritize.py --mode genome-multiple --out genomic_multiplicity
python prioritize.py --mode genome-missing --out genomic_absence
python prioritize.py --mode combined --weights alternative=2,deficit=1,consensus=1 --out custom
python prioritize.py --probe 225799_at --out CYTOR

All scores sort descending with probe ID as a deterministic tie-breaker. Tied cases have no evidence-based ordering. The modes were developed during inspection of the retained cases and are exploratory inspection tools. Combined mode is a user-defined weighted sum; the demonstration defaults to alternative+deficit+ambiguity, whose correlated components can count related evidence more than once. Weights and filters are recorded in run.json.

COUNTING
For each probe set, reference, orientation, mismatch cutoff and gene ID, count distinct oligo IDs. Multiple transcripts and positions do not multiply an oligo's support. Each oligo can support several genes; gene counts must not be added to count unique oligos. Retain gene IDs separately even if symbols coincide, including patch records. These are annotation-record counts, not independent biological loci. The intended group is the archived candidate symbol with LINC00152/CYTOR normalized as aliases; other symbols/unnamed IDs are alternative annotations. All intended identifiers are pooled by unique oligo ID. The retained archives do not provide a fully harmonized cross-reference gene-ID mapping. Cross-reference concordance therefore uses shared explicit symbols only and the intersection of supporting oligo IDs; a score of zero is not evidence that loci disagree. All reference-specific IDs remain available.

Alternative score: largest exact original-orientation count for any one alternative gene ID.
Deficit score: 11 minus the intended exact original-orientation count, excluding references with absent target records.
Ambiguity score: largest intersection of intended-supporting oligos and oligos supporting one alternative gene ID.
Consensus score: largest number of the same oligos supporting the same named alternative in both references (symbol unions within each reference, then intersection across references).
Orientation score: exact reverse-complement intended-supporting oligos not supported in original orientation.
Genome-multiple score: oligos with more than one contiguous genomic site.
Genome-missing score: oligos with no contiguous genomic site, kept separate from multiplicity because splicing can cause missing contiguous hits.
With --reference both, reference-specific scores use the larger value; this is an exploratory union view, not cross-reference agreement. Consensus always uses both references; genome scores use the single shared genomic audit. --max-intended uses the minimum across selected references, including missing-target cases; inspect target_represented before interpreting this filter. --min-alternative uses the maximum. --min-consensus is a separate two-reference requirement. Orientation '+' is the supplied sequence; '-' is its reverse complement, not a direct statement about assay direction.

FILES
inputs/*_hits.csv: original oligo-level transcript hit records, including sequence, reference IDs, orientation, mismatch count and transcript positions.
inputs/*_genome.csv: original oligo sequences and genomic site counts.
inputs/annotation.csv: required columns extracted from the prior 114-row annotation record for labels, coverage and regression checks.
input_manifest.json: archived source paths relative to the study project and SHA256 hashes. No personal absolute paths are required.
results/gene_support.csv: every observed gene-support combination, including intended and alternative genes, both orientations, exact and at-most-one-substitution evidence, supporting oligo and transcript IDs and available transcript locations. Zero-match genes have no hit row; intended zero counts are explicit in probe_metrics.csv.
results/probe_metrics.csv: all 114 probe sets and all ranking components, including intended support and alternative support separately for each reference.
results/support_distribution.csv: full 11-to-1 exact alternative-support distribution. One probe may occur in several support bins if it matches several genes. Probe-gene pairs and distinct probe sets are separate columns.
results/all_rankings.csv: all 114 ranks for every mode under the default both-reference view, with default weights for combined mode.
results/ranking.csv and report.html: the selected view; the self-contained HTML exposes sequences, oligo IDs, gene IDs, transcript IDs and locations without a server.

NEW DATA
Prepare the same named input files and pass --inputs DIRECTORY. Empty sample/reference files still require headers. Provide the complete eleven-oligo inventory for each set, complete hit tables, and an annotation.csv with candidate labels, reference coverage, reference counts and category columns matching the supplied schema. New annotation counts are independently checked against the hits on every run; mismatches fail explicitly. Do not append new sets only to a ranking table. Obtain new sequence matches using the archived sequence/genome engine and document reference versions and header mappings first. Unmodified RefSeq FASTA is not a drop-in engine input. This aggregation package contains the study-specific reference mappings as retained hit metadata, not a new reference download service.

INTERPRETATION
Complete alternative support with no intended support warrants immediate qualification of an exclusive intended-gene claim. Complete shared support warrants reporting ambiguity. Partial alternative support requires inspecting which oligos and loci are involved; it is not a vote for relabeling. Reference-specific overlaps and reference coverage gaps require annotation/locus review. Independent probes, updated reference checks and a sensitivity analysis excluding unresolved sets can inform analyses requiring reliable gene-level attribution. Experimental specificity and effects on biological conclusions require separate evidence.
