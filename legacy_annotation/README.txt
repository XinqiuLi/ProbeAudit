Archived annotation regeneration and selection records

This directory retains the scientific inputs and programs from version 1.0.0,
https://doi.org/10.5281/zenodo.23273760. From this directory run:
python regenerate_annotation.py
python verify_annotation_r4.py

The first program rebuilds the 114-row annotation table; the second checks the
screening conditions independently. These operations use retained summaries.
The retained_sequence_code/ directory contains the sequence engine, reference
manifest and the scientific portions of the October 5 extension protocol.
Historical submission instructions have been omitted from the protocol copy;
the original published version remains available at the DOI above. Personal
absolute paths in selection_manifest.json have been reduced to basenames;
source hashes remain unchanged.

eligibility_frame.csv contains the 2,209 eligible sets and their fixed hashes;
selected_probes.csv contains the first 100. The qualifying reference headers
are provided for the 96 selected target symbols. Reference release identifiers
and external source checksums determine the relevant historical inputs.
