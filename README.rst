===========
Clinvar-TSV
===========

This project allows to download ClinVar files and convert them into easy-to-use TSV files.
The project is based on the from by `macarthurlab/clinvar <https://github.com/macarthur-lab/clinvar>`_.
As the latter is Python 2 only and unmaintained, we decided tot reimplement large parts of it using Python 3.

The code in this repository differs to the original code mainly in:

- Makes available a Python executable `clinvar_tsv` that performs the download and conversion.
- Workflow execution is done using Snakemake.
- No by-allele grouping is done.
- No annotation using ExAC or gnomAD is done.
- It is guaranteed that the clinical significance array has the same length as the review status array.
- The code has been refactored to consist of more smaller functions which hopefully improves maintainability.
