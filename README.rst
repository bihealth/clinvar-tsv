.. image:: https://coveralls.io/repos/github/bihealth/clinvar-tsv/badge.svg?branch=main
   :target: https://coveralls.io/github/bihealth/clinvar-tsv?branch=main

===========
Clinvar-TSV
===========

The code in this repository allows to first download,t hen convert ClinVar XML files into TSV files (one for b37 and b38).
The TSV files will contain one entry for each ClinVar ``<ReferenceClinVarAssertion>`` entry with important information extracted from ClinVar.
The code is used by `bihealth/varfish-db-downloader <https://github.com/bihealth/varfish-db-downloader/>`_.

- `clinvar-tsv on PyPi <https://pypi.org/project/clinvar-tsv/>`__
- `clinvar-tsv on bioconda <http://bioconda.github.io/recipes/clinvar-tsv/README.html>`__

--------
Overview
--------

Users usually run the tool by calling ``clinvar_tsv main``.

::

    clinvar_tsv main \
        --cores 2 \
        --b37-path hs37d5.fa \
        --b38-path hs38.fa

This will call a Snakemake workflow that will in turn do the following

1. Download the latest ClinVar XML file to the ``downloads/`` directory using ``wget``.
2. Parse the XML file and convert it into a "raw" TSV file in ``parsed`` for each the 37 and 38 release with ``clinvar_tsv parse_xml``.
   This file contains one record for each ClinVar VCV record.
3. Sort this file by coordinate and VCV ID using Unix ``sort``, and finally...
4. Merge the lines in the resulting TSV file (for each genome build) by VCV ID and produce aggregate summaries for each VCV.

There are two summaries:

- ``summary_clinvar_*`` -- which merges record which attempts to imitate the `approach taken by ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/>`__
- ``summary_paranoid_*`` -- which considers all assessment as equally important, whether the reporter provided assessment criteria or not

----------
References
----------

Documentation in ClinVar:

- https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
- https://www.ncbi.nlm.nih.gov/clinvar/docs/help/
- https://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/
