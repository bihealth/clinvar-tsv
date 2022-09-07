"""Parse Clinvar XML

Code taken from MacArthur lab.

Reference on clinvar XML tag:

- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/clinvar_submission.xsd
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README
"""

import json
import xml.etree.ElementTree as ET

import binning
import cattr
from logzero import logger
import tqdm

from clinvar_tsv.common import ClinVarSet, DateTimeEncoder, as_pg_list

TSV_HEADER = "\t".join(
    (
        "release",
        "chromosome",
        "start",
        "end",
        "bin",
        "reference",
        "alternative",
        "variation_type",
        "symbols",
        "hgnc_ids",
        "vcv",
        "rcv",
        "review_status",
        "gold_stars",
        "pathogenicity",
        "origin",
        "details",
    )
)


class ClinvarParser:
    """Helper class for parsing Clinvar XML"""

    def __init__(self, input_file, out_b37, out_b38, max_rcvs=None):
        #: ``file``-like object to load the XML from
        self.input = input_file
        #: ``file``-like object to write GRCh37 variants to
        self.out_b37 = out_b37
        #: ``file``-like object to write GRCh38 variants to
        self.out_b38 = out_b38
        #: Number of processed RCVs, used for progress disable.
        self.rcvs = 0
        #: Largest number of rcvs to process out (for testing only)
        self.max_rcvs = max_rcvs

    def run(self):
        logger.info("Parsing elements...")
        out_files = {"GRCh37": self.out_b37, "GRCh38": self.out_b38}
        for out_file in out_files.values():
            print(TSV_HEADER, file=out_file)
        with tqdm.tqdm(unit="rcvs") as progress:
            for event, elem in ET.iterparse(self.input):
                if elem.tag == "ClinVarSet" and event == "end":
                    self.rcvs += 1
                    clinvar_set = ClinVarSet.from_element(elem)
                    if clinvar_set.ref_cv_assertion.observed_in:
                        origin = clinvar_set.ref_cv_assertion.observed_in.origin
                    else:
                        origin = "."
                    for genotype_set in clinvar_set.ref_cv_assertion.genotype_sets:
                        for measure_set in genotype_set.measure_sets:
                            for measure in measure_set.measures:
                                for build, location in measure.sequence_locations.items():
                                    if build not in out_files:
                                        continue
                                    elif location.ref is not None and location.alt is not None:
                                        if len(location.ref) == 1 and len(location.alt) == 1:
                                            variation_type = "snv"
                                        elif len(location.ref) == len(location.alt):
                                            variation_type = "mnv"
                                        else:
                                            variation_type = "indel"
                                        row = [
                                            build,
                                            location.chrom,
                                            location.start,
                                            location.stop,
                                            binning.assign_bin(location.start - 1, location.stop),
                                            location.ref,
                                            location.alt,
                                            variation_type,
                                            as_pg_list(measure.symbols),
                                            as_pg_list(measure.hgnc_ids),
                                            measure_set.accession,
                                            clinvar_set.ref_cv_assertion.clinvar_accession,
                                            clinvar_set.ref_cv_assertion.review_status,
                                            clinvar_set.ref_cv_assertion.gold_stars,
                                            clinvar_set.ref_cv_assertion.pathogenicity,
                                            origin,
                                            json.dumps(
                                                cattr.unstructure(clinvar_set), cls=DateTimeEncoder
                                            )
                                            .replace(r"\"", "'")
                                            .replace('"', '"""'),
                                        ]
                                        print("\t".join(map(str, row)), file=out_files[build])
                        progress.update()
                    elem.clear()
                if self.max_rcvs and self.rcvs >= self.max_rcvs:
                    logger.info("Breaking out after processing %d RCVs (as configured)", self.rcvs)
                    break
        logger.info("Done parsing elements")
