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


def try_map(values, converter, catch_exc):
    """Try to convert the given iterable of converters ignoring exceptions up to the last one"""
    for i, value in enumerate(values):
        try:
            return converter(value)
        except catch_exc:
            if i + 1 == len(values):
                raise  # re-raise


class ClinvarParser:
    """Helper class for parsing Clinvar XML"""

    def __init__(
        self, input_file, out_b37_small, out_b37_sv, out_b38_small, out_b38_sv, max_rcvs=None
    ):
        #: ``file``-like object to load the XML from
        self.input = input_file
        #: ``file``-like object to write GRCh37 small variants to
        self.out_b37_small = out_b37_small
        #: ``file``-like object to write GRCh37 SVs to
        self.out_b37_sv = out_b37_sv
        #: ``file``-like object to write GRCh38 small variants to
        self.out_b38_small = out_b38_small
        #: ``file``-like object to write GRCh38 SVs to
        self.out_b38_sv = out_b38_sv
        #: Number of processed RCVs, used for progress disable.
        self.rcvs = 0
        #: Largest number of rcvs to process out (for testing only)
        self.max_rcvs = max_rcvs

    def run(self):  # noqa: C901
        logger.info("Parsing elements...")
        out_files = {
            "GRCh37": {"small": self.out_b37_small, "sv": self.out_b37_sv},
            "GRCh38": {"small": self.out_b38_small, "sv": self.out_b38_sv},
        }
        for d in out_files.values():
            for out_file in d.values():
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
                                    if location.ref is None or location.alt is None:
                                        if measure.measure_type in (
                                            "single nucleotide variant",
                                            "Variation",
                                        ):
                                            continue  # skip, just a small variant without a coordinate
                                        try:
                                            start = try_map(
                                                (
                                                    location.start,
                                                    location.outer_start,
                                                    location.inner_start,
                                                ),
                                                int,
                                                TypeError,
                                            )
                                            stop = try_map(
                                                (
                                                    location.stop,
                                                    location.outer_stop,
                                                    location.inner_stop,
                                                ),
                                                int,
                                                TypeError,
                                            )
                                        except TypeError:
                                            logger.info(
                                                "Cannot determine location from %s", measure
                                            )
                                            continue
                                        row = [
                                            build,
                                            location.chrom,
                                            start,
                                            stop,
                                            binning.assign_bin(start - 1, stop),
                                            location.ref or ".",
                                            location.alt or ".",
                                            measure.measure_type,
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
                                        print("\t".join(map(str, row)), file=out_files[build]["sv"])
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
                                        print(
                                            "\t".join(map(str, row)), file=out_files[build]["small"]
                                        )
                        progress.update()
                    elem.clear()
                if self.max_rcvs and self.rcvs >= self.max_rcvs:
                    logger.info("Breaking out after processing %d RCVs (as configured)", self.rcvs)
                    break
        logger.info("Done parsing elements")
