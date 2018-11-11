"""Parse Clinvar XML

Code taken from MacArthur lab.

Reference on clinvar XML tag:

- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/clinvar_submission.xsd
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README
"""

from collections import defaultdict
import logging
import re
import xml.etree.ElementTree as ET

from .exceptions import XmlParseException, NoSequenceLocation


#: Regular expression for grepping out PubMed mentions.  ``group(1)`` will be all the text after the word PubMed or
#: # PMID.
MENTIONS_PUBMED_REGEX = "(?:PubMed|PMID)(.*)"

#: Regular expression for grepping out PubMedIDs. ``group(1)`` will be the first PubMed ID, ``group(2)`` will be all
#: remaining text.
EXTRACT_PUBMED_ID_REGEX = "[^0-9]+([0-9]+)[^0-9](.*)"

#: Header to write to the TSV file.
HEADER = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "start",
    "stop",
    "strand",
    "variation_type",
    "variation_id",
    "rcv",
    "scv",
    "allele_id",
    "symbol",
    "hgvs_c",
    "hgvs_p",
    "molecular_consequence",
    "clinical_significance",
    "clinical_significance_ordered",
    "pathogenic",
    "likely_pathogenic",
    "uncertain_significance",
    "likely_benign",
    "benign",
    "review_status",
    "review_status_ordered",
    "last_evaluated",
    "all_submitters",
    "submitters_ordered",
    "all_traits",
    "all_pmids",
    "inheritance_modes",
    "age_of_onset",
    "prevalence",
    "disease_mechanism",
    "origin",
    "xrefs",
    "dates_ordered",
    "multi",
]


def replace_semicolons(s, replace_with=":"):
    """Replace semicolons with colons."""
    return s.replace(";", replace_with)


def remove_newlines_and_tabs(s):
    """Remove newline and tab characters."""
    return re.sub("[\t\n\r]", " ", s)


class ClinvarParser:
    """Helper class for parsing Clinvar XML"""

    def __init__(self, input_file, output_single, output_multi, genome_build, max_rows=None):
        #: ``file``-like object to load the XML from
        self.input = input_file
        #: ``file``-like object to write single variants to
        self.single = output_single
        #: ``file``-like object to write multi variants to (haplotypes, compound het.)
        self.multi = output_multi
        #: The genome build, one of ``"GRCh37"`` or ``"GRCh38"``.
        self.genome_build = genome_build
        #: Number of written rows, used for progress disable.
        self.written_rows = 0
        #: Counter for skipped items
        self.skipped_counter = {"missing SequenceLocation": 0}
        #: Counters for single
        self.single_counter = 0
        #: Counters for multi
        self.multi_counter = 0
        #: Largest number of rows to write out (for testing only)
        self.max_rows = max_rows

    def _write_row(self, output, row):
        if self.genome_build == "GRCh38":
            # Adjust contig name to be the one in hs38.
            if not row["chrom"].startswith("chr"):
                row["chrom"] = "chr" + row["chrom"]
            if row["chrom"] == "chrMT":
                row["chrom"] = "chrM"
        output.write("\t".join([str(row[field]) for field in HEADER]))
        output.write("\n")

    def _write_header(self, output):
        output.write("\t".join(HEADER))
        output.write("\n")

    def run(self):
        logging.info("Writing headers...")
        self._write_header(self.single)
        if self.multi.name != self.single.name and self.multi is not self.single:
            self._write_header(self.multi)
        logging.info("Parsing elements...")
        for event, elem in ET.iterparse(self.input):
            if elem.tag == "ClinVarSet" and event == "end":
                self._process_clinvar_set(elem)
            if self.max_rows and self.written_rows >= self.max_rows:
                logging.info(
                    "Breaking out after writing %d rows (as configured)", self.written_rows
                )
                break
        logging.info("Done parsing elements")

    def _process_clinvar_set(self, clinvar_set):
        """Process one ClinVarSet element

        Yields `(output_file, record)`, can yield multiple in case of comp. het. variants.
        """
        # Find the one ClinVarAccession tag corresponding to the reference accession RCV*
        rcv = clinvar_set.find("./ReferenceClinVarAssertion/ClinVarAccession")
        if rcv.attrib.get("Type") != "RCV":  # pragma: no cover
            raise XmlParseException(
                "Assumed to find RCV record but found type: %s" % rcv.attrib.get("Type")
            )

        row = {"rcv": rcv.attrib.get("Acc")}

        # There should only be one ReferenceClinVarAssertion
        ref_clinvar_assertions = clinvar_set.findall(".//ReferenceClinVarAssertion")
        if len(ref_clinvar_assertions) > 1:  # pragma: no cover
            raise XmlParseException("Assumed to find only one ReferenceClinVarAssertion")

        # Process all MeasureSet elements.  There can be multiple ones in the case of het. comp. variants.
        # A good example for this is RCV000201320.
        measure_sets = ref_clinvar_assertions[0].findall(".//MeasureSet")
        rows = []
        for measure_set in measure_sets:
            rows += list(self._process_measure_set({**row}, clinvar_set, measure_set))
        if rows:
            if len(rows) > 1:
                self.multi_counter += 1
                output_file = self.multi
            else:
                self.single_counter += 1
                output_file = self.single
            for row in rows:
                self._write_row(output_file, {**row, "multi": int(len(rows) > 1)})
            self.written_rows += 1
            logging.info(
                "%d entries completed %s, %d total",
                self.written_rows,
                ", ".join("%s skipped due to %s" % (v, k) for k, v in self.skipped_counter.items()),
                self.written_rows + sum(self.skipped_counter.values()),
            )

    def _process_measure_set(self, row, clinvar_set, measure_set):
        """Process the given measure set."""
        # Extend row with variation ID and type.
        row = {
            **row,
            "variation_id": measure_set.attrib.get("ID"),
            "variation_type": measure_set.attrib.get("Type"),
            "scv": ";".join(sorted(set(self._yield_scvs_in_clinvar_set(clinvar_set)))),
            "all_pmids": ";".join(sorted(set(self._yield_all_pmids_in_clinvar_set(clinvar_set)))),
            **self._find_submitters_in_clinvar_set(clinvar_set),
            **self._find_review_significance_last_evaluted_in_clinvar_set(clinvar_set),
            **self._find_ordered_review_status_and_significance(clinvar_set),
            **self._find_disease_infos(clinvar_set),
            "symbol": self._get_symbol_from_measure_set(measure_set),
        }
        measures = measure_set.findall(".//Measure")
        for i, measure in enumerate(measures):
            try:
                yield from self._process_measure(row, clinvar_set, measure_set, measure)
            except NoSequenceLocation:  # pragma: no cover
                pass  # swallow

    def _process_measure(self, row, clinvar_set, measure_set, measure):
        symbol = row["symbol"] or self._get_symbol_from_measure(measure)
        row = {
            **row,
            "symbol": symbol,
            "allele_id": measure.attrib.get("ID"),
            **self._get_genomic_location(measure, symbol),
            **self._get_molecular_consequence(measure),
        }
        yield row

    def _yield_all_pmids_in_clinvar_set(self, clinvar_set):
        # Find all the Citation nodes, and get the PMIDs out of them.
        for citation in clinvar_set.findall(".//Citation"):
            for id_node in citation.findall(".//ID"):
                if id_node.attrib.get("Source") == "PubMed":
                    yield id_node.text

        # Now find the Comment nodes, regex your way through the comments and extract anything that appears to be a
        # PMID.
        for comment in clinvar_set.findall(".//Comment"):
            mentions_pubmed = re.search(MENTIONS_PUBMED_REGEX, comment.text)
            if mentions_pubmed is not None and mentions_pubmed.group(1) is not None:
                remaining_text = mentions_pubmed.group(1)
                while True:
                    pubmed_id_extraction = re.search(EXTRACT_PUBMED_ID_REGEX, remaining_text)
                    if pubmed_id_extraction is None:
                        break
                    elif pubmed_id_extraction.group(1) is not None:
                        yield pubmed_id_extraction.group(1)
                        if pubmed_id_extraction.group(2) is not None:
                            remaining_text = pubmed_id_extraction.group(2)

    def _yield_scvs_in_clinvar_set(self, clinvar_set):
        # Find all SCV accession numbers from all additional ClinVarAssertion tags.
        scv_numbers = []
        for scv in clinvar_set.findall(".//ClinVarAssertion/ClinVarAccession"):
            if scv.attrib.get("Type") == "SCV":
                yield scv.attrib.get("Acc")

    def _find_submitters_in_clinvar_set(self, clinvar_set):
        # Now find any/all submitters.
        submitters_ordered = []
        for submitter_node in clinvar_set.findall(".//ClinVarSubmissionID"):
            if submitter_node.attrib is not None and "submitter" in submitter_node.attrib:
                submitters_ordered.append(submitter_node.attrib["submitter"].replace(";", ","))
        # Field all_submitters will get deduplicated while submitters_ordered won't
        return {
            "submitters_ordered": ";".join(submitters_ordered),
            "all_submitters": ";".join(set(submitters_ordered)),
        }

    def _find_review_significance_last_evaluted_in_clinvar_set(self, clinvar_set):
        result = {"review_status": "", "clinical_significance": "", "last_evaluated": "0000-00-00"}
        clinical_significance = clinvar_set.find(
            ".//ReferenceClinVarAssertion/ClinicalSignificance"
        )
        if clinical_significance.find(".//ReviewStatus") is not None:
            result["review_status"] = clinical_significance.find(".//ReviewStatus").text
        if clinical_significance.find(".//Description") is not None:
            result["clinical_significance"] = clinical_significance.find(".//Description").text
        if clinical_significance.attrib.get("DateLastEvaluated") is not None:
            result["last_evaluated"] = clinical_significance.attrib.get(
                "DateLastEvaluated", "0000-00-00"
            )
        return result

    def _find_ordered_review_status_and_significance(self, clinvar_set):
        """Get ordered version of review status and significance as well as the pathogenicity counts."""
        review_status_ordered = []
        clinical_significance_ordered = []
        for sig_elem in clinvar_set.findall(".//ClinVarAssertion/ClinicalSignificance"):
            review_status = sig_elem.find("./ReviewStatus")
            if review_status is not None:
                review_status_ordered.append(review_status.text)
            else:  # pragma: no cover
                review_status_ordered.append("")
            significance = sig_elem.find("./Description")
            if significance is not None:
                clinical_significance_ordered.append(significance.text)
            else:  # pragma: no cover
                clinical_significance_ordered.append("")
        assert len(review_status_ordered) == len(clinical_significance_ordered)
        dates_ordered = [
            x.attrib.get("DateLastEvaluated", "0000-00-00")
            for x in clinvar_set.findall(".//ClinVarAssertion/ClinicalSignificance")
            if x is not None
        ]
        keys = (
            "pathogenic",
            "likely_pathogenic",
            "uncertain_significance",
            "benign",
            "likely_benign",
        )
        return {
            "review_status_ordered": ";".join(review_status_ordered),
            "clinical_significance_ordered": ";".join(clinical_significance_ordered),
            **{key: clinical_significance_ordered.count(key) for key in keys},
            "dates_ordered": ";".join(dates_ordered),
        }

    @staticmethod
    def _disease_attribute_to_column(attribute_type):
        if attribute_type == "ModeOfInheritance":
            return "inheritance_modes"
        elif attribute_type in ("age of onset", "prevalence", "disease mechanism"):
            return attribute_type.replace(" ", "_")
        else:
            return None

    def _find_disease_infos(self, clinvar_set):
        result = {
            "all_traits": [],
            "inheritance_modes": [],
            "age_of_onset": [],
            "prevalence": [],
            "disease_mechanism": [],
            "origin": [],
            "xrefs": [],
        }
        for trait_set in clinvar_set.findall(".//TraitSet"):
            disease_name_nodes = trait_set.findall(".//Name/ElementValue")
            trait_values = [
                disease_name_node.text
                for disease_name_node in disease_name_nodes
                if disease_name_node.attrib is not None
                and disease_name_node.attrib.get("Type") == "Preferred"
            ]
            result["all_traits"] += trait_values

            for attribute_node in trait_set.findall(".//AttributeSet/Attribute"):
                attribute_type = attribute_node.attrib.get("Type")
                column_name = self.__class__._disease_attribute_to_column(attribute_type)
                column_value = attribute_node.text.strip()
                if column_name and column_value:
                    result[column_name].append(column_value)

            # Put all the cross references one column, it may contains NCBI gene ID, conditions ID in disease
            # databases.
            for xref_node in trait_set.findall(".//XRef"):
                xref_db = xref_node.attrib.get("DB")
                xref_id = xref_node.attrib.get("ID")
                result["xrefs"].append("%s:%s" % (xref_db, xref_id))
        for origin in clinvar_set.findall(".//ReferenceClinVarAssertion/ObservedIn/Sample/Origin"):
            result["origin"].append(origin.text)
        # Transform the result
        result = {
            key: value if key == "all_traits" else sorted(set(value))
            for key, value in result.items()
        }
        result = {key: ";".join(map(replace_semicolons, value)) for key, value in result.items()}
        return result

    def _get_symbol_from_measure_set(self, measure_set):
        result = ""
        var_name = measure_set.find(".//Name/ElementValue").text
        if var_name is not None:
            match = re.search(r"\(([A-Za-z0-9]+)\)", var_name)
            if match is not None:
                result = match.group(1)
        return result

    def _get_symbol_from_measure(self, measure):
        for symbol in measure.findall(".//Symbol") or []:
            if symbol.find("ElementValue").attrib.get("Type") == "Preferred":
                return symbol.find("ElementValue").text
        return ""

    def _get_genomic_location(self, measure, symbol):
        # Get SequenceLocation tag, if any.
        for sequence_location in measure.findall(".//SequenceLocation"):
            if sequence_location.attrib.get("Assembly") == self.genome_build:
                attribs = [
                    sequence_location.attrib.get(key)
                    for key in ("Chr", "start", "referenceAllele", "alternateAllele")
                ]
                if all([x is not None for x in attribs]):
                    genomic_location = sequence_location
                    # break after finding the first non-empty locating matching our genome build
                    break
        else:  # pragma: no cover
            self.skipped_counter["missing SequenceLocation"] += 1
            logging.debug("Skipping Measure because it has no location")
            raise NoSequenceLocation("No sequence location in Measure")

        # Get strand value
        strand = ""
        for measure_relationship in measure.findall(".//MeasureRelationship"):
            if measure_relationship.find(".//Symbol/ElementValue").text != symbol:
                continue
            for sequence_location in measure_relationship.findall(".//SequenceLocation"):
                if (
                    "Strand" in sequence_location.attrib
                    and genomic_location.attrib["Accession"]
                    == sequence_location.attrib["Accession"]
                ):
                    strand = sequence_location.attrib["Strand"]
                    break

        # Build result
        return {
            "chrom": genomic_location.attrib["Chr"],
            "pos": genomic_location.attrib["start"],
            "ref": genomic_location.attrib["referenceAllele"],
            "alt": genomic_location.attrib["alternateAllele"],
            "start": genomic_location.attrib["start"],
            "stop": genomic_location.attrib["stop"],
            "strand": strand,
        }

    def _get_molecular_consequence(self, measure):
        molecular_consequence = set()
        hgvs_c = ""
        hgvs_p = ""
        attribute_set = measure.findall("./AttributeSet")
        for attribute_node in attribute_set:
            attribute_type = attribute_node.find("./Attribute").attrib.get("Type")
            attribute_value = attribute_node.find("./Attribute").text
            # Find hgvs_c
            if attribute_type == "HGVS, coding, RefSeq" and "c." in attribute_value:
                hgvs_c = attribute_value
            # Find hgvs_p
            if attribute_type == "HGVS, protein, RefSeq" and "p." in attribute_value:
                hgvs_p = attribute_value
            # Aggregate all molecular consequences
            if attribute_type == "MolecularConsequence":
                for xref in attribute_node.findall(".//XRef"):
                    if xref.attrib.get("DB") == "RefSeq":
                        # print xref.attrib.get('ID'), attribute_value
                        molecular_consequence.add(
                            ":".join([xref.attrib.get("ID"), attribute_value])
                        )
        return {
            "molecular_consequence": remove_newlines_and_tabs(
                ";".join(map(replace_semicolons, molecular_consequence))
            ),
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
        }
