"""Parse Clinvar XML

Code taken from MacArthur lab.

Reference on clinvar XML tag:

- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/clinvar_submission.xsd
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README
"""

import datetime
import json
from uuid import UUID, uuid4
import typing
import re
import xml.etree.ElementTree as ET
from dateutil.parser import isoparse
from itertools import chain

from .exceptions import XmlParseException, NoSequenceLocation

import attr
import cattr
from logzero import logger


class UUIDEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, UUID):
            # if the obj is uuid, we simply return the value of uuid
            return obj.hex
        elif isinstance(obj, datetime.date):
            return obj.strftime("%Y-%m-%d")
        else:
            return json.JSONEncoder.default(self, obj)


#: Regular expression for grepping out PubMed mentions.  ``group(1)`` will be all the text after the word PubMed or
#: # PMID.
MENTIONS_PUBMED_REGEX = "(?:PubMed|PMID)(.*)"

#: Regular expression for grepping out PubMedIDs. ``group(1)`` will be the first PubMed ID, ``group(2)`` will be all
#: remaining text.
EXTRACT_PUBMED_ID_REGEX = "[^0-9]+([0-9]{1,10})[^0-9](.*)"

#: Header to write to the TSV file.
HEADER = [
    "release",
    "chromosome",
    "position",
    "reference",
    "alternative",
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


def convert_postgres_tsv(val):
    """Convert value into PostgreSQL TSV value."""

    def escape_str(x):
        if isinstance(x, str):
            return '"""%s"""' % x
        elif x is None:
            return "NULL"
        else:
            return str(x)

    if isinstance(val, list):
        return "{%s}" % ",".join([escape_str(v) for v in val])
    elif val is None:
        return ""
    else:
        return str(val)


TClinicalSignificance = typing.TypeVar("ClinicalSignificance")
TReferenceClinVarAssertion = typing.TypeVar("ReferenceClinVarAssertion")
TClinVarSet = typing.TypeVar("ClinVarSet")


@attr.s(frozen=True, auto_attribs=True)
class ClinicalSignificance:
    """Represent clinical significance."""

    #: Date of last evaluation.
    date_evaluated: datetime.date
    #: Significance review status.
    review_status: str
    #: Significance description.
    description: str
    #: Comments.
    comments: typing.Tuple[str, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @staticmethod
    def from_element(element: ET.Element) -> TClinicalSignificance:
        return ClinicalSignificance(
            date_evaluated=isoparse(element.attrib.get("DateLastEvaluated")),
            review_status=element.find("ReviewStatus").text,
            description=element.find("Description").text,
            comments=tuple(elem.text for elem in element.findall("./Comment")),
        )


@attr.s(frozen=True, auto_attribs=True)
class ObservedDataDescription:
    """Relevant information from ObservedData[@Type='Description']."""

    #: Optional description text.
    description: typing.Optional[str]
    #: PubMed IDs.
    pubmed_ids: typing.Tuple[int]
    #: OMIM IDs.
    omim_ids: typing.Tuple[int]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @staticmethod
    def from_element(element: ET.Element):
        description = element.find("./Attribute[@Type='Description']").text
        if description == "not provided":
            return None
        else:
            pubmed_ids = [
                int(elem.text) for elem in element.findall("./Citation/ID[@Source='PubMed']")
            ]
            omim_ids = [
                int(elem.attrib.get("ID")) for elem in element.findall("./XRef[@Type='MIM']")
            ]
            return ObservedDataDescription(
                description=description, pubmed_ids=tuple(pubmed_ids), omim_ids=tuple(omim_ids),
            )


@attr.s(frozen=True, auto_attribs=True)
class ObservedIn:
    """Relevant part of ObservedIn."""

    #: The origin of the sample.
    origin: str
    #: The species of the sample.
    species: str
    #: Affected state
    affected_status: str
    #: Optional observation info.
    observed_data_description: typing.Optional[ObservedDataDescription]
    #: Comments.
    comments: typing.Tuple[str, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        if element is None:
            return None
        for elem in element.findall("./ObservedData"):
            if elem.find("./Attribute[@Type='Description']") is not None:
                observed_data_description = ObservedDataDescription.from_element(elem)
                break
        else:  # no break above
            observed_data_description = None
        return ObservedIn(
            origin=element.find("./Sample/Origin").text,
            species=element.find("./Sample/Species").text,
            affected_status=element.find("./Sample/AffectedStatus").text,
            observed_data_description=observed_data_description,
            comments=tuple(elem.text for elem in element.findall("./Comment")),
        )


def mapply(f, x):
    if x is None:
        return x
    else:
        return f(x)


@attr.s(frozen=True, auto_attribs=True)
class SequenceLocation:
    """The relevant information from a SequenceLocation."""

    assembly: str
    chrom: str
    chrom_acc: str
    start: typing.Optional[int]
    stop: typing.Optional[int]
    outer_start: typing.Optional[int]
    outer_stop: typing.Optional[int]
    inner_start: typing.Optional[int]
    inner_stop: typing.Optional[int]
    ref: typing.Optional[str]
    alt: typing.Optional[str]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        return SequenceLocation(
            assembly=element.attrib.get("Assembly"),
            chrom=element.attrib.get("Assembly"),
            chrom_acc=element.attrib.get("Accession"),
            start=mapply(int, element.attrib.get("start")),
            stop=mapply(int, element.attrib.get("stop")),
            outer_start=mapply(int, element.attrib.get("outerStart")),
            outer_stop=mapply(int, element.attrib.get("ousterStop")),
            inner_start=mapply(int, element.attrib.get("innerStart")),
            inner_stop=mapply(int, element.attrib.get("innerStop")),
            ref=element.attrib.get("referenceAlleleVCF"),
            alt=element.attrib.get("alternateAlleleVCF"),
        )


@attr.s(frozen=True, auto_attribs=True)
class Measure:
    """Represent the relevant informatino from a Measure."""

    measure_type: str
    sequence_locations: typing.Dict[str, SequenceLocation]
    comments: typing.Tuple[str, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        return Measure(
            measure_type=element.attrib.get("Type"),
            sequence_locations={
                elem.attrib.get("Assembly"): SequenceLocation.from_element(elem)
                for elem in chain(
                    element.findall("./SequenceLocation"),
                    element.findall("./MeasureRelationship/SequenceLocation"),
                    )
            },
            comments=tuple(
                elem.text
                for elem in chain(
                    element.findall("./MeasureRelationship/Comment"), element.findall("./Comment")
                )
            ),
        )


@attr.s(frozen=True, auto_attribs=True)
class Trait:
    """Represent the relevant information from a Trait."""

    preferred_name: typing.Optional[str]
    alternate_names: typing.Tuple[str, ...]

    @classmethod
    def from_element(cls, element: ET.Element):
        elem = element.find("./Name/ElementValue[@Type='Preferred']")
        preferred_name = elem.text if elem is not None else None
        alternate_names = [
            elem.text for elem in element.findall("./Name/ElementValue[@Type='Alternate']")
        ]
        return Trait(preferred_name=preferred_name, alternate_names=tuple(alternate_names),)


@attr.s(frozen=True, auto_attribs=True)
class TraitSet:
    """Represent the relevant information from a TraitSet."""

    #: Value of the "Type" attribute.
    set_type: str
    #: Numeric identifier for the ClinVarSet.
    id_no: typing.Optional[int]
    #: The traits in the set.
    traits: typing.Tuple[Trait, ...]

    @classmethod
    def from_element(cls, element: ET.Element):
        traits = [Trait.from_element(elem) for elem in element.findall("./Trait")]
        return TraitSet(
            set_type=element.attrib.get("Type"),
            id_no=mapply(int, element.attrib.get("ID")),
            traits=tuple(traits),
        )


@attr.s(frozen=True, auto_attribs=True)
class MeasureSet:
    """Represent the relevant information from a MeasureSet."""

    set_type: str
    accession: str
    measures: typing.Tuple[Measure, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        measures = [Measure.from_element(elem) for elem in element.findall("./Measure")]
        return MeasureSet(
            set_type=element.attrib.get("Type"),
            accession=element.attrib.get("Acc"),
            measures=tuple(measures),
        )


@attr.s(frozen=True, auto_attribs=True)
class GenotypeSet:
    """Represents a genotype observation in ClinVar.

    NB: we introduce dummy sets even for non-compound variants.
    """

    set_type: str
    accession: str
    measure_sets: typing.Tuple[MeasureSet, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        if element.tag == "GenotypeSet":
            measure_sets = [
                MeasureSet.from_element(elem) for elem in element.findall("./MeasureSet")
            ]
            return GenotypeSet(
                set_type=element.attrib.get("Type"),
                accession=element.attrib.get("Acc"),
                measure_sets=tuple(measure_sets),
            )
        else:
            assert element.tag == "MeasureSet"
            return GenotypeSet(
                set_type=element.attrib.get("Type"),
                accession=element.attrib.get("Acc"),
                measure_sets=(MeasureSet.from_element(element),),
            )


@attr.s(frozen=True, auto_attribs=True)
class ReferenceClinVarAssertion:
    """Represent the relevant parts of a ReferenceClinVarAssertion."""

    #: Numeric identifier for the ClinVarSet.
    id_no: id
    #: Record status
    record_status: str
    #: Date of creation.
    date_created: datetime.date
    #: Date of last update.
    date_updated: datetime.date

    #: The accession number.
    clinvar_accession: str
    #: The version of the record.
    version_no: int
    #: Description where the variant was observed.
    observed_in: typing.Optional[ObservedIn]
    #: Genotype sets
    genotype_sets: typing.Tuple[GenotypeSet, ...]
    #: Trait sets.
    trait_sets: typing.Tuple[TraitSet, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        gts_elem = element.findall("./GenotypeSet")
        if gts_elem:
            genotype_sets = [GenotypeSet.from_element(elem) for elem in gts_elem]
        else:
            genotype_sets = [
                GenotypeSet.from_element(elem) for elem in element.findall("./MeasureSet")
            ]
        trait_sets = [TraitSet.from_element(elem) for elem in element.findall("./TraitSet")]
        return ReferenceClinVarAssertion(
            id_no=int(element.attrib.get("ID")),
            record_status=element.find("RecordStatus").text,
            date_created=isoparse(element.attrib.get("DateCreated")),
            date_updated=isoparse(element.attrib.get("DateLastUpdated")),
            clinvar_accession=element.find("ClinVarAccession").attrib.get("Acc"),
            version_no=int(element.find("ClinVarAccession").attrib.get("Version")),
            observed_in=ObservedIn.from_element(element.find("ObservedIn")),
            genotype_sets=tuple(genotype_sets),
            trait_sets=tuple(trait_sets),
        )


@attr.s(frozen=True, auto_attribs=True)
class ClinVarAssertion:
    """Represent the relevant parts of a ClinVarAssertion."""

    #: Numeric identifier for the ClinVarSet.
    id_no: id
    #: Record status
    record_status: str
    #: Date of submission.
    submitter_date: typing.Optional[datetime.date]

    #: The accession number.
    clinvar_accession: str
    #: The version of the record.
    version_no: int
    #: Description where the variant was observed.
    observed_in: typing.Optional[ObservedIn]
    #: Genotype sets
    genotype_sets: typing.Tuple[GenotypeSet, ...]
    #: Trait sets.
    trait_sets: typing.Tuple[TraitSet, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @classmethod
    def from_element(cls, element: ET.Element):
        submitter_date = element.find("ClinVarSubmissionID").attrib.get("submitterDate")
        gts_elem = element.findall("./GenotypeSet")
        if gts_elem:
            genotype_sets = [GenotypeSet.from_element(elem) for elem in gts_elem]
        else:
            genotype_sets = [
                GenotypeSet.from_element(elem) for elem in element.findall("./MeasureSet")
            ]
        trait_sets = [TraitSet.from_element(elem) for elem in element.findall("./TraitSet")]
        return ClinVarAssertion(
            id_no=int(element.attrib.get("ID")),
            record_status=element.find("RecordStatus").text,
            submitter_date=isoparse(submitter_date) if submitter_date else None,
            clinvar_accession=element.find("ClinVarAccession").attrib.get("Acc"),
            version_no=int(element.find("ClinVarAccession").attrib.get("Version")),
            observed_in=ObservedIn.from_element(element.find("ObservedIn")),
            genotype_sets=tuple(genotype_sets),
            trait_sets=tuple(trait_sets),
        )


@attr.s(frozen=True, auto_attribs=True)
class ClinVarSet:
    """Represent the relevant parts of a ClinVarSet."""

    #: Numeric identifier for the ClinVarSet.
    id_no: id
    #: Record status
    record_status: str
    #: Record title
    title: str
    #: The ReferenceClinVarAssertion, if any.
    ref_cv_assertion: typing.Optional[ReferenceClinVarAssertion]
    #: The ClinVarAssertion objects, if any.
    cv_assertions: typing.Tuple[ClinVarAssertion, ...]

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @staticmethod
    def from_element(element: ET.Element) -> TClinVarSet:
        return ClinVarSet(
            id_no=int(element.attrib.get("ID")),
            record_status=element.find("RecordStatus").text,
            title=element.find("Title").text,
            ref_cv_assertion=ReferenceClinVarAssertion.from_element(
                element.find("ReferenceClinVarAssertion")
            ),
            cv_assertions=tuple(
                ClinVarAssertion.from_element(element)
                for element in element.findall("ClinVarAssertion")
            ),
        )


@attr.s(frozen=True, auto_attribs=True)
class ReleaseSet:
    """Root tag representation."""

    #: The release date.
    release_date: datetime.date

    #: Auto-generate UUID for the object.
    uuid: UUID = attr.ib(factory=uuid4)

    @staticmethod
    def from_element(element: ET.Element, uuid=None):
        return ReleaseSet(release_date=isoparse(element.attrib.get("Dated")), uuid=uuid or uuid4())


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

    def run(self):
        logger.info("Parsing elements...")
        release_uuid = uuid4()
        for event, elem in ET.iterparse(self.input):
            if elem.tag == "ReleaseSet" and event == "end":
                release_set = ReleaseSet.from_element(elem, uuid=release_uuid)
                print(json.dumps(cattr.unstructure(release_set), cls=UUIDEncoder, indent="  "))
                print("--")
            elif elem.tag == "ClinVarSet" and event == "end":
                clinvar_set = ClinVarSet.from_element(elem)
                print(json.dumps(cattr.unstructure(clinvar_set), cls=UUIDEncoder, indent="  "))
                print("--")
                elem.clear()
                self.written_rows += 1
            if self.max_rows and self.written_rows >= self.max_rows:
                logger.info("Breaking out after writing %d rows (as configured)", self.written_rows)
                break
        logger.info("Done parsing elements")
