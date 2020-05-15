"""Parse Clinvar XML

Code taken from MacArthur lab.

Reference on clinvar XML tag:

- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/clinvar_submission.xsd
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README
"""

# TODO:
#  - extract ACMG rating
#  - extract how assertion level
#  - add star fields

import datetime
import json
import typing
import re
import xml.etree.ElementTree as ET
from dateutil.parser import isoparse
from itertools import chain

import attr
import cattr
import binning
from logzero import logger
import tqdm

#: Mapping from review status to gold stars.
GOLD_STAR_MAP = {
    "no assertion provided": 0,
    "no assertion criteria provided": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting interpretations": 1,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}


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
        "variation_id",
        "rcv",
        "gold_stars",
        "review_status",
        "pathogenicity",
        "details",
    )
)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (datetime.datetime, datetime.date, datetime.time)):
            return obj.isoformat()
        return super().default(obj)


def replace_semicolons(s, replace_with=":"):
    """Replace semicolons with colons."""
    return s.replace(";", replace_with)


def remove_newlines_and_tabs(s):
    """Remove newline and tab characters."""
    return re.sub("[\t\n\r]", " ", s)


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

    @staticmethod
    def from_element(element: ET.Element) -> TClinicalSignificance:
        return ClinicalSignificance(
            date_evaluated=mapply(isoparse, element.attrib.get("DateLastEvaluated")),
            review_status=element.find("ReviewStatus").text,
            description=element.find("Description").text,
            comments=tuple(elem.text for elem in element.findall("./Comment")),
        )


@attr.s(frozen=True, auto_attribs=True)
class ObservedDataDescription:
    """Relevant information from ObservedData/Attribute[@Type='Description']."""

    #: Optional description text.
    description: typing.Optional[str]
    #: PubMed IDs.
    pubmed_ids: typing.Tuple[int]
    #: OMIM IDs.
    omim_ids: typing.Tuple[int]

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
                description=description, pubmed_ids=tuple(pubmed_ids), omim_ids=tuple(omim_ids)
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

    @classmethod
    def from_element(cls, element: ET.Element):
        if "positionVCF" in element.attrib:
            ref = element.attrib["referenceAlleleVCF"]
            alt = element.attrib["alternateAlleleVCF"]
            start = int(element.attrib["positionVCF"])
            stop = start + len(ref) - 1
        else:
            start = mapply(int, element.attrib.get("start"))
            stop = mapply(int, element.attrib.get("stop"))
            ref = element.attrib.get("referenceAlleleVCF")
            alt = element.attrib.get("referenceAlleleVCF")
        return SequenceLocation(
            assembly=element.attrib.get("Assembly"),
            chrom=element.attrib.get("Chr"),
            chrom_acc=element.attrib.get("Accession"),
            start=start,
            stop=stop,
            outer_start=mapply(int, element.attrib.get("outerStart")),
            outer_stop=mapply(int, element.attrib.get("ousterStop")),
            inner_start=mapply(int, element.attrib.get("innerStart")),
            inner_stop=mapply(int, element.attrib.get("innerStop")),
            ref=ref,
            alt=alt,
        )


@attr.s(frozen=True, auto_attribs=True)
class Measure:
    """Represent the relevant informatino from a Measure."""

    measure_type: str
    sequence_locations: typing.Dict[str, SequenceLocation]
    comments: typing.Tuple[str, ...]

    @classmethod
    def from_element(cls, element: ET.Element):
        return Measure(
            measure_type=element.attrib.get("Type"),
            sequence_locations={
                elem.attrib.get("Assembly"): SequenceLocation.from_element(elem)
                for elem in element.findall("./SequenceLocation")
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
        return Trait(preferred_name=preferred_name, alternate_names=tuple(alternate_names))


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
    #: Clinical significance entries.
    clin_sigs: typing.Tuple[ClinicalSignificance]

    #: Number of gold stars show on ClinVar.
    gold_stars: int
    #: Review status
    review_status: str
    #: Assertion of pathogenicity.
    pathogenicity: str

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
        clin_sigs = [
            ClinicalSignificance.from_element(elem)
            for elem in element.findall("./ClinicalSignificance")
        ]

        review_status = "no assertion criteria provided"
        pathogenicity = "uncertain significance"
        gold_stars = 0
        for clin_sig in clin_sigs:
            review_status = clin_sig.review_status
            pathogenicity = clin_sig.description.lower()
            gold_stars = GOLD_STAR_MAP[review_status]

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
            clin_sigs=tuple(clin_sigs),
            gold_stars=gold_stars,
            review_status=review_status,
            pathogenicity=pathogenicity,
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

    @staticmethod
    def from_element(element: ET.Element):
        return ReleaseSet(release_date=isoparse(element.attrib.get("Dated")))


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
                                            variation_type = "snv"
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
                                            clinvar_set.ref_cv_assertion.id_no,
                                            clinvar_set.ref_cv_assertion.clinvar_accession,
                                            clinvar_set.ref_cv_assertion.gold_stars,
                                            clinvar_set.ref_cv_assertion.review_status,
                                            clinvar_set.ref_cv_assertion.pathogenicity,
                                            json.dumps(
                                                cattr.unstructure(clinvar_set), cls=DateTimeEncoder
                                            ).replace('"', '"""'),
                                        ]
                                        print("\t".join(map(str, row)), file=out_files[build])
                        progress.update()
                    elem.clear()
                if self.max_rcvs and self.rcvs >= self.max_rcvs:
                    logger.info("Breaking out after processing %d RCVs (as configured)", self.rcvs)
                    break
        logger.info("Done parsing elements")
