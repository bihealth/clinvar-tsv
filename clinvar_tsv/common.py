"""Common code for clinvar-tsv"""

import datetime
import enum
from itertools import chain
import json
import typing
import xml.etree.ElementTree as ET

import attr
import cattr
from dateutil.parser import isoparse
from logzero import logger

cattr.register_structure_hook(
    datetime.datetime, lambda dt_str, _: isoparse(dt_str) if dt_str else None
)
cattr.register_structure_hook(
    datetime.date, lambda dt_str, _: isoparse(dt_str).date() if dt_str else None
)
cattr.register_structure_hook(
    datetime.time, lambda dt_str, _: isoparse(dt_str).time() if dt_str else None
)


def as_pg_list(vals: typing.Iterable[str]) -> str:
    """Convert to Postgres TSV list of strings."""
    return "{%s}" % (",".join(map(json.dumps, vals)))


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


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (datetime.datetime, datetime.date, datetime.time)):
            return obj.isoformat()
        return super().default(obj)  # pragma: no cover


_ClinicalSignificance = typing.TypeVar("ClinicalSignificance")
_ReferenceClinVarAssertion = typing.TypeVar("ReferenceClinVarAssertion")
_ClinVarAssertion = typing.TypeVar("ClinVarAssertion")
_ClinVarSet = typing.TypeVar("ClinVarSet")
_ClinvarAssertion = typing.TypeVar("ClinvarAssertion")


@attr.s(frozen=True, auto_attribs=True)
class ClinicalSignificance:
    """Represent clinical significance."""

    #: Date of last evaluation.
    date_evaluated: datetime.date
    #: Significance review status.
    review_status: str
    #: Significance description.
    description: typing.Optional[str]
    #: Comments.
    comments: typing.Tuple[str, ...]

    @staticmethod
    def from_element(element: ET.Element) -> _ClinicalSignificance:
        description = element.find("Description")
        return ClinicalSignificance(
            date_evaluated=_mapply(isoparse, element.attrib.get("DateLastEvaluated")),
            review_status=element.find("ReviewStatus").text,
            description=None if description is None else description.text,
            comments=tuple(elem.text for elem in element.findall("./Comment")),
        )


@attr.s(frozen=True, auto_attribs=True)
class ObservedDataDescription:
    """Relevant information from ObservedData/Attribute[@Type='Description']."""

    #: Optional description text.
    description: typing.Optional[str]
    #: PubMed IDs.
    pubmed_ids: typing.Tuple[int, ...]
    #: OMIM IDs.
    omim_ids: typing.Tuple[int, ...]

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


def _mapply(f, x):
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
            start = _mapply(int, element.attrib.get("start"))
            stop = _mapply(int, element.attrib.get("stop"))
            ref = element.attrib.get("referenceAlleleVCF")
            alt = element.attrib.get("referenceAlleleVCF")
        return SequenceLocation(
            assembly=element.attrib.get("Assembly"),
            chrom=element.attrib.get("Chr"),
            chrom_acc=element.attrib.get("Accession"),
            start=start,
            stop=stop,
            outer_start=_mapply(int, element.attrib.get("outerStart")),
            outer_stop=_mapply(int, element.attrib.get("ousterStop")),
            inner_start=_mapply(int, element.attrib.get("innerStart")),
            inner_stop=_mapply(int, element.attrib.get("innerStop")),
            ref=ref,
            alt=alt,
        )


@attr.s(frozen=True, auto_attribs=True)
class Measure:
    """Represent the relevant informatino from a Measure."""

    measure_type: str
    symbols: typing.Tuple[str, ...]
    hgnc_ids: typing.Tuple[str, ...]
    sequence_locations: typing.Dict[str, SequenceLocation]
    comments: typing.Tuple[str, ...]

    @classmethod
    def from_element(cls, element: ET.Element):
        symbols = [
            elem.text for elem in element.findall('.//Symbol/ElementValue[@Type="Preferred"]')
        ]
        hgnc_ids = [elem.attrib.get("ID") for elem in element.findall('.//XRef[@DB="HGNC"]')]
        return Measure(
            measure_type=element.attrib.get("Type"),
            symbols=tuple(sorted(set(symbols))),
            hgnc_ids=tuple(sorted(set(hgnc_ids))),
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
    #: Numeric id_no for the ClinVarSet.
    id_no: typing.Optional[int]
    #: The traits in the set.
    traits: typing.Tuple[Trait, ...]

    @classmethod
    def from_element(cls, element: ET.Element):
        traits = [Trait.from_element(elem) for elem in element.findall("./Trait")]
        return TraitSet(
            set_type=element.attrib.get("Type"),
            id_no=_mapply(int, element.attrib.get("ID")),
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

    #: Numeric id_no for the ClinVarSet.
    id_no: int
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
    clin_sigs: typing.Tuple[ClinicalSignificance, ...]

    #: Number of gold stars shown on ClinVar.
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
            if clin_sig.description is not None:
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

    #: Numeric id_no for the ClinVarSet.
    id_no: int
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
    #: Clinical significance entries.
    clin_sigs: typing.Tuple[ClinicalSignificance, ...]

    #: Review status
    review_status: str
    #: Assertion of pathogenicity.
    pathogenicity: str

    @classmethod
    def from_element(cls, element: ET.Element) -> _ClinvarAssertion:
        submitter_date = element.find("ClinVarSubmissionID").attrib.get("submitterDate")
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
        for clin_sig in clin_sigs:
            if clin_sig.description is not None:
                review_status = clin_sig.review_status
                pathogenicity = clin_sig.description.lower()

        return ClinVarAssertion(
            id_no=int(element.attrib.get("ID")),
            record_status=element.find("RecordStatus").text,
            submitter_date=isoparse(submitter_date) if submitter_date else None,
            clinvar_accession=element.find("ClinVarAccession").attrib.get("Acc"),
            version_no=int(element.find("ClinVarAccession").attrib.get("Version")),
            observed_in=ObservedIn.from_element(element.find("ObservedIn")),
            genotype_sets=tuple(genotype_sets),
            trait_sets=tuple(trait_sets),
            clin_sigs=tuple(clin_sigs),
            review_status=review_status,
            pathogenicity=pathogenicity,
        )


@attr.s(frozen=True, auto_attribs=True)
class ClinVarSet:
    """Represent the relevant parts of a ClinVarSet."""

    #: Numeric id_no for the ClinVarSet.
    id_no: int
    #: Record status
    record_status: str
    #: Record title
    title: str
    #: The ReferenceClinVarAssertion, if any.
    ref_cv_assertion: typing.Optional[ReferenceClinVarAssertion]
    #: The ClinVarAssertion objects, if any.
    cv_assertions: typing.Tuple[ClinVarAssertion, ...]

    @staticmethod
    def from_element(element: ET.Element) -> _ClinVarSet:
        et_title = element.find("Title")
        return ClinVarSet(
            id_no=int(element.attrib.get("ID")),
            record_status=element.find("RecordStatus").text,
            title=et_title.text if et_title else "[NO TITLE]",
            ref_cv_assertion=ReferenceClinVarAssertion.from_element(
                element.find("ReferenceClinVarAssertion")
            ),
            cv_assertions=tuple(
                ClinVarAssertion.from_element(element)
                for element in element.findall("ClinVarAssertion")
            ),
        )


@attr.s(frozen=True, auto_attribs=True)
class VariationClinVarRecord:
    """Aggregated information of multiple ``ClinVarSet`` records."""

    #: Genome build/release name.
    release: str
    #: Chromosome/contig name.
    chromosome: str
    #: Start position (1-based).
    start: int
    #: End position (1-based).
    end: int
    #: UCSC bin.
    bin: int
    #: Reference style (VCF file).
    reference: str
    #: Alternative allele (VCF style).
    alternative: str
    #: Either "snv" or "indel".
    variation_type: str
    #: VCV ClinVar accession.
    clinvar_accession: str


@attr.s(frozen=True, auto_attribs=True)
class ReleaseSet:
    """Root tag representation."""

    #: The release date.
    release_date: datetime.date

    @staticmethod
    def from_element(element: ET.Element):
        return ReleaseSet(release_date=isoparse(element.attrib.get("Dated")))


_PATHOGENICITY_LABELS = {
    -2: (
        "benign",
        # below: other
        "no known pathogenicity",
        "non-pathogenic",
        "poly",
    ),
    -1: (
        "likely benign",
        # below: other
        "probable-non-pathogenic",
        "probably not pathogenic",
        "protective",
        "suspected benign",
    ),
    0: (
        "uncertain significance",
        # below: other
        "affects",
        "association",
        "association not found",
        "cancer",
        "confers sensitivity",
        "drug response",
        "drug response",
        "drug-response",
        "histocompatibility",
        "not provided",
        "other",
        "protective",
        "risk factor",
        "uncertain",
        "unknown",
        "untested",
        "variant of unknown significance",
        "associated with leiomyomas",
    ),
    1: (
        "likely pathogenic",
        # below: other
        "affects",
        "association",
        "confers sensitivity",
        "conflicting interpretations of pathogenicity",
        "probable-pathogenic",
        "probably pathogenic",
        "risk factor",
        "suspected pathogenic",
    ),
    2: (
        "pathogenic",
        # below: other
        "moderate",
        "mut",
        "pathologic",
    ),
}


class Pathogenicity(enum.Enum):
    BENIGN = -2
    LIKELY_BENIGN = -1
    UNCERTAIN = 0
    LIKELY_PATHOGENIC = 1
    PATHOGENIC = 2

    def __lt__(self, other):
        if self.__class__ == other.__class__:
            return self.value < other.value
        else:
            raise TypeError(
                f"Incompatible types {self.__class__} vs. {other.__class__}"
            )  # pragma: no cover

    @classmethod
    def from_label(cls, label, variant_id):
        for val in cls:
            if label in _PATHOGENICITY_LABELS[val.value]:
                return val
        # We sometimes see "likely pathogenic - $something"
        for val in cls:
            for pl in _PATHOGENICITY_LABELS[val.value]:
                if label.startswith(pl):
                    return val
        logger.warning("Invalid label %s for variant %s", label, variant_id)  # pragma: no cover
        return cls.UNCERTAIN

    def label(self):
        return _PATHOGENICITY_LABELS[self.value][0]


_REVIEW_STATUS_LABELS = {
    1: ("conflicting interpretations", "conflicting interpretations of pathogenicity"),
    2: ("criteria provided",),
    3: ("multiple submitters",),
    4: ("no assertion criteria provided",),
    5: ("no assertion provided",),
    6: ("no conflicts",),
    7: ("practice guideline",),
    8: ("reviewed by expert panel",),
    9: ("single submitter",),
}


class ReviewStatus(enum.Enum):
    CONFLICTING_INTERPRETATIONS = 1
    CRITERIA_PROVIDED = 2
    MULTIPLE_SUBMITTERS = 3
    NO_ASSERTION_CRITERIA_PROVIDED = 4
    NO_ASSERTION_PROVIDED = 5
    NO_CONFLICTS = 6
    PRACTICE_GUIDELINE = 7
    EXPERT_PANEL = 8
    SINGLE_SUBMITTER = 9

    def __lt__(self, other):
        if self.__class__ == other.__class__:
            return self.value < other.value
        else:
            raise TypeError(
                f"Incompatible types {self.__class__} vs. {other.__class__}"
            )  # pragma: no cover

    @classmethod
    def from_label(cls, label, variant_id):
        for val in cls:
            if label in _REVIEW_STATUS_LABELS[val.value]:
                return val
        raise ValueError(f"Invalid label: {label} for {variant_id}")  # pragma: no cover

    def label(self):
        return _REVIEW_STATUS_LABELS[self.value][0]
