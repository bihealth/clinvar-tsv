"""Merging of normalized ClinVar TSV files."""

import datetime
import enum
import json
import re
import typing

import attr
import cattr
from dateutil.parser import isoparse

from .parse_clinvar_xml import ClinVarSet, DateTimeEncoder, VariationClinVarRecord, as_pg_list

HEADER_OUT = (
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
    "point_rating",
    "pathogenicity",
    "review_status",
    "pathogenicity_summary",
    "details",
)


cattr.register_structure_hook(
    datetime.datetime, lambda dt_str, _: isoparse(dt_str) if dt_str else None
)
cattr.register_structure_hook(
    datetime.date, lambda dt_str, _: isoparse(dt_str).date() if dt_str else None
)
cattr.register_structure_hook(
    datetime.time, lambda dt_str, _: isoparse(dt_str).time() if dt_str else None
)

_ASSERTION_LABELS = {
    -2: ("benign",),
    -1: (
        "likely benign",
        # below: other
        "protective",
    ),
    0: (
        "uncertain significance",
        # below: other
        "association not found",
        "drug response",
        "not provided",
        "other",
    ),
    1: (
        "likely pathogenic",
        # below: other
        "affects",
        "association",
        "conflicting interpretations of pathogenicity",
        "risk factor",
        "confers sensitivity",
    ),
    2: ("pathogenic",),
}


class ClinVarAssertion(enum.Enum):
    BENIGN = -2
    LIKELY_BENIGN = -1
    UNCERTAIN = 0
    LIKELY_PATHOGENIC = 1
    PATHOGENIC = 2

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented

    @classmethod
    def from_label(cls, label):
        for val in cls:
            if label in _ASSERTION_LABELS[val.value]:
                return val
        raise ValueError("Invalid label: %s" % label)

    def label(self):
        return _ASSERTION_LABELS[self.value][0]


def most_extreme(assertions: typing.Iterable[ClinVarAssertion]):
    if min(assertions).value < 0 and max(assertions).value > 0:
        return max(assertions)
    elif min(assertions).value < 0:
        return min(assertions)
    else:
        return max(assertions)


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


class ClinVarReviewStatus(enum.Enum):
    CONFLICTING_INTERPRETATIONS = 1
    CRITERIA_PROVIDED = 2
    MULTIPLE_SUBMITTERS = 3
    NO_ASSERTION_CRITERIA_PROVIDED = 4
    NO_ASSERTION_PROVIDED = 5
    NO_CONFLICTS = 6
    PRACTICE_GUIDELINE = 7
    EXPERT_PANEL = 8
    SINGLE_SUBMITTER = 9

    @classmethod
    def from_label(cls, label):
        for val in cls:
            if label in _REVIEW_STATUS_LABELS[val.value]:
                return val
        raise ValueError("Invalid label: %s" % label)

    def label(self):
        return _REVIEW_STATUS_LABELS[self.value][0]


@attr.s(frozen=True, auto_attribs=True)
class ReviewedAssertion:
    review_statuses: typing.Tuple[ClinVarReviewStatus, ...]
    assertions: typing.Tuple[ClinVarAssertion, ...]

    def review_status_label(self):
        return ", ".join([s.label() for s in self.review_statuses])

    def assertions_max_label(self):
        return most_extreme(self.assertions).label()

    def assertions_summary_label(self):
        if len(self.assertions) == 1:
            return "; ".join([a.label() for a in self.assertions])
        else:
            tmp = {}
            assertions = tuple(sorted(self.assertions, reverse=True))
            for assertion in assertions:
                tmp.setdefault(assertion, 0)
                tmp[assertion] += 1
            return "; ".join(["%s(%d)" % (a.label(), tmp[a]) for a in assertions])

    def from_multiple(self):
        return (
            ClinVarReviewStatus.MULTIPLE_SUBMITTERS in self.review_statuses
            or self.num_assertions() > 1
        )

    def has_conflicting(self):
        return ClinVarReviewStatus.CONFLICTING_INTERPRETATIONS in self.review_statuses or (
            (self.num_assertions() > 1)
            and (self.min_assertion().value < 0)
            and (self.max_assertion().value > 0)
        )

    def min_assertion(self):
        return min(self.assertions)

    def max_assertion(self):
        return min(self.assertions)

    def num_assertions(self):
        return len(self.assertions)

    @classmethod
    def from_clinvar_set(cls, elem: ClinVarSet):
        return ReviewedAssertion(
            review_statuses=tuple(
                map(
                    ClinVarReviewStatus.from_label,
                    re.split(r", ?|/", elem.ref_cv_assertion.review_status),
                )
            ),
            assertions=tuple(
                map(
                    ClinVarAssertion.from_label,
                    re.split(r", ?|/", elem.ref_cv_assertion.pathogenicity),
                )
            ),
        )

    @classmethod
    def combine(cls, elems):
        review_statuses = (
            ClinVarReviewStatus.MULTIPLE_SUBMITTERS,
            ClinVarReviewStatus.NO_CONFLICTS,
        )
        assertions = []
        for elem in elems:
            assertions += list(elem.assertions)
        return ReviewedAssertion(
            review_statuses=review_statuses, assertions=tuple(sorted(set(assertions)))
        )


def is_germline(elem: ClinVarSet):
    try:
        return elem.ref_cv_assertion.observed_in.origin != "somatic"
    except AttributeError:
        return True  # assume germline if origin not given


def assertion_provided(elem: ClinVarSet):
    return "no assertion provided" not in elem.ref_cv_assertion.pathogenicity


def three_points(
    reviewed_assertions: typing.List[ReviewedAssertion],
) -> typing.Optional[ReviewedAssertion]:
    """Extract three point ReviewedAssertion of possible."""
    for review_status in (ClinVarReviewStatus.PRACTICE_GUIDELINE, ClinVarReviewStatus.EXPERT_PANEL):
        for reviewed_assertion in reviewed_assertions:
            if review_status in reviewed_assertion.review_statuses:
                return ReviewedAssertion(
                    review_statuses=(review_status,),
                    assertions=(most_extreme(reviewed_assertion.assertions),),
                )
    return None


def two_points(
    reviewed_assertions: typing.List[ReviewedAssertion],
) -> typing.Optional[ReviewedAssertion]:
    """Extract two point ReviewedAssertion of possible."""
    if not any(map(ReviewedAssertion.from_multiple, reviewed_assertions)):
        return None  # is single => no two points
    elif any(map(ReviewedAssertion.has_conflicting, reviewed_assertions)):
        return None
    else:
        return ReviewedAssertion.combine(reviewed_assertions)


def one_point(
    reviewed_assertions: typing.List[ReviewedAssertion],
) -> typing.Optional[ReviewedAssertion]:
    if (
        len(reviewed_assertions) == 1
        and ClinVarReviewStatus.MULTIPLE_SUBMITTERS not in reviewed_assertions[0].review_statuses
    ):  # single submitter
        return reviewed_assertions[0]
    else:  # multiple submitters, must have conflicts
        combined = ReviewedAssertion.combine(reviewed_assertions)
        assert reviewed_assertions[0].has_conflicting
        return combined


def zero_points(chunk) -> typing.Optional[ReviewedAssertion]:
    reviewed_assertions = list(map(ReviewedAssertion.from_clinvar_set, chunk))
    combined = ReviewedAssertion.combine(reviewed_assertions)
    if any(map(assertion_provided, filter(is_germline, chunk))):
        return None  # at least one germline with assertion => 1+ stars
    elif any(filter(is_germline, chunk)):  # all germline: no assertion
        return ReviewedAssertion(
            review_statuses=(ClinVarReviewStatus.NO_ASSERTION_PROVIDED,),
            assertions=combined.assertions,
        )
    else:  # all somatic
        return combined


def merge_and_write(valss, chunk, out_tsv):
    vals = valss[0]
    # Try to assign zero point assertion (all somatic or all germline have
    # no pathogenicity assigned.  If we have 1+ points then try three star
    # assignments with only germline variants and those with an assigned assertion.
    summary = zero_points(chunk)
    if summary:
        points = "0"
    else:
        points = "."  # counted as not in clinvar
    # Attempt to obtain three point assertions.
    if not summary:
        # Remove variants of somatic origin.
        chunk = list(filter(assertion_provided, filter(is_germline, chunk)))
        # Get assertion/pathogenicity pairs.
        reviewed_assertions = list(map(ReviewedAssertion.from_clinvar_set, chunk))
        summary = three_points(reviewed_assertions)
        if summary:
            points = "3"
    # If this fails, go on and try to find one point assertions.
    if not summary:
        summary = two_points(reviewed_assertions)
        if summary:
            points = "2"
    if not summary:
        summary = one_point(reviewed_assertions)
        if summary:
            points = "1"
    # Merge symbols & HGNC IDs.
    symbols, hgnc_ids = [], []
    for one_vals in valss:
        if one_vals["symbols"] != "{}":
            symbols += list(map(json.loads, one_vals["symbols"][1:-1].split(",")))
        if one_vals["hgnc_ids"] != "{}":
            hgnc_ids += list(map(json.loads, one_vals["hgnc_ids"][1:-1].split(",")))
    # Write out record.
    print(
        "\t".join(
            [
                vals["release"],
                vals["chromosome"],
                vals["start"],
                vals["end"],
                vals["bin"],
                vals["reference"],
                vals["alternative"],
                vals["variation_type"],
                as_pg_list(sorted(set(symbols))),
                as_pg_list(sorted(set(hgnc_ids))),
                vals["vcv"],
                points,
                summary.assertions_max_label(),
                summary.review_status_label(),
                summary.assertions_summary_label(),
                json.dumps([cattr.unstructure(entry) for entry in chunk], cls=DateTimeEncoder)
                .replace(r"\"", "'")
                .replace('"', '"""'),
            ]
        ),
        file=out_tsv,
    )


def merge_tsvs(in_tsv, out_tsv):
    header_in = in_tsv.readline().strip().split("\t")
    print("\t".join(HEADER_OUT), file=out_tsv)

    prev_vals = None
    chunk = []
    valss = []
    while True:
        line = in_tsv.readline().strip().replace('"""', '"')
        if not line:
            break
        vals = dict(zip(header_in, line.split("\t")))
        if prev_vals and vals["vcv"] != prev_vals["vcv"]:  # write chunk and start new one
            merge_and_write(valss, chunk, out_tsv)
            chunk = []
            valss = []
        prev_vals = vals
        obj = json.loads(vals["details"])
        chunk.append(cattr.structure(obj, ClinVarSet))
        valss.append(vals)
    if prev_vals:  # write final chunk
        merge_and_write(valss, chunk, out_tsv)
