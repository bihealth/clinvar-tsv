"""Merging of normalized ClinVar TSV files."""

import itertools
import json
import re
import typing

import attr
import cattr

from clinvar_tsv.common import ClinVarSet, DateTimeEncoder, Pathogenicity, ReviewStatus, as_pg_list

HEADER_OUT = (
    "release",
    "chromosome",
    "start",
    "end",
    "bin",
    "reference",
    "alternative",
    "clinvar_version",
    "variation_type",
    "symbols",
    "hgnc_ids",
    "vcv",
    "summary_clinvar_review_status_label",
    "summary_clinvar_pathogenicity_label",
    "summary_clinvar_pathogenicity",
    "summary_clinvar_gold_stars",
    "summary_paranoid_review_status_label",
    "summary_paranoid_pathogenicity_label",
    "summary_paranoid_pathogenicity",
    "summary_paranoid_gold_stars",
    "details",
)


_ReviewedPathogenicity = typing.TypeVar("_ReviewedPathogenicity")


@attr.s(frozen=True, auto_attribs=True)
class ReviewedPathogenicity:
    review_statuses: typing.Tuple[ReviewStatus, ...]
    pathogenicities: typing.Tuple[Pathogenicity, ...]

    def review_status_label(self):
        return ", ".join([s.label() for s in self.review_statuses])

    def pathogenicity_label(self) -> str:
        if ReviewStatus.CONFLICTING_INTERPRETATIONS in self.review_statuses:
            result = "conflicting interpretations of pathogenicity - "
            counts = {}
            for pathogenicity in sorted(self.pathogenicities):
                counts.setdefault(pathogenicity, 0)
                counts[pathogenicity] += 1
            result += "; ".join("%s (%d)" % (k.label(), v) for k, v in counts.items())
            return result
        else:
            pathogenicities = list(sorted(set(self.pathogenicities)))
            return " / ".join(a.label() for a in pathogenicities)

    def pathogenicity_list(self, *, all_on_conflicts: bool) -> typing.List[str]:
        if (
            ReviewStatus.CONFLICTING_INTERPRETATIONS in self.review_statuses
            and not all_on_conflicts
        ):
            return [Pathogenicity.UNCERTAIN.label()]
        else:
            return [pathogenicity.label() for pathogenicity in sorted(set(self.pathogenicities))]

    def gold_stars(self):
        if ReviewStatus.PRACTICE_GUIDELINE in self.review_statuses:
            return 4
        elif ReviewStatus.EXPERT_PANEL in self.review_statuses:
            return 3
        elif (
            ReviewStatus.CRITERIA_PROVIDED in self.review_statuses
            and ReviewStatus.CRITERIA_PROVIDED in self.review_statuses
            and ReviewStatus.MULTIPLE_SUBMITTERS in self.review_statuses
            and ReviewStatus.NO_CONFLICTS in self.review_statuses
        ):
            return 2
        elif (
            ReviewStatus.CRITERIA_PROVIDED in self.review_statuses
            and ReviewStatus.CONFLICTING_INTERPRETATIONS in self.review_statuses
        ):
            return 1
        elif (
            ReviewStatus.CRITERIA_PROVIDED in self.review_statuses
            and ReviewStatus.SINGLE_SUBMITTER in self.review_statuses
        ):
            return 1
        else:
            return 0

    @classmethod
    def from_clinvar_set(cls, elem: ClinVarSet) -> typing.Tuple[_ReviewedPathogenicity, ...]:
        result = []
        for cv_assertion in elem.cv_assertions:
            result.append(
                ReviewedPathogenicity(
                    review_statuses=tuple(
                        map(
                            lambda label: ReviewStatus.from_label(label, elem.id_no),
                            re.split(r", ?|/", cv_assertion.review_status),
                        )
                    ),
                    pathogenicities=tuple(
                        map(
                            lambda label: Pathogenicity.from_label(label, elem.id_no),
                            re.split(r", ?|/", cv_assertion.pathogenicity),
                        )
                    ),
                )
            )
        return tuple(result)

    @classmethod
    def combine(cls, elems: typing.Iterable[_ReviewedPathogenicity]) -> _ReviewedPathogenicity:
        def sign(x):
            if x < 0:
                return -1
            elif x > 0:
                return 1
            else:
                return 0

        all_pathogenicities = []
        for elem in elems:
            all_pathogenicities += elem.pathogenicities
        all_pathogenicities = list(sorted(all_pathogenicities))

        all_statuses = []
        for elem in elems:
            all_statuses += elem.review_statuses
        all_statuses = list(sorted(all_statuses))

        multiple_submitters = len(elems) > 1

        if all_statuses == [ReviewStatus.PRACTICE_GUIDELINE]:
            review_statuses = [ReviewStatus.PRACTICE_GUIDELINE]
        elif all_statuses == [ReviewStatus.EXPERT_PANEL]:
            review_statuses = [ReviewStatus.EXPERT_PANEL]
        else:
            if multiple_submitters:
                review_statuses = [ReviewStatus.MULTIPLE_SUBMITTERS]
            else:
                review_statuses = [ReviewStatus.SINGLE_SUBMITTER]
            if ReviewStatus.CRITERIA_PROVIDED in all_statuses:
                review_statuses.append(ReviewStatus.CRITERIA_PROVIDED)
            else:
                review_statuses.append(ReviewStatus.NO_ASSERTION_CRITERIA_PROVIDED)

        signs = {sign(assertion.value) for assertion in set(all_pathogenicities)}
        if multiple_submitters:
            if len(signs) == 1:
                review_statuses.append(ReviewStatus.NO_CONFLICTS)
            else:
                review_statuses.append(ReviewStatus.CONFLICTING_INTERPRETATIONS)

        if 1 in signs and 0 in signs:
            assertions = [assertion for assertion in all_pathogenicities if assertion.value >= 0]
        elif 1 in signs and -1 in signs:
            assertions = [assertion for assertion in all_pathogenicities if assertion.value != 0]
        elif 0 in signs and -1 in signs:
            assertions = [assertion for assertion in all_pathogenicities if assertion.value <= 0]
        else:
            assertions = []
            for assertion in all_pathogenicities:
                if assertion not in assertions:
                    assertions.append(assertion)

        return ReviewedPathogenicity(
            review_statuses=tuple(review_statuses), pathogenicities=tuple(sorted(assertions))
        )


def rp_stratum(rp: ReviewedPathogenicity) -> int:
    """Get stratum from ClinVarAssertion record."""
    if ReviewStatus.PRACTICE_GUIDELINE in rp.review_statuses:
        return 4
    elif ReviewStatus.EXPERT_PANEL in rp.review_statuses:
        return 3
    elif ReviewStatus.CRITERIA_PROVIDED in rp.review_statuses:
        return 2
    else:
        return 1


def summarize(
    chunk: typing.List[ClinVarSet], *, stratify_by_review_status: bool
) -> ReviewedPathogenicity:
    """Create summary ``ReviewedPathogenicity`` from ``chunk``.

    If ``stratify_by_review_status`` then only the highest stratum will be used (e.g., assessments
    with assertion criteria beat those without).  Otherwise, all will be considered.
    """
    # Obtain list of ReviewedAssertion objects
    rps = list(itertools.chain(*map(ReviewedPathogenicity.from_clinvar_set, chunk)))
    # Process, possibly stratified by review status
    stratified = {}
    for rp in rps:
        stratum = rp_stratum(rp) if stratify_by_review_status else 0
        stratified.setdefault(stratum, [])
        stratified[stratum].append(rp)
    highest_stratum = stratified[max(stratified.keys())]

    return ReviewedPathogenicity.combine(highest_stratum)


def merge_and_write(clinvar_version, valss, chunk, out_tsv):
    # Summarize chunks in clinvar and paranoid way.
    clinvar_summary = summarize(chunk, stratify_by_review_status=True)
    paranoid_summary = summarize(chunk, stratify_by_review_status=False)

    # Concatenate symbols & HGNC IDs.
    vals = valss[0]
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
                clinvar_version,
                vals["variation_type"],
                as_pg_list(sorted(set(symbols))),
                as_pg_list(sorted(set(hgnc_ids))),
                vals["vcv"],
                clinvar_summary.review_status_label(),
                clinvar_summary.pathogenicity_label(),
                as_pg_list(clinvar_summary.pathogenicity_list(all_on_conflicts=False)),
                str(clinvar_summary.gold_stars()),
                paranoid_summary.review_status_label(),
                paranoid_summary.pathogenicity_label(),
                as_pg_list(paranoid_summary.pathogenicity_list(all_on_conflicts=True)),
                str(paranoid_summary.gold_stars()),
                json.dumps([cattr.unstructure(entry) for entry in chunk], cls=DateTimeEncoder)
                .replace(r"\"", "'")
                .replace('"', '"""'),
            ]
        ),
        file=out_tsv,
    )


def merge_tsvs(clinvar_version, in_tsv, out_tsv):
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
        # if vals["vcv"] not in ("VCV000210112", "VCV000243036"):
        #     continue
        if prev_vals and vals["vcv"] != prev_vals["vcv"]:  # write chunk and start new one
            merge_and_write(clinvar_version, valss, chunk, out_tsv)
            chunk = []
            valss = []
        prev_vals = vals
        obj = json.loads(vals["details"])
        chunk.append(cattr.structure(obj, ClinVarSet))
        valss.append(vals)
    if prev_vals:  # write final chunk
        merge_and_write(clinvar_version, valss, chunk, out_tsv)
