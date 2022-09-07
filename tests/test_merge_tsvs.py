import contextlib

import pytest  # noqa

from clinvar_tsv.merge_tsvs import (
    Pathogenicity,
    ReviewedPathogenicity,
    ReviewStatus,
    merge_tsvs,
    summarize,
)


def test_summarize_stratified_single(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            )
        ],
        stratify_by_review_status=True,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.SINGLE_SUBMITTER,
            ReviewStatus.CRITERIA_PROVIDED,
        ),
        pathogenicities=(Pathogenicity.LIKELY_BENIGN,),
    )
    assert result == expected

    assert result.review_status_label() == "single submitter, criteria provided"
    assert result.pathogenicity_label() == "likely benign"
    assert result.gold_stars() == 1


def test_summarize_unstratified_single(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            )
        ],
        stratify_by_review_status=False,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.SINGLE_SUBMITTER,
            ReviewStatus.CRITERIA_PROVIDED,
        ),
        pathogenicities=(Pathogenicity.LIKELY_BENIGN,),
    )
    assert result == expected

    assert result.review_status_label() == "single submitter, criteria provided"
    assert result.pathogenicity_label() == "likely benign"
    assert result.gold_stars() == 1


def test_summarize_stratified_lb_lp(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="no assertion criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="no assertion criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=0,
                    review_status="no assertion criteria provided, single submitter",
                    pathogenicity="likely pathogenic",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="no assertion criteria provided, single submitter",
                        pathogenicity="likely pathogenic",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=False,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.MULTIPLE_SUBMITTERS,
            ReviewStatus.NO_ASSERTION_CRITERIA_PROVIDED,
            ReviewStatus.CONFLICTING_INTERPRETATIONS,
        ),
        pathogenicities=(Pathogenicity.LIKELY_BENIGN, Pathogenicity.LIKELY_PATHOGENIC),
    )
    assert result == expected

    assert (
        result.review_status_label()
        == "multiple submitters, no assertion criteria provided, conflicting interpretations"
    )
    assert (
        result.pathogenicity_label()
        == "conflicting interpretations of pathogenicity - likely benign (1); likely pathogenic (1)"
    )
    assert result.gold_stars() == 0


def test_summarize_unstratified_lb_lp(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=0,
                    review_status="no assertion criteria provided, single submitter",
                    pathogenicity="likely pathogenic",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="no assertion criteria provided, single submitter",
                        pathogenicity="likely pathogenic",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=False,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.MULTIPLE_SUBMITTERS,
            ReviewStatus.CRITERIA_PROVIDED,
            ReviewStatus.CONFLICTING_INTERPRETATIONS,
        ),
        pathogenicities=(Pathogenicity.LIKELY_BENIGN, Pathogenicity.LIKELY_PATHOGENIC),
    )
    assert result == expected

    assert (
        result.review_status_label()
        == "multiple submitters, criteria provided, conflicting interpretations"
    )
    assert (
        result.pathogenicity_label()
        == "conflicting interpretations of pathogenicity - likely benign (1); likely pathogenic (1)"
    )
    assert result.gold_stars() == 1


def test_summarize_stratified_lb_vus(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=0,
                    review_status="no assertion criteria provided, single submitter",
                    pathogenicity="uncertain significance",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="no assertion criteria provided, single submitter",
                        pathogenicity="uncertain significance",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=True,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.SINGLE_SUBMITTER,
            ReviewStatus.CRITERIA_PROVIDED,
        ),
        pathogenicities=(Pathogenicity.LIKELY_BENIGN,),
    )
    assert result == expected

    assert result.review_status_label() == "single submitter, criteria provided"
    assert result.pathogenicity_label() == "likely benign"
    assert result.gold_stars() == 1


def test_summarize_unstratified_lb_vus(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=0,
                    review_status="no assertion criteria provided, single submitter",
                    pathogenicity="uncertain significance",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="no assertion criteria provided, single submitter",
                        pathogenicity="uncertain significance",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=False,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.MULTIPLE_SUBMITTERS,
            ReviewStatus.CRITERIA_PROVIDED,
            ReviewStatus.CONFLICTING_INTERPRETATIONS,
        ),
        pathogenicities=(Pathogenicity.LIKELY_BENIGN, Pathogenicity.UNCERTAIN),
    )
    assert result == expected

    assert (
        result.review_status_label()
        == "multiple submitters, criteria provided, conflicting interpretations"
    )
    assert (
        result.pathogenicity_label()
        == "conflicting interpretations of pathogenicity - likely benign (1); uncertain significance (1)"
    )
    assert result.gold_stars() == 1


def test_summarize_unstratified_lb_b(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely benign",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="no assertion criteria provided, single submitter",
                    pathogenicity="benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="no assertion criteria provided, single submitter",
                        pathogenicity="benign",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=False,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.MULTIPLE_SUBMITTERS,
            ReviewStatus.CRITERIA_PROVIDED,
            ReviewStatus.NO_CONFLICTS,
        ),
        pathogenicities=(
            Pathogenicity.BENIGN,
            Pathogenicity.LIKELY_BENIGN,
        ),
    )
    assert result == expected

    assert result.review_status_label() == "multiple submitters, criteria provided, no conflicts"
    assert result.pathogenicity_label() == "benign / likely benign"
    assert result.gold_stars() == 2


def test_summarize_unstratified_lp_p(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="likely pathogenic",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="likely pathogenic",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="pathogenic",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="pathogenic",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=True,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.MULTIPLE_SUBMITTERS,
            ReviewStatus.CRITERIA_PROVIDED,
            ReviewStatus.NO_CONFLICTS,
        ),
        pathogenicities=(
            Pathogenicity.LIKELY_PATHOGENIC,
            Pathogenicity.PATHOGENIC,
        ),
    )
    assert result == expected

    assert result.review_status_label() == "multiple submitters, criteria provided, no conflicts"
    assert result.pathogenicity_label() == "likely pathogenic / pathogenic"
    assert result.gold_stars() == 2


def test_summarize_unstratified_vus_vus(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="uncertain significance",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="uncertain significance",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="uncertain significance",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="uncertain significance",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=True,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(
            ReviewStatus.MULTIPLE_SUBMITTERS,
            ReviewStatus.CRITERIA_PROVIDED,
            ReviewStatus.NO_CONFLICTS,
        ),
        pathogenicities=(Pathogenicity.UNCERTAIN,),
    )
    assert result == expected

    assert result.review_status_label() == "multiple submitters, criteria provided, no conflicts"
    assert result.pathogenicity_label() == "uncertain significance"
    assert result.gold_stars() == 2


def test_summarize_expert_panel(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=3,
                    review_status="reviewed by expert panel",
                    pathogenicity="likely pathogenic",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="reviewed by expert panel",
                        pathogenicity="likely pathogenic",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=1,
                    review_status="criteria provided, single submitter",
                    pathogenicity="benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="criteria provided, single submitter",
                        pathogenicity="benign",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=True,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(ReviewStatus.EXPERT_PANEL,),
        pathogenicities=(Pathogenicity.LIKELY_PATHOGENIC,),
    )
    assert result == expected

    assert result.review_status_label() == "reviewed by expert panel"
    assert result.pathogenicity_label() == "likely pathogenic"
    assert result.gold_stars() == 3


def test_summarize_practice_guideline(cvs_factory, rcva_factory, cva_factory):
    result = summarize(
        [
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=4,
                    review_status="practice guideline",
                    pathogenicity="likely pathogenic",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="practice guideline",
                        pathogenicity="likely pathogenic",
                    ),
                ),
            ),
            cvs_factory(
                ref_cv_assertion=rcva_factory(
                    gold_stars=3,
                    review_status="reviewed by expert panel",
                    pathogenicity="benign",
                ),
                cv_assertions=(
                    cva_factory(
                        review_status="reviewed by expert panel",
                        pathogenicity="benign",
                    ),
                ),
            ),
        ],
        stratify_by_review_status=True,
    )
    expected = ReviewedPathogenicity(
        review_statuses=(ReviewStatus.PRACTICE_GUIDELINE,),
        pathogenicities=(Pathogenicity.LIKELY_PATHOGENIC,),
    )
    assert result == expected

    assert result.review_status_label() == "practice guideline"
    assert result.pathogenicity_label() == "likely pathogenic"
    assert result.gold_stars() == 4


def test_merge_tsvs_spta1(tmpdir):
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/parsed-74722873.37.tsv", "rt"))
        out37 = stack.push((tmpdir / "merged-74722873.37.tsv").open("wt"))
        merge_tsvs("VER", inputf, out37)

    with (tmpdir / "merged-74722873.37.tsv").open("rt") as inputf:
        lines = inputf.readlines()
        assert len(lines) == 2
        assert lines[0].split("\t") == [
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
            "details\n",
        ]
        assert lines[1].split("\t")[:20] == [
            "GRCh37",
            "9",
            "13150558",
            "13150558",
            "685",
            "A",
            "C",
            "VER",
            "snv",
            '{"MPDZ"}',
            '{"HGNC:7208"}',
            "VCV000376860",
            "multiple submitters, criteria provided, conflicting interpretations",
            (
                "conflicting interpretations of pathogenicity - benign (1); "
                "uncertain significance (1)"
            ),
            '{"uncertain significance"}',
            "1",
            "multiple submitters, criteria provided, conflicting interpretations",
            (
                "conflicting interpretations of pathogenicity - benign (1); uncertain "
                "significance (1)"
            ),
            '{"benign","uncertain significance"}',
            "1",
        ]
