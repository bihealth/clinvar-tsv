import datetime

import factory
from pytest_factoryboy import register

from clinvar_tsv.common import (
    GOLD_STAR_MAP,
    ClinVarAssertion,
    ClinVarSet,
    ReferenceClinVarAssertion,
)


class RcvaFactory(factory.Factory):
    class Meta:
        model = ReferenceClinVarAssertion

    id_no = factory.Sequence(lambda n: 10_000 + n)
    record_status = "current"
    date_created = factory.LazyFunction(datetime.date.today)
    date_updated = factory.LazyFunction(datetime.date.today)

    clinvar_accession = factory.Sequence(lambda n: "RCV%06d" % (10_000 + n))
    version_no = 1
    observed_in = None
    genotype_sets = ()
    trait_sets = ()
    clin_sigs = ()

    gold_stars = factory.Iterator(GOLD_STAR_MAP.values())
    review_status = factory.Iterator(GOLD_STAR_MAP.keys())
    pathogenicity = factory.Iterator(
        ("pathogenic", "likely pathogenic", "uncertain significance", "likely benign", "benign")
    )


class CvaFactory(factory.Factory):
    class Meta:
        model = ClinVarAssertion

    id_no = factory.Sequence(lambda n: 10_000 + n)
    record_status = "current"
    submitter_date = factory.LazyFunction(datetime.date.today)

    clinvar_accession = factory.Sequence(lambda n: "SCV%06d" % (10_000 + n))
    version_no = 1
    observed_in = None
    genotype_sets = ()
    trait_sets = ()
    clin_sigs = ()

    review_status = factory.Iterator(GOLD_STAR_MAP.keys())
    pathogenicity = factory.Iterator(
        ("pathogenic", "likely pathogenic", "uncertain significance", "likely benign", "benign")
    )


class CvsFactory(factory.Factory):
    class Meta:
        model = ClinVarSet

    id_no = factory.Sequence(lambda n: 10_000 + n)
    record_status = "current"
    title = factory.LazyAttribute(lambda o: f"ClinVarSet %s{o.id_no}")
    ref_cv_assertion = factory.SubFactory(RcvaFactory)

    @factory.lazy_attribute
    def cv_assertions(self):
        return (CvaFactory(), CvaFactory())


register(CvaFactory)
register(RcvaFactory)
register(CvsFactory)
