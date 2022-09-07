import contextlib

import pytest  # noqa

from clinvar_tsv.parse_clinvar_xml import ClinvarParser


def test_parse_74722873(tmpdir):
    """Test with record seen as problematic before"""
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/clinvar-74722873.xml", "rt"))
        out37 = stack.push((tmpdir / "out37.tsv").open("wt"))
        out38 = stack.push((tmpdir / "out38.tsv").open("wt"))
        parser = ClinvarParser(inputf, out37, out38)
        parser.run()

    with (tmpdir / "out37.tsv").open("rt") as input37:
        lines = input37.readlines()
        assert len(lines) == 2
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh37",
            "9",
            "13150558",
            "13150558",
            "685",
            "A",
            "C",
            "snv",
            '{"MPDZ"}',
            '{"HGNC:7208"}',
            "VCV000376860",
            "RCV000430591",
        ]

    with (tmpdir / "out38.tsv").open("rt") as input38:
        lines = input38.readlines()
        assert len(lines) == 2
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh38",
            "9",
            "13150559",
            "13150559",
            "685",
            "A",
            "C",
            "snv",
            '{"MPDZ"}',
            '{"HGNC:7208"}',
            "VCV000376860",
            "RCV000430591",
        ]


def test_parse_74722873_in_context(tmpdir):
    """Test with record seen as problematic before, this time in contet"""
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/clinvar-in-context-74722873.xml", "rt"))
        out37 = stack.push((tmpdir / "out37.tsv").open("wt"))
        out38 = stack.push((tmpdir / "out38.tsv").open("wt"))
        parser = ClinvarParser(inputf, out37, out38)
        parser.run()

    with (tmpdir / "out37.tsv").open("rt") as input37:
        lines = input37.readlines()
        assert len(lines) == 71
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh37",
            "14",
            "105246551",
            "105246551",
            "1387",
            "C",
            "T",
            "snv",
            '{"AKT1"}',
            '{"HGNC:391"}',
            "VCV000013983",
            "RCV000430173",
        ]
        assert lines[2].split("\t")[:12] == [
            "GRCh37",
            "17",
            "41245624",
            "41245624",
            "899",
            "C",
            "G",
            "snv",
            '{"BRCA1"}',
            '{"HGNC:1100"}',
            "VCV000054401",
            "RCV000430192",
        ]

    with (tmpdir / "out38.tsv").open("rt") as input38:
        lines = input38.readlines()
        assert len(lines) == 71
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh38",
            "14",
            "104780214",
            "104780214",
            "1384",
            "C",
            "T",
            "snv",
            '{"AKT1"}',
            '{"HGNC:391"}',
            "VCV000013983",
            "RCV000430173",
        ]
        assert lines[2].split("\t")[:12] == [
            "GRCh38",
            "17",
            "43093607",
            "43093607",
            "913",
            "C",
            "G",
            "snv",
            '{"BRCA1"}',
            '{"HGNC:1100"}',
            "VCV000054401",
            "RCV000430192",
        ]


def test_parse_spta1(tmpdir):
    """Test with record seen as problematic before"""
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/clinvar-spta1.xml", "rt"))
        out37 = stack.push((tmpdir / "out37.tsv").open("wt"))
        out38 = stack.push((tmpdir / "out38.tsv").open("wt"))
        parser = ClinvarParser(inputf, out37, out38)
        parser.run()

    with (tmpdir / "out37.tsv").open("rt") as input37:
        lines = input37.readlines()
        assert len(lines) == 3
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh37",
            "1",
            "158624528",
            "158624528",
            "1795",
            "G",
            "T",
            "snv",
            '{"SPTA1"}',
            '{"HGNC:11272"}',
            "VCV000012846",
            "RCV000013699",
        ]
        assert lines[2].split("\t")[:12] == [
            "GRCh37",
            "1",
            "158624528",
            "158624528",
            "1795",
            "G",
            "T",
            "snv",
            '{"SPTA1"}',
            '{"HGNC:11272"}',
            "VCV000012846",
            "RCV000251633",
        ]

    with (tmpdir / "out38.tsv").open("rt") as input38:
        lines = input38.readlines()
        assert len(lines) == 3
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh38",
            "1",
            "158654738",
            "158654738",
            "1795",
            "G",
            "T",
            "snv",
            '{"SPTA1"}',
            '{"HGNC:11272"}',
            "VCV000012846",
            "RCV000013699",
        ]
        assert lines[2].split("\t")[:12] == [
            "GRCh38",
            "1",
            "158654738",
            "158654738",
            "1795",
            "G",
            "T",
            "snv",
            '{"SPTA1"}',
            '{"HGNC:11272"}',
            "VCV000012846",
            "RCV000251633",
        ]
