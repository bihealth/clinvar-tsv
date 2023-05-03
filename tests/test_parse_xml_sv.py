import contextlib

import pytest  # noqa

from clinvar_tsv.parse_clinvar_xml import ClinvarParser


def test_parse_92148661(tmpdir):
    """Test with record seen as problematic before"""
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/clinvar-92148661.xml", "rt"))
        out37_small = stack.push((tmpdir / "out37.small.tsv").open("wt"))
        out37_sv = stack.push((tmpdir / "out37.sv.tsv").open("wt"))
        out38_small = stack.push((tmpdir / "out38.small.tsv").open("wt"))
        out38_sv = stack.push((tmpdir / "out38.sv.tsv").open("wt"))
        parser = ClinvarParser(inputf, out37_small, out37_sv, out38_small, out38_sv)
        parser.run()

    with (tmpdir / "out37.sv.tsv").open("rt") as input37:
        lines = input37.readlines()
        assert len(lines) == 2
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh37",
            "1",
            "24171566",
            "24194858",
            "769",
            ".",
            ".",
            "Deletion",
            '{"FUCA1"}',
            '{"HGNC:4006"}',
            "VCV000000682",
            "RCV000000717",
        ]

    with (tmpdir / "out38.sv.tsv").open("rt") as input38:
        lines = input38.readlines()
        assert len(lines) == 2
        assert lines[0].split("\t")[:2] == ["release", "chromosome"]
        assert lines[1].split("\t")[:12] == [
            "GRCh38",
            "1",
            "23845077",
            "23868290",
            "95",
            ".",
            ".",
            "Deletion",
            '{"FUCA1"}',
            '{"HGNC:4006"}',
            "VCV000000682",
            "RCV000000717",
        ]
