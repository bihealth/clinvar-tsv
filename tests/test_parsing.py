import contextlib

import pytest

from clinvar_tsv.parse_clinvar_xml import ClinvarParser


def test_parse_74722873(tmpdir):
    """Test with record seen as problematic before"""
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/clinvar-74722873.xml", "rt"))
        out37 = stack.push((tmpdir / "out37.tsv").open("wt"))
        out38 = stack.push((tmpdir / "out38.tsv").open("wt"))
        parser = ClinvarParser(inputf, out37, out38)
        res = parser.run()


def test_parse_74722873_in_context(tmpdir):
    """Test with record seen as problematic before, this time in contet"""
    with contextlib.ExitStack() as stack:
        inputf = stack.push(open("tests/data/clinvar-in-context-74722873.xml", "rt"))
        out37 = stack.push((tmpdir / "out37.tsv").open("wt"))
        out38 = stack.push((tmpdir / "out38.tsv").open("wt"))
        parser = ClinvarParser(inputf, out37, out38)
        parser.run()
