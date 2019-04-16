"""Tests for ``parse_clinvar_xml`` program."""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

import os.path

import pytest

from clinvar_tsv.parse_clinvar_xml import HEADER, ClinvarParser

# --------------------------------------------------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------------------------------------------------


@pytest.fixture
def input_clinvar_small_del():
    """Return handle to opened ``clinvar_small_del.xml``"""
    with open(
        os.path.join(os.path.dirname(__file__), "data/clinvar_small_del.xml"), "rt"
    ) as inputf:
        yield inputf


@pytest.fixture
def input_clinvar_snv():
    """Return handle to opened ``input_clinvar_snv.xml``"""
    with open(os.path.join(os.path.dirname(__file__), "data/clinvar_snv.xml"), "rt") as inputf:
        yield inputf


@pytest.fixture
def input_clinvar_two_measure_sets():
    """Return handle to opened ``input_clinvar_two_measure_sets.xml``"""
    with open(
        os.path.join(os.path.dirname(__file__), "data/clinvar_two_measure_sets.xml"), "rt"
    ) as inputf:
        yield inputf


@pytest.fixture
def output_single(tmpdir):
    """Return handle to temp file for single output"""
    with open(str(tmpdir.join("output_single.tsv")), "w+t") as outputf:
        yield outputf


@pytest.fixture
def output_multi(tmpdir):
    """Return handle to temp file for multi output"""
    with open(str(tmpdir.join("output_multi.tsv")), "w+t") as outputf:
        yield outputf


# --------------------------------------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------------------------------------


def assert_header_written(output_single, output_multi):
    """Assert that the output has bene written to both output single and multi."""
    # Flush output files and seek to start
    for handle in (output_single, output_multi):
        handle.flush()
        handle.seek(0)
    # Read in header of output lines
    header = output_single.readline().strip().split("\t")
    assert header == HEADER
    if output_single is not output_multi and output_single.name != output_multi.name:
        header = output_multi.readline().strip().split("\t")
        assert header == HEADER


def load_tsv(handle):
    lines = handle.readlines()
    result = []
    for line in lines:
        while line[-1] in "\r\n":
            line = line[:-1]
        result.append(line.split("\t"))
    return result


# --------------------------------------------------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------------------------------------------------


def test_parse_clinvar_snv_grch37(input_clinvar_snv, output_single, output_multi):
    """Test case of simple SNV for GRCh37"""
    # Parse input file
    parser = ClinvarParser(input_clinvar_snv, output_single, output_multi, "GRCh37")
    parser.run()
    # Check output file
    assert_header_written(output_single, output_multi)
    # Check subsequent
    lines_single = load_tsv(output_single)
    assert len(lines_single) == 1
    assert lines_single[0][:6] == ["GRCh37", "13", "95228658", "T", "C", "95228658"]
    assert lines_single[0][-1:] == ["0"]
    # Check clinical significance
    assert lines_single[0][19:24] == ['0', '1', '0', '0', '0']
    lines_multi = load_tsv(output_multi)
    assert lines_multi == []


def test_parse_clinvar_snv_grch38(input_clinvar_snv, output_single, output_multi):
    """Test case of simple SNV for GRCh38"""
    # Parse input file
    parser = ClinvarParser(input_clinvar_snv, output_single, output_multi, "GRCh38")
    parser.run()
    # Check output file
    assert_header_written(output_single, output_multi)
    # Check subsequent
    lines_single = load_tsv(output_single)
    assert len(lines_single) == 1
    assert lines_single[0][:6] == ["GRCh38", "chr13", "94576404", "T", "C", "94576404"]
    assert lines_single[0][-1:] == ["0"]
    # Check clinical significance
    assert lines_single[0][19:24] == ['0', '1', '0', '0', '0']
    lines_multi = load_tsv(output_multi)
    assert lines_multi == []


def test_parse_clinvar_two_measure_sets_grch37(
    input_clinvar_two_measure_sets, output_single, output_multi
):
    """Test case of two measure sets for GRCh37.

    This was a case that the original clinvar xml-to-tsv script would ignore.
    """
    # Parse input file
    parser = ClinvarParser(input_clinvar_two_measure_sets, output_single, output_multi, "GRCh37")
    parser.run()
    # Check output file
    assert_header_written(output_single, output_multi)
    # Check subsequent
    lines_single = load_tsv(output_single)
    assert lines_single == []
    lines_multi = load_tsv(output_multi)
    assert len(lines_multi) == 2
    assert lines_multi[0][:6] == ["GRCh37", "7", "103276772", "C", "T", "103276772"]
    assert lines_multi[0][-1:] == ["1"]
    assert lines_multi[1][:6] == ["GRCh37", "7", "103132416", "A", "C", "103132416"]
    assert lines_multi[1][-1:] == ["1"]


def test_parse_clinvar_two_measure_sets_grch38(
    input_clinvar_two_measure_sets, output_single, output_multi
):
    """Test case of two measure sets for GRCh38.

    This was a case that the original clinvar xml-to-tsv script would ignore.
    """
    # Parse input file
    parser = ClinvarParser(input_clinvar_two_measure_sets, output_single, output_multi, "GRCh38")
    parser.run()
    # Check output file
    assert_header_written(output_single, output_multi)
    # Check subsequent
    lines_single = load_tsv(output_single)
    assert lines_single == []
    lines_multi = load_tsv(output_multi)
    assert len(lines_multi) == 2
    assert lines_multi[0][:6] == ["GRCh38", "chr7", "103636325", "C", "T", "103636325"]
    assert lines_multi[0][-1:] == ["1"]
    assert lines_multi[1][:6] == ["GRCh38", "chr7", "103491969", "A", "C", "103491969"]
    assert lines_multi[1][-1:] == ["1"]
