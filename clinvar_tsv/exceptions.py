"""Exceptions for ``clinvar-tsv``."""


class ClinvarTsvException(Exception):
    """Exception raised by ``clinvar-tsv``"""


class XmlParseException(Exception):
    """Raised on problems parsing the XML"""


class NoSequenceLocation(Exception):
    """Raised on problems parsing the XML"""
