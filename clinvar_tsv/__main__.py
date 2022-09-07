"""Command line interface for clinvar-tsv"""

import argparse
import gzip
import logging
import os
import os.path
import sys

import snakemake

from clinvar_tsv import __version__

from . import merge_tsvs, normalize, parse_clinvar_xml


def run_inspect(args):
    kwargs = {
        "snakefile": os.path.join(os.path.dirname(__file__), "Snakefile"),
        "workdir": args.work_dir,
        "summary": True,
        "verbose": True,
    }
    snakemake.snakemake(**kwargs)


def run_main(args):
    """Entry point for the main program."""
    kwargs = {
        "snakefile": os.path.join(os.path.dirname(__file__), "Snakefile"),
        "workdir": args.work_dir,
        "config": {
            "b37_path": args.b37_path,
            "b38_path": args.b38_path,
            "debug": args.debug,
            "clinvar_version": args.clinvar_version,
        },
        "printshellcmds": True,
        "verbose": True,
        "force_incomplete": True,
        "cores": 8,
    }
    return not snakemake.snakemake(**kwargs)


def open_maybe_gzip(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)


def run_parse_xml(args):
    """Parse XML file."""
    with open(args.output_b37, "wt") as output_b37:
        with open(args.output_b38, "wt") as output_b38:
            parser = parse_clinvar_xml.ClinvarParser(
                input_file=open_maybe_gzip(args.clinvar_xml, "rb"),
                out_b37=output_b37,
                out_b38=output_b38,
                max_rcvs=args.max_rcvs,
            )
            parser.run()


def run_normalize_tsv(args):
    with open(args.input_tsv, "rt") as input_tsv:
        with open(args.output_tsv, "wt") as output_tsv:
            normalize.normalize_tab_delimited_file(input_tsv, output_tsv, args.reference)


def run_merge_tsvs(args):
    with open(args.input_tsv, "rt") as input_tsv:
        with open(args.output_tsv, "wt") as output_tsv:
            merge_tsvs.merge_tsvs(args.clinvar_version, input_tsv, output_tsv)


def run(args):
    """Entry point after parsing command line arguments"""
    logging.basicConfig(level=logging.INFO)
    return args.func(args)


def main(argv=None):
    """Main entry point before parsing command line arguments."""
    parser = argparse.ArgumentParser("clinvar-tsv")

    parser.add_argument(
        "--version", action="version", version="%(prog)s {version}".format(version=__version__)
    )

    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = "command"

    # -----------------------------------------------------------------------
    # Command: inspect
    # -----------------------------------------------------------------------

    parser_inspect = subparsers.add_parser("inspect", help="Show files to be created")
    parser_inspect.set_defaults(func=run_inspect)
    parser_inspect.add_argument("--work-dir", default=os.getcwd(), help="Path to working directory")

    # -----------------------------------------------------------------------
    # Command: main
    # -----------------------------------------------------------------------

    parser_main = subparsers.add_parser("main", help="Run the full process pipeline")
    parser_main.add_argument(
        "--b37-path", required=True, help="Path to GRCh37 FAI-indexed FASTA file."
    )
    parser_main.add_argument(
        "--b38-path", required=True, help="Path to GRCh38 FAI-indexed FASTA file."
    )
    parser_main.add_argument("--work-dir", default=os.getcwd(), help="Path to working directory")
    parser_main.add_argument("--cores", default=1, type=int, help="Number of cores to use")
    parser_main.add_argument(
        "--debug", default=False, action="store_true", help="Enables debugging helps"
    )
    parser_main.add_argument(
        "--clinvar-version", required=True, help="String to put as clinvar version"
    )
    parser_main.set_defaults(func=run_main)

    # -----------------------------------------------------------------------
    # Command: parse_xml
    # -----------------------------------------------------------------------

    parser_parse_xml = subparsers.add_parser("parse_xml", help="Parse the Clinvar XML")
    parser_parse_xml.add_argument("--clinvar-xml", required=True, help="Path to Clinvar XML file.")
    parser_parse_xml.add_argument(
        "--output-b37", required=True, help="Output path for GRCh37 file."
    )
    parser_parse_xml.add_argument(
        "--output-b38", required=True, help="Output path for GRCh38 file."
    )
    parser_parse_xml.add_argument(
        "--max-rcvs", required=False, type=int, help="Maximal number of RCV records to process."
    )
    parser_parse_xml.set_defaults(func=run_parse_xml)

    # -----------------------------------------------------------------------
    # Command: normalize_tsv
    # -----------------------------------------------------------------------

    parser_normalize_tsv = subparsers.add_parser("normalize_tsv", help="Parse the Clinvar XML")
    parser_normalize_tsv.add_argument(
        "--reference", required=True, help="Path to reference FASTA file"
    )
    parser_normalize_tsv.add_argument("--input-tsv", required=True, help="Path to input TSV file.")
    parser_normalize_tsv.add_argument(
        "--output-tsv", required=True, help="Path to output TSV file."
    )
    parser_normalize_tsv.set_defaults(func=run_normalize_tsv)

    # -----------------------------------------------------------------------
    # Command: merge_tsvs
    # -----------------------------------------------------------------------

    parser_merge_tsvs = subparsers.add_parser(
        "merge_tsvs", help="Merge TSV file (result: one per VCV)"
    )
    parser_merge_tsvs.add_argument("--input-tsv", required=True, help="Path to input TSV file.")
    parser_merge_tsvs.add_argument("--output-tsv", required=True, help="Path to output TSV file.")
    parser_merge_tsvs.add_argument(
        "--clinvar-version", required=True, help="String to put as clinvar version"
    )
    parser_merge_tsvs.set_defaults(func=run_merge_tsvs)

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
