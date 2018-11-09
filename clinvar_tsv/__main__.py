"""Command line interface for clinvar-tsv"""

import argparse
import gzip
import logging
import os
import os.path
import sys

import snakemake

from . import parse_clinvar_xml


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
        "config": {"b37_path": args.b37_path, "b38_path": args.b38_path},
        "printshellcmds": True,
        "verbose": True,
    }
    return snakemake.snakemake(**kwargs)


def open_maybe_gzip(path, mode):
    try:
        return gzip.open(path, mode)
    except ValueError:
        return open(path, mode)


def run_parse_xml(args):
    """Parse XML file."""
    build = {"b37": "GRCh37", "b38": "GRCh38"}
    with open_maybe_gzip(args.output_single, "wt") as single_f:
        with open_maybe_gzip(args.output_multi, "wt") as multi_f:
            parser = parse_clinvar_xml.ClinvarParser(
                open_maybe_gzip(args.clinvar_xml, "rb"),
                single_f,
                multi_f,
                genome_build=build[args.genome_build],
            )
            parser.run()


def run_normalize_tsv(args):
    touch(args.output_tsv)


def run(args):
    """Entry point after parsing command line arguments"""
    logging.basicConfig(level=logging.INFO)
    return args.func(args)


def main(argv=None):
    """Main entry point before parsing command line arguments."""
    parser = argparse.ArgumentParser("clinvar-tsv")

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
    parser_main.set_defaults(func=run_main)

    # -----------------------------------------------------------------------
    # Command: parse_xml
    # -----------------------------------------------------------------------

    parser_parse_xml = subparsers.add_parser("parse_xml", help="Parse the Clinvar XML")
    parser_parse_xml.add_argument("--clinvar-xml", required=True, help="Path to Clinvar XML file.")
    parser_parse_xml.add_argument(
        "--genome-build",
        required=True,
        help="The genome build this variant is for.",
        choices=("b37", "b38"),
    )
    parser_parse_xml.add_argument(
        "--output-single", required=True, help="Output path for single TSV file."
    )
    parser_parse_xml.add_argument(
        "--output-multi", required=True, help="Output path to multi TSV file."
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

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
