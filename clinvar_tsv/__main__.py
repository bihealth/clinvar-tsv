"""Command line interface for clinvar-tsv"""

import argparse
import logging
import os.path
import sys

from .download import FtpDownloader


#: Host name of NCBI FTP.
NCBI_FTP_HOST = "ftp.ncbi.nlm.nih.gov"

#: Path to Clinvar XML file
PATH_CLINVAR_XML = "/pub/clinvar/xml"
#: Name of Clinvar XML file
NAME_CLINVAR_XML = "ClinVarFullRelease_00-latest.xml.gz"

#: Path to Clinvar summary TSV file
PATH_VARIANT_SUMMARY_TXT = "/pub/clinvar/tab_delimited"
#: Name of Clinvar summary TSV file
NAME_VARIANT_SUMMARY_TXT = "variant_summary.txt.gz"


def run_download(args):
    """Run the download process."""
    downloader = FtpDownloader(NCBI_FTP_HOST)
    downloader.download_if_changed(
        os.path.join(PATH_CLINVAR_XML, NAME_CLINVAR_XML),
        os.path.join(args.download_dir, NAME_CLINVAR_XML),
    )
    downloader.download_if_changed(
        os.path.join(PATH_VARIANT_SUMMARY_TXT, NAME_VARIANT_SUMMARY_TXT),
        os.path.join(args.download_dir, NAME_VARIANT_SUMMARY_TXT),
    )


def run(args):
    """Entry point after parsing command line arguments"""
    logging.basicConfig(level=logging.DEBUG)
    return args.func(args)


def main(argv=None):
    """Main entry point before parsing command line arguments."""
    parser = argparse.ArgumentParser("clinvar-tsv")
    parser.add_argument(
        "--download-dir", required=True, help="Directory to use for downloading files into."
    )
    # parser.add_argument(
    #     "--b37-path", required=True, help="Path to GRCh37 FAI-indexed FASTA file."
    # )
    # parser.add_argument(
    #     "--b38-path", required=True, help="Path to GRCh38 FAI-indexed FASTA file."
    # )

    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = "command"

    parser_download = subparsers.add_parser(
        "download", help="Download files from NCBI into --download-dir"
    )
    parser_download.set_defaults(func=run_download)

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
