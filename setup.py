#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import find_packages, setup

import versioneer

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


with open("README.md") as readme_file:
    readme = readme_file.read()

with open("CHANGELOG.md") as history_file:
    history = history_file.read()

test_requirements = parse_requirements("requirements/test.txt")
install_requirements = parse_requirements("requirements/base.txt")

setup(
    name="clinvar-tsv",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=("Python 3 library for accessing and managing BioMedical sheets"),
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    author="Manuel Holtgrewe",
    author_email="manuel.holtgrewe@bihealth.de",
    url="https://github.com/bihealth/clinvar-tsv",
    packages=find_packages(),
    package_dir={"clinvar_tsv": "clinvar_tsv"},
    entry_points={"console_scripts": ["clinvar_tsv = clinvar_tsv.__main__:main"]},
    include_package_data=True,
    install_requires=install_requirements,
    license="MIT license",
    zip_safe=False,
    keywords="clinvar",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    test_suite="tests",
    tests_require=test_requirements,
)
