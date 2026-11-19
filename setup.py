#!/usr/bin/env python3
"""
Setup script for adjusted-identity package.
"""

from setuptools import setup, find_packages
import os
import setuptools

# The license is declared here rather than in pyproject.toml so that one source
# tree builds under both old and new setuptools. Versions from 77 on take an
# SPDX expression and have deprecated the older pyproject.toml form
# (license = {file = "LICENSE"}); 75.3, the last version that runs on Python
# 3.8, rejects an SPDX string in pyproject.toml and only knows the plain
# "license" field. Once Python 3.8 support is dropped this can become
# license = "BSD-2-Clause" in pyproject.toml.
if int(setuptools.__version__.split(".")[0]) >= 77:
    license_args = {"license_expression": "BSD-2-Clause"}
else:
    license_args = {"license": "BSD-2-Clause"}

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

setup(
    name="adjusted-identity",
    version="0.2.11",
    author="Josh Walker",
    author_email="joshowalker@yahoo.com",
    description="Adjusted Identity Calculator for DNA Sequences with MycoBLAST-style preprocessing and MSA support",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/joshuaowalker/adjusted-identity",
    **license_args,
    packages=find_packages(exclude=["tests", "tests.*"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.8",
    install_requires=[
        "edlib>=1.3.9",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
        ],
    },
    keywords="bioinformatics dna sequence alignment mycology identity barcode MycoBLAST ITS",
    project_urls={
        "Bug Reports": "https://github.com/joshuaowalker/adjusted-identity/issues",
        "Source": "https://github.com/joshuaowalker/adjusted-identity",
    },
)