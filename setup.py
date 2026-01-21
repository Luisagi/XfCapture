#!/usr/bin/env python3
"""
Setup script for XfCapture
"""

from setuptools import setup, find_packages
from pathlib import Path

setup(
    name="xf_capture",
    version="0.0.2",
    author="Luis F. Arias-Giraldo",
    description=
    "This pipeline processes paired-end sequencing data from \
    X. fastidiosa targeted sequence capture enrichment (Xf-TSCE) ",
    url="https://github.com/luisagi/XfCapture",
    license="MIT",
    
    # Package configuration
    packages=find_packages(),
    py_modules=["xf_capture"],
    
    # Include non-Python files
    include_package_data=True,
    package_data={
        "": [
            "Snakefile",
            "config.yaml",
            "envs/*.yaml",
            "data/*.tar.gz",
            "utils/*.py",
            "utils/*.R",
        ],
    },
    
    # Dependencies
    python_requires=">=3.8",
    install_requires=[
        "pyyaml>=5.0",
        "snakemake>=9.0",
    ],
    
    # Entry points for command-line scripts
    entry_points={
        "console_scripts": [
            "xf_capture=xf_capture:main",
        ],
    }
)
