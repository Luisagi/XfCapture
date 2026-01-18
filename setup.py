#!/usr/bin/env python3
"""
Setup script for XfCapture
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="xf_capture",
    version="0.0.2",
    author="Luis F. Arias-Giraldo",
    description="Xylella fastidiosa Capture Sequencing Analysis Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
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
        "snakemake>=7.0",
    ],
    
    # Entry points for command-line scripts
    entry_points={
        "console_scripts": [
            "xf_capture=xf_capture:main",
        ],
    },
    
    # Classifiers
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    
    # Keywords
    keywords="xylella fastidiosa genomics snakemake pipeline bioinformatics",
)
