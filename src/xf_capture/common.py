#!/usr/bin/env python3
"""
Common utilities and shared functions for XfCapture.
"""

import shutil


def check_conda() -> str | None:
    """
    Check if conda or mamba is available.

    Returns:
        str: Command name ('mamba' or 'conda') if found, None otherwise
    """
    for cmd in ["mamba", "conda"]:
        if shutil.which(cmd):
            return cmd
    return None


def check_snakemake() -> bool:
    """
    Check if snakemake is available.

    Returns:
        bool: True if snakemake is found, False otherwise
    """
    return shutil.which("snakemake") is not None


def print_version() -> None:
    """Print XfCapture version information."""
    from importlib.metadata import version, PackageNotFoundError

    try:
        pkg_version = version("xf_capture")
        print(f"XfCapture version {pkg_version}")
    except PackageNotFoundError:
        print("XfCapture version unknown (development mode)")


def format_bytes(bytes_size: int) -> str:
    """
    Format bytes into human-readable string.

    Args:
        bytes_size: Size in bytes

    Returns:
        str: Formatted size (e.g., "1.5 GB")
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024.0
    return f"{bytes_size:.1f} PB"
