#!/usr/bin/env python3

import subprocess
import sys
import yaml
from pathlib import Path

from xf_capture.setup import get_default_workflow_dir

# Package paths
PACKAGE_DIR = Path(__file__).resolve().parent
SNAKEFILE = PACKAGE_DIR / "workflows" / "Snakefile"

# Phase 3 rules, in dependency order. Used by --force to rerun the whole
# multi-phylo subgraph without touching Phase 1/2 rules (fastp, kraken2, etc.).
MULTI_PHYLO_RULES = [
    "prepare_multi_alignment_files",
    "align_and_trim_multi_genes",
    "multi_phylogenetic_analysis",
    "plot_multi_tree",
    "multi_phylogeny_phase",
]


def generate_multiphylo_config(output_dir: str, samples_list: str, iqtree_threads: int = 8) -> Path:
    """
    Extend an existing run_config.yaml (from a prior `xf_capture run`) with
    the multiphylo section, without modifying the original file.

    Args:
        output_dir: Directory from a prior `xf_capture run` (must contain run_config.yaml)
        samples_list: Path to the --samples-list text file
        iqtree_threads: Number of threads for the joint IQ-TREE job (overrides the
            value stored by the original `xf_capture run` for this analysis only)

    Returns:
        Path to the generated multiphylo_config.yaml
    """
    output_path = Path(output_dir).resolve()
    run_config_path = output_path / "run_config.yaml"

    if not run_config_path.exists():
        print(f"[Error] No run_config.yaml found in: {output_path}")
        print("[Error] multi-phylo requires an output directory from a completed 'xf_capture run'.")
        sys.exit(1)

    samples_list_path = Path(samples_list).resolve()
    if not samples_list_path.exists():
        print(f"[Error] Samples list file does not exist: {samples_list_path}")
        sys.exit(1)

    with open(run_config_path, 'r') as f:
        config = yaml.safe_load(f)

    config["multiphylo"] = {
        "samples_list": str(samples_list_path),
    }
    config.setdefault("threads", {})["iqtree"] = iqtree_threads

    config_path = output_path / "multiphylo_config.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return config_path


def run_multiphylo(
    output_dir: str,
    samples_list: str,
    workflow_dir: str | None = None,
    cores: int = 16,
    iqtree_threads: int = 8,
    alignment_jobs: int = 4,
    force: bool = False,
    extra_args: list | None = None,
) -> None:
    """
    Run the joint multi-sample phylogenetic analysis (Phase 3) over an
    existing XfCapture results directory.

    Args:
        output_dir: Directory from a prior `xf_capture run` (Phase 1 must be complete)
        samples_list: Path to a text file with one sample name per line
        workflow_dir: Path to workflow directory (from setup command)
        cores: Total number of CPU cores available to Snakemake
        iqtree_threads: Number of threads for the joint IQ-TREE job
        alignment_jobs: Number of parallel MAFFT/ClipKit jobs
        force: Force Snakemake to rerun the whole multi-phylo subgraph
            (prepare/align/tree/plot/summary) even if outputs already exist
        extra_args: Additional arguments to pass to Snakemake
    """
    print("\n" + "="*70)
    print("[Multi-Phylo] Starting joint multi-sample phylogenetic analysis")
    print("="*70 + "\n")

    output_path = Path(output_dir).resolve()
    if not output_path.exists():
        print(f"[Error] Output directory does not exist: {output_path}")
        sys.exit(1)

    if workflow_dir is None:
        workflow_dir = get_default_workflow_dir()
        if workflow_dir:
            print(f"[Multi-Phylo] Using saved workflow directory: {workflow_dir}")

    print(f"[Multi-Phylo] Output directory: {output_path}")
    print(f"[Multi-Phylo] Samples list:     {Path(samples_list).resolve()}")

    config_path = generate_multiphylo_config(output_dir, samples_list, iqtree_threads=iqtree_threads)
    print(f"[Multi-Phylo] Config file:      {config_path}")

    conda_prefix = None
    if workflow_dir:
        conda_prefix = Path(workflow_dir) / "conda_envs"
        if conda_prefix.exists():
            print(f"[Multi-Phylo] Conda prefix:     {conda_prefix}")

    if not SNAKEFILE.exists():
        print(f"[Error] Snakefile not found: {SNAKEFILE}")
        print("[Error] Please ensure the package was installed correctly.")
        sys.exit(1)

    snakemake_cmd = [
        "snakemake",
        "multi_phylogeny_phase",
        "--snakefile", str(SNAKEFILE),
        "--configfile", str(config_path),
        "--cores", str(cores),
        "--use-conda",
    ]

    if conda_prefix and conda_prefix.exists():
        snakemake_cmd.extend(["--conda-prefix", str(conda_prefix)])

    # --resources and --forcerun both take nargs='+' and must not be followed
    # by anything else positional (e.g. a target rule name) or it gets
    # swallowed as a resource/forced rule.
    snakemake_cmd.extend([
        "--resources",
        f"alignment_jobs={alignment_jobs}",
        "iqtree_jobs=1",
    ])

    if force:
        print("[Multi-Phylo] --force: rerunning the full multi-phylo subgraph")
        snakemake_cmd.extend(["--forcerun"] + MULTI_PHYLO_RULES)

    if extra_args:
        snakemake_cmd.extend(extra_args)

    print(f"\n[Multi-Phylo] Executing Snakemake target 'multi_phylogeny_phase'")
    print("="*70)

    try:
        result = subprocess.run(snakemake_cmd)
        if result.returncode != 0:
            print("\n[Error] Multi-phylo analysis failed")
            sys.exit(1)
    except KeyboardInterrupt:
        print("\n[Info] Multi-phylo analysis interrupted by user")
        sys.exit(130)

    print("\n" + "="*70)
    print("[Multi-Phylo] Joint phylogenetic analysis completed successfully!")
    print("="*70)
