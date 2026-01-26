#!/usr/bin/env python3

import subprocess
import sys
import yaml
from pathlib import Path

from xf_capture.setup import get_default_workflow_dir
from xf_capture.common import check_snakemake, check_conda

# Package paths
PACKAGE_DIR = Path(__file__).resolve().parent
SNAKEFILE = PACKAGE_DIR / "workflows" / "Snakefile"
RESOURCES_DIR = PACKAGE_DIR / "resources"


def generate_config(
    input_dir: str,
    output_dir: str,
    workflow_dir: str | None = None,
    kraken_db: str | None = None,
    k2_mapping_memory: bool = False
) -> Path:
    """
    Generate a config file for the Snakemake run.

    When workflow_dir is provided, uses the reference sequences and databases
    from that directory. Otherwise, requires manual specification of paths.

    Args:
        input_dir: Directory containing input FASTQ files
        output_dir: Directory for pipeline output files
        workflow_dir: Path to workflow directory (created with setup)
        kraken_db: Path to Kraken2 database (overrides workflow_dir setting)
        k2_mapping_memory: Enable Kraken2 memory mapping mode

    Returns:
        Path to the generated config file
    """
    # Base configuration
    config = {
        "directories": {
            "fastq_dir": str(Path(input_dir).resolve()),
            "output_dir": str(Path(output_dir).resolve()),
        },
        "references": {
            "probes": "/path/to/probes.fasta",
            "xf_genomes": "/path/to/xf_genomes",
        },
        "kraken2": {
            "database": kraken_db if kraken_db else "/path/to/kraken2/database",
            "memory_mapping": k2_mapping_memory,
        }
    }

    # If workflow_dir is provided, load existing configuration and override paths
    if workflow_dir:
        workflow_config_path = Path(workflow_dir) / "xf_capture_config.yaml"
        if workflow_config_path.exists():
            with open(workflow_config_path, 'r') as f:
                wf_config = yaml.safe_load(f)

            if "xf_capture" in wf_config:
                # Override Kraken2 database path
                config["kraken2"]["database"] = wf_config["xf_capture"].get(
                    "kraken2_database", config["kraken2"]["database"]
                )
                # Override reference sequences paths from workflow_dir
                if "probes" in wf_config["xf_capture"]:
                    config["references"]["probes"] = wf_config["xf_capture"]["probes"]
                if "xf_genomes" in wf_config["xf_capture"]:
                    config["references"]["xf_genomes"] = wf_config["xf_capture"]["xf_genomes"]

    # Override with explicit kraken_db if provided
    if kraken_db:
        config["kraken2"]["database"] = kraken_db

    # Save config to output directory
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    config_path = output_path / "run_config.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return config_path


def run_pipeline(
    input_dir: str,
    output_dir: str,
    workflow_dir: str | None = None,
    kraken_db: str | None = None,
    cores: int = 16,
    kraken_jobs: int = 1,
    alignment_jobs: int = 4,
    iqtree_jobs: int = 2,
    auto: bool = True,
    k2_mapping_memory: bool = False,
    extra_args: list | None = None,
) -> None:
    """
    Run the XfCapture pipeline.

    This function coordinates the Snakemake workflow execution in two phases:
    1. Phase 1: QC, taxonomy, probes reconstruction, and MLST
    2. Phase 2: Phylogenetic analysis (only for successful samples)

    Args:
        input_dir: Directory containing paired-end FASTQ files
        output_dir: Directory where results will be written
        workflow_dir: Path to workflow directory (from setup command)
        kraken_db: Path to Kraken2 database (overrides workflow_dir)
        cores: Total number of CPU cores to use
        kraken_jobs: Number of parallel Kraken2 jobs
        alignment_jobs: Number of parallel alignment jobs
        iqtree_jobs: Number of parallel IQ-TREE jobs
        auto: Continue to Phase 2 automatically without prompting
        k2_mapping_memory: Enable Kraken2 memory mapping mode
        extra_args: Additional arguments to pass to Snakemake
    """
    print("\n" + "="*70)
    print("[Run] Starting XfCapture pipeline")
    print("="*70 + "\n")

    # Check prerequisites
    if not check_snakemake():
        print("[Error] Snakemake not found. Please install Snakemake first.")
        print("[Error] Installation: conda install -c bioconda snakemake")
        sys.exit(1)

    # If workflow_dir not provided, try to load from user configuration
    if workflow_dir is None:
        workflow_dir = get_default_workflow_dir()
        if workflow_dir:
            print(f"[Run] Using saved workflow directory: {workflow_dir}")
        else:
            print("[Warning] No workflow directory configured.")
            print("[Warning] Run 'xf_capture setup --dir <path>' first,")
            print("[Warning] or specify --workflow-dir manually.")

    # Validate input directory
    input_path = Path(input_dir).resolve()
    if not input_path.exists():
        print(f"[Error] Input directory does not exist: {input_path}")
        sys.exit(1)

    # Generate config file
    print(f"[Run] Input directory:  {input_path}")
    print(f"[Run] Output directory: {Path(output_dir).resolve()}")

    config_path = generate_config(
        input_dir=input_dir,
        output_dir=output_dir,
        workflow_dir=workflow_dir,
        kraken_db=kraken_db,
        k2_mapping_memory=k2_mapping_memory
    )
    print(f"[Run] Config file:      {config_path}")

    # Determine conda prefix
    conda_prefix = None
    if workflow_dir:
        conda_prefix = Path(workflow_dir) / "conda_envs"
        if conda_prefix.exists():
            print(f"[Run] Conda prefix:     {conda_prefix}")

    # Verify Snakefile exists
    if not SNAKEFILE.exists():
        print(f"[Error] Snakefile not found: {SNAKEFILE}")
        print("[Error] Please ensure the package was installed correctly.")
        sys.exit(1)

    # Build snakemake command
    snakemake_cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--configfile", str(config_path),
        "--cores", str(cores),
        "--use-conda",
        "--resources",
        f"kraken_jobs={kraken_jobs}",
        f"alignment_jobs={alignment_jobs}",
        f"iqtree_jobs={iqtree_jobs}"
    ]

    if conda_prefix and conda_prefix.exists():
        snakemake_cmd.extend(["--conda-prefix", str(conda_prefix)])

    if extra_args:
        snakemake_cmd.extend(extra_args)

    print(f"\n[Run] Executing Snakemake with {cores} cores")
    print("="*70)

    # -------------------------------------------------------------------------
    # Phase 1: Main pipeline up to MLST
    # -------------------------------------------------------------------------
    print("\n[Phase 1] Running pipeline up to MLST + Checkpoints...\n")

    try:
        result = subprocess.run(snakemake_cmd)
        if result.returncode != 0:
            print("\n[Error] Phase 1 failed")
            sys.exit(1)
    except KeyboardInterrupt:
        print("\n[Info] Pipeline interrupted by user")
        sys.exit(130)

    print("\n[Phase 1] Completed. Checking successful samples...")

    # -------------------------------------------------------------------------
    # Check for successful samples
    # -------------------------------------------------------------------------
    output_path = Path(output_dir).resolve()
    checkpoint_dir = output_path / "03.probes_reconstruction" / "checkpoint"

    if checkpoint_dir.exists():
        checkpoint_files = list(checkpoint_dir.glob("*_reconstruction_check.txt"))
        successful = 0
        for f in checkpoint_files:
            content = f.read_text()
            if "SUCCESS:" in content:
                successful += 1

        print(f"[Run] Total samples:      {len(checkpoint_files)}")
        print(f"[Run] Successful samples: {successful}")

        if successful == 0:
            print("\n[Info] No successful samples. Pipeline terminated.")
            return

        # Ask to continue with phylogeny (unless auto mode)
        if not auto:
            try:
                response = input("\n[Run] Continue with phylogenetic analysis? (y/N): ")
                if response.lower() not in ['y', 'yes']:
                    print("[Info] Pipeline terminated by user choice.")
                    return
            except EOFError:
                print("\n[Info] Non-interactive mode. Use --auto to continue automatically.")
                return

        # -------------------------------------------------------------------------
        # Phase 2: Phylogenetic analysis
        # -------------------------------------------------------------------------
        print("\n[Phase 2] Running phylogenetic analysis...\n")

        phylo_cmd = snakemake_cmd.copy()
        phylo_cmd.append("phylogeny_phase")

        try:
            result = subprocess.run(phylo_cmd)
            if result.returncode != 0:
                print("\n[Error] Phase 2 failed")
                sys.exit(1)
        except KeyboardInterrupt:
            print("\n[Info] Pipeline interrupted by user")
            sys.exit(130)

    print("\n" + "="*70)
    print("[Run] XfCapture pipeline completed successfully!")
    print("="*70)
