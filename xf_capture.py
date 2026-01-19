#!/usr/bin/env python3
"""
XfCapture - A Snakemake-based pipeline for Xylella fastidiosa genomic analysis

Usage:
    xf_capture setup --dir <directory>
    xf_capture run --input_dir <dir> --output_dir <dir> [options]

Examples:
    # Setup the workflow (install conda environments and download databases)
    xf_capture setup --dir /path/to/workflow_data

    # Run the pipeline
    xf_capture run --input_dir /path/to/fastq --output_dir /path/to/results
"""

import argparse
import os
import sys
import subprocess
import shutil
import yaml
from pathlib import Path


# Get the directory where xf_capture is installed
SCRIPT_DIR = Path(__file__).resolve().parent
SNAKEFILE = SCRIPT_DIR / "Snakefile"
CONFIG_TEMPLATE = SCRIPT_DIR / "config.yaml"
ENVS_DIR = SCRIPT_DIR / "envs"
REFERENCE_SEQS_DIR = SCRIPT_DIR / "reference_seqs"

# User configuration directory (persistent settings)
USER_CONFIG_DIR = Path.home() / ".xf_capture"
USER_CONFIG_FILE = USER_CONFIG_DIR / "config.yaml"


def print_header():
    """Print the XfCapture header."""
    header = """
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║   ██╗  ██╗███████╗ ██████╗ █████╗ ██████╗ ████████╗██╗   ██╗██████╗  ║
║   ╚██╗██╔╝██╔════╝██╔════╝██╔══██╗██╔══██╗╚══██╔══╝██║   ██║██╔══██╗ ║
║    ╚███╔╝ █████╗  ██║     ███████║██████╔╝   ██║   ██║   ██║██████╔╝ ║
║    ██╔██╗ ██╔══╝  ██║     ██╔══██║██╔═══╝    ██║   ██║   ██║██╔══██╗ ║
║   ██╔╝ ██╗██║     ╚██████╗██║  ██║██║        ██║   ╚██████╔╝██║  ██║ ║
║   ╚═╝  ╚═╝╚═╝      ╚═════╝╚═╝  ╚═╝╚═╝        ╚═╝    ╚═════╝ ╚═╝  ╚═╝ ║
║                                                                      ║
║   Xylella fastidiosa Capture Sequencing Analysis Pipeline            ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
    """
    print(header)


def check_conda():
    """Check if conda/mamba is available."""
    for cmd in ["mamba", "conda"]:
        if shutil.which(cmd):
            return cmd
    return None


def check_snakemake():
    """Check if snakemake is available."""
    return shutil.which("snakemake") is not None


def save_user_config(workflow_dir: str):
    """
    Save the workflow directory path to user configuration file.
    This allows subsequent runs to find the workflow without --workflow_dir.
    """
    USER_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    
    config_data = {
        "workflow_dir": str(Path(workflow_dir).resolve())
    }
    
    with open(USER_CONFIG_FILE, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)
    
    return USER_CONFIG_FILE


def get_default_workflow_dir() -> str:
    """
    Get the default workflow directory from user configuration.
    Returns None if not configured.
    """
    if not USER_CONFIG_FILE.exists():
        return None
    
    try:
        with open(USER_CONFIG_FILE, 'r') as f:
            config = yaml.safe_load(f)
        
        if config and "workflow_dir" in config:
            workflow_path = Path(config["workflow_dir"])
            if workflow_path.exists():
                return str(workflow_path)
            else:
                print(f"[Warning] Configured workflow_dir does not exist: {workflow_path}")
                return None
    except Exception as e:
        print(f"[Warning] Failed to read user config: {e}")
        return None
    
    return None


def setup_workflow(workflow_dir: str):
    """
    Setup the XfCapture workflow:
    - Create directory structure
    - Extract reference sequences
    - Install conda environments
    - Download required databases (Kraken2)
    
    Skips components that are already installed.
    """
    print_header()
    print("\n[Setup] Starting XfCapture workflow setup...\n")
    
    workflow_path = Path(workflow_dir).resolve()
    
    # Create directory structure
    dirs_to_create = [
        workflow_path / "conda_envs",
        workflow_path / "databases" / "kraken2",
        workflow_path / "reference_seqs",
    ]
    
    for dir_path in dirs_to_create:
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"[Setup] Created directory: {dir_path}")
    
    # -------------------------------------------------------------------------
    # Extract reference sequences from bundled tar.gz
    # -------------------------------------------------------------------------
    print("\n" + "="*70)
    print("[Setup] Reference sequences...")
    print("="*70)
    
    ref_seqs_tarball = SCRIPT_DIR / "data" / "ref_seqs.tar.gz"
    ref_seqs_dest = workflow_path / "reference_seqs"
    xf_genomes_dir = ref_seqs_dest / "xf_genomes"
    probes_file = ref_seqs_dest / "probes.fasta"
    
    # Check if reference sequences are already extracted
    if probes_file.exists() and xf_genomes_dir.exists() and any(xf_genomes_dir.glob("*.fna")):
        print(f"[Setup] ✓ Reference sequences already extracted. Skipping.")
    else:
        if not ref_seqs_tarball.exists():
            print(f"[Error] Reference sequences archive not found: {ref_seqs_tarball}")
            print("[Error] Please ensure the package was installed correctly.")
            sys.exit(1)
        
        print(f"[Setup] Source archive: {ref_seqs_tarball}")
        print(f"[Setup] Destination:    {ref_seqs_dest}")
        
        try:
            # Extract the main tar.gz
            import tarfile
            import gzip
            
            with tarfile.open(ref_seqs_tarball, "r:gz") as tar:
                tar.extractall(path=ref_seqs_dest)
            print("[Setup] ✓ Reference sequences extracted")
            
            # Decompress all .gz files in xf_genomes/
            if xf_genomes_dir.exists():
                gz_files = list(xf_genomes_dir.glob("*.gz"))
                if gz_files:
                    print(f"[Setup] Decompressing {len(gz_files)} genome files...")
                    
                    for gz_file in gz_files:
                        output_file = gz_file.with_suffix("")  # Remove .gz extension
                        with gzip.open(gz_file, 'rb') as f_in:
                            with open(output_file, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        gz_file.unlink()  # Remove the .gz file after extraction
                    
                    print(f"[Setup] ✓ {len(gz_files)} genome files decompressed")
            
        except Exception as e:
            print(f"[Error] Failed to extract reference sequences: {e}")
            sys.exit(1)
    
    # -------------------------------------------------------------------------
    # Install conda environments using Snakemake
    # -------------------------------------------------------------------------
    print("\n" + "="*70)
    print("[Setup] Conda environments...")
    print("="*70)
    
    # Check for snakemake
    if not check_snakemake():
        print("\n[Error] Snakemake not found. Please install Snakemake first.")
        sys.exit(1)
    
    conda_prefix = workflow_path / "conda_envs"
    print(f"[Setup] Installing conda environments to: {conda_prefix}")
    
    # Create a temporary config file for Snakemake to parse the Snakefile
    # The Snakefile requires a config.yaml with specific keys to be parsed
    import tempfile
    temp_config = {
        'directories': {
            'output_dir': str(workflow_path),
            'fastq_dir': str(workflow_path)  # Dummy path, not used for env creation
        },
        'references': {
            'xf_genomes': str(ref_seqs_dest / "xf_genomes")
        }
    }
    
    temp_config_file = workflow_path / "temp_config.yaml"
    try:
        with open(temp_config_file, 'w') as f:
            yaml.dump(temp_config, f, default_flow_style=False)
        
        # Use snakemake to create all conda environments at once
        snakemake_cmd = [
            "snakemake",
            "--snakefile", str(SNAKEFILE),
            "--configfile", str(temp_config_file),
            "--use-conda",
            "--conda-prefix", str(conda_prefix),
            "--conda-create-envs-only",
            "--cores", "1",
            "--quiet"
        ]
        
        print(f"[Setup] Running: {' '.join(snakemake_cmd)}")
        print("[Setup] This may take a while...")
        
        # Redirect stdout/stderr to suppress Snakefile parsing output
        result = subprocess.run(
            snakemake_cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print("[Setup] ✓ Conda environments installed successfully")
    except subprocess.CalledProcessError as e:
        print(f"[Warning] Failed to install conda environments: {e}")
        if e.stderr:
            print(f"[Warning] Error details: {e.stderr}")
        print("[Setup] You may need to run the pipeline with --use-conda to create them.")
    finally:
        # Clean up temporary config file
        if temp_config_file.exists():
            temp_config_file.unlink()
    
    # -------------------------------------------------------------------------
    # Download Kraken2 database
    # -------------------------------------------------------------------------
    print("\n" + "="*70)
    print("[Setup] Kraken2 database (PlusPF 8GB)...")
    print("="*70)
    
    kraken_db_path = workflow_path / "databases" / "kraken2"
    
    # Check if Kraken2 database is already downloaded (look for hash.k2d file)
    kraken_hash_file = kraken_db_path / "hash.k2d"
    
    if kraken_hash_file.exists():
        print(f"[Setup] ✓ Kraken2 database already exists at: {kraken_db_path}")
        print("[Setup] Skipping download.")
    else:
        # Download using wget
        kraken2_db_url = "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20251015.tar.gz"
        
        download_cmd = f"""
        cd {kraken_db_path} && \\
        wget -c {kraken2_db_url} -O k2_pluspf_08gb.tar.gz && \\
        tar -xzf k2_pluspf_08gb.tar.gz && \\
        rm k2_pluspf_08gb.tar.gz
        """
        
        print(f"[Setup] Database will be downloaded to: {kraken_db_path}")
        print("[Setup] This may take a while depending on your internet connection...")
        
        try:
            subprocess.run(download_cmd, shell=True, check=True)
            print("[Setup] ✓ Kraken2 database downloaded successfully")
        except subprocess.CalledProcessError as e:
            print(f"[Warning] Failed to download Kraken2 database: {e}")
            print("[Setup] You can manually download the database later.")
    
    # Create configuration file for this setup
    config_file = workflow_path / "xf_capture_config.yaml"
    config_data = {
        "xf_capture": {
            "workflow_dir": str(workflow_path),
            "conda_prefix": str(conda_prefix),
            "kraken2_database": str(kraken_db_path),
            "probes": str(ref_seqs_dest / "probes.fasta"),
            "mlst": str(ref_seqs_dest / "mlst.fasta"),
            "xf_genomes": str(ref_seqs_dest / "xf_genomes"),
        }
    }
    
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)
    
    print(f"\n[Setup] Configuration saved to: {config_file}")
    
    # Save workflow_dir to user configuration for future runs
    user_config_file = save_user_config(workflow_path)
    print(f"[Setup] User config saved to:   {user_config_file}")
    
    # Print summary
    print("\n" + "="*70)
    print("[Setup] XfCapture setup completed!")
    print("="*70)
    print(f"""
Summary:
  - Reference sequences: {ref_seqs_dest}
  - Conda environments:  {conda_prefix}
  - Kraken2 database:    {kraken_db_path}
  - Workflow config:     {config_file}
  - User config:         {user_config_file}

To run the pipeline, simply use:
  xf_capture run --input_dir <fastq_dir> --output_dir <results_dir>

The workflow directory is now saved and will be used automatically.
You can override it with --workflow_dir if needed.
    """)


def generate_config(input_dir: str, output_dir: str, workflow_dir: str = None, 
                    kraken_db: str = None, threads: int = 8) -> Path:
    """
    Generate a temporary config file for the run.
    
    When workflow_dir is provided, uses the reference sequences and databases
    from that directory. Otherwise, falls back to the package's bundled references.
    """
    # Default configuration using package's bundled references
    config = {
        "directories": {
            "fastq_dir": str(Path(input_dir).resolve()),
            "output_dir": str(Path(output_dir).resolve()),
        },
        "references": {
            "probes": str(REFERENCE_SEQS_DIR / "probes.fasta"),
            "xf_genomes": str(REFERENCE_SEQS_DIR / "xf_genomes"),
        },
        "kraken2": {
            "database": kraken_db if kraken_db else "/path/to/kraken2/database",
            "memory_mapping": False,
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


def run_pipeline(input_dir: str, output_dir: str, workflow_dir: str = None,
                 kraken_db: str = None, cores: int = 16, 
                 kraken_jobs: int = 1, alignment_jobs: int = 5,
                 iqtree_jobs: int = 3, auto: bool = False,
                 dry_run: bool = False, extra_args: list = None):
    """
    Run the XfCapture pipeline.
    """
    print_header()
    print("\n[Run] Starting XfCapture pipeline...\n")
    
    # Check prerequisites
    if not check_snakemake():
        print("[Error] Snakemake not found. Please install Snakemake first.")
        sys.exit(1)
    
    # If workflow_dir not provided, try to load from user configuration
    if workflow_dir is None:
        workflow_dir = get_default_workflow_dir()
        if workflow_dir:
            print(f"[Run] Using saved workflow directory: {workflow_dir}")
        else:
            print("[Warning] No workflow directory configured.")
            print("[Warning] Run 'xf_capture set_workflow --dir <path>' first,")
            print("[Warning] or specify --workflow_dir manually.")
    
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
        threads=cores
    )
    print(f"[Run] Config file:      {config_path}")
    
    # Determine conda prefix
    conda_prefix = None
    if workflow_dir:
        conda_prefix = Path(workflow_dir) / "conda_envs"
        if conda_prefix.exists():
            print(f"[Run] Conda prefix:     {conda_prefix}")
    
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
        f"iqtree_jobs={iqtree_jobs}",
    ]
    
    if conda_prefix and conda_prefix.exists():
        snakemake_cmd.extend(["--conda-prefix", str(conda_prefix)])
    
    if dry_run:
        snakemake_cmd.append("--dry-run")
    
    if extra_args:
        snakemake_cmd.extend(extra_args)
    
    print(f"\n[Run] Executing: {' '.join(snakemake_cmd)}\n")
    print("="*70)
    
    # Phase 1: Main pipeline up to MLST
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
    
    # Check for successful samples
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
        if not auto and not dry_run:
            try:
                response = input("\n[Run] Continue with phylogenetic analysis? (y/N): ")
                if response.lower() not in ['y', 'yes']:
                    print("[Info] Pipeline terminated by user choice.")
                    return
            except EOFError:
                print("\n[Info] Non-interactive mode. Use --auto to continue automatically.")
                return
        
        # Phase 2: Phylogenetic analysis
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


def main():
    """Main entry point for XfCapture CLI."""
    parser = argparse.ArgumentParser(
        prog="xf_capture",
        description="XfCapture - Xylella fastidiosa Capture Sequencing Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Setup the workflow (install environments and databases)
  xf_capture setup --dir /path/to/workflow_data

  # Run the pipeline
  xf_capture run --input_dir /path/to/fastq --output_dir /path/to/results

  # Run with custom options
  xf_capture run --input_dir ./data --output_dir ./results \\
             --workflow_dir /path/to/workflow_data --cores 16 --auto
        """
    )
    
    # Create subparsers for commands
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # -------------------------------------------------------------------------
    # setup subcommand
    # -------------------------------------------------------------------------
    setup_parser = subparsers.add_parser(
        "setup",
        help="Setup workflow: install conda environments and download databases.",
        description="Setup the XfCapture workflow by installing conda environments "
                    "and downloading required databases."
    )
    
    setup_parser.add_argument(
        "--dir",
        required=True,
        metavar="PATH",
        help="Directory where conda environments and databases will be installed."
    )
    
    # -------------------------------------------------------------------------
    # run subcommand
    # -------------------------------------------------------------------------
    run_parser = subparsers.add_parser(
        "run",
        help="Run the XfCapture pipeline.",
        description="Execute the XfCapture analysis pipeline on FASTQ files."
    )
    
    run_parser.add_argument(
        "--input_dir",
        required=True,
        metavar="PATH",
        help="Directory containing input FASTQ files (paired-end reads)."
    )
    
    run_parser.add_argument(
        "--output_dir",
        required=True,
        metavar="PATH",
        help="Directory for pipeline output files."
    )
    
    run_parser.add_argument(
        "--workflow_dir",
        metavar="PATH",
        help="Path to the workflow directory (created with set_workflow). "
             "Contains conda environments and databases."
    )
    
    run_parser.add_argument(
        "--kraken_db",
        metavar="PATH",
        help="Path to Kraken2 database (overrides workflow_dir setting)."
    )
    
    run_parser.add_argument(
        "--cores",
        type=int,
        default=10,
        help="Number of CPU cores to use (default: 10)."
    )
    
    run_parser.add_argument(
        "--kraken_jobs",
        type=int,
        default=1,
        help="Maximum concurrent Kraken2 jobs (default: 1)."
    )
    
    run_parser.add_argument(
        "--alignment_jobs",
        type=int,
        default=5,
        help="Maximum concurrent alignment jobs (default: 5)."
    )
    
    run_parser.add_argument(
        "--iqtree_jobs",
        type=int,
        default=3,
        help="Maximum concurrent IQ-TREE jobs (default: 3)."
    )
    
    run_parser.add_argument(
        "--auto",
        action="store_true",
        help="Automatically continue with phylogenetic analysis without prompting."
    )
    
    run_parser.add_argument(
        "--dry_run", "-n",
        action="store_true",
        help="Perform a dry run (show what would be done without executing)."
    )
    
    # Parse arguments
    args, extra_args = parser.parse_known_args()
    
    # Show help if no command provided
    if args.command is None:
        print_header()
        parser.print_help()
        sys.exit(0)
    
    # Handle setup command
    if args.command == "setup":
        setup_workflow(args.dir)
        return
    
    # Handle run command
    if args.command == "run":
        run_pipeline(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            workflow_dir=args.workflow_dir,
            kraken_db=args.kraken_db,
            cores=args.cores,
            kraken_jobs=args.kraken_jobs,
            alignment_jobs=args.alignment_jobs,
            iqtree_jobs=args.iqtree_jobs,
            auto=args.auto,
            dry_run=args.dry_run,
            extra_args=extra_args if extra_args else None
        )


if __name__ == "__main__":
    main()
