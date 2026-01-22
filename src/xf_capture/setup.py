#!/usr/bin/env python3

from pathlib import Path
import sys
import shutil
import subprocess
import yaml
import tarfile
import gzip
import hashlib

# Rutas internas del paquete
PACKAGE_DIR = Path(__file__).resolve().parent
WORKFLOWS_DIR = PACKAGE_DIR / "workflows"
RESOURCES_DIR = PACKAGE_DIR / "resources"

# Configuración del usuario
USER_CONFIG_DIR = Path.home() / ".xf_capture"
USER_CONFIG_FILE = USER_CONFIG_DIR / "config.yaml"


# Kraken2 database configurations per size. Each entry provides the
# tarball filename, the download URL and the MD5 file URL.
KRAKEN2_DATABASES = {
    "8GB": {
        "tarball": "k2_pluspf_08_GB_20251015.tar.gz",
        "url": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20251015.tar.gz",
        "md5_url": "https://genome-idx.s3.amazonaws.com/kraken/pluspf_08_GB_20251015/pluspf_08_GB.md5",
    },
    "16GB": {
        "tarball": "k2_pluspf_16_GB_20251015.tar.gz",
        "url": "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20251015.tar.gz",
        "md5_url": "https://genome-idx.s3.amazonaws.com/kraken/pluspf_16_GB_20251015/pluspf_16_GB.md5",
    },
}


def save_user_config(workflow_dir: Path) -> Path:
    """
    Save the workflow directory path to user configuration file.
    This allows subsequent runs to find the workflow without --workflow-dir.
    """
    USER_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    data = {"workflow_dir": str(workflow_dir)}
    with open(USER_CONFIG_FILE, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False)
    return USER_CONFIG_FILE


def get_default_workflow_dir() -> str | None:
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


def verify_md5(file_path: Path, expected_md5: str) -> bool:
    """Verify MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest() == expected_md5


def download_kraken2_database(kraken_db_path: Path, db_key: str = "8GB") -> None:
    """
    Download Kraken2 PlusPF 8GB database with MD5 verification and retry logic.
    """
    # Check if database is already complete
    kraken_hash_file = kraken_db_path / "hash.k2d"
    kraken_taxo_file = kraken_db_path / "taxo.k2d"
    kraken_opts_file = kraken_db_path / "opts.k2d"

    db_complete = (kraken_hash_file.exists() and
                   kraken_taxo_file.exists() and
                   kraken_opts_file.exists())

    if db_complete:
        print(f"[Setup] ✓ Kraken2 database already exists at: {kraken_db_path}")
        print("[Setup] Skipping download.")
        return

    # Database URLs
    # Select database config
    db_cfg = KRAKEN2_DATABASES.get(db_key)
    if not db_cfg:
        print(f"[Error] Unknown Kraken2 DB key: {db_key}")
        return

    kraken2_db_url = db_cfg["url"]
    kraken2_md5_url = db_cfg["md5_url"]
    tarball_name = db_cfg["tarball"]
    tarball_file = kraken_db_path / tarball_name

    print(f"[Setup] Database will be downloaded to: {kraken_db_path}")
    print("[Setup] This may take a while depending on your internet connection...")

    max_retries = 3
    for attempt in range(1, max_retries + 1):
        try:
            print(f"[Setup] Download attempt {attempt}/{max_retries}...")

            # Download MD5 checksum file
            print("[Setup] Downloading MD5 checksum...")
            md5_file = kraken_db_path / "expected.md5"
            subprocess.run(
                f"wget -q {kraken2_md5_url} -O {md5_file}",
                shell=True, check=True
            )

            # Parse expected MD5 (format: "md5sum  filename").
            # Choose the checksum matching the tarball filename. If no exact
            # match is found, fall back to any .tar.gz entry or the last line.
            expected_md5 = None
            md5_text = md5_file.read_text()
            for line in md5_text.splitlines():
                parts = line.strip().split()
                if len(parts) < 1:
                    continue
                # Some md5 files use two columns: <md5> <filename>
                if len(parts) >= 2:
                    fname = parts[1]
                else:
                    # If only one token, skip
                    continue

                # Normalize filename (strip ./ or paths)
                fname_norm = fname.replace("./", "").split("/")[-1]
                if fname_norm == tarball_name or fname_norm.endswith(tarball_name):
                    expected_md5 = parts[0]
                    break

            if expected_md5 is None:
                # fallback: pick any entry that ends with .tar.gz
                for line in md5_text.splitlines():
                    parts = line.strip().split()
                    if len(parts) >= 2 and parts[1].endswith('.tar.gz'):
                        expected_md5 = parts[0]
                        break

            if expected_md5 is None:
                # final fallback: use first token of last non-empty line
                for line in reversed(md5_text.splitlines()):
                    if line.strip():
                        expected_md5 = line.strip().split()[0]
                        break

            print(f"[Setup] Expected MD5: {expected_md5}")

            # Download database (resume if partial download exists)
            print("[Setup] Downloading database (this may take a while)...")
            # Download database (resume if partial download exists).
            # Save using the canonical tarball filename so MD5 checks match.
            print(f"[Setup] Downloading {tarball_name} ...")
            subprocess.run(
                f"cd {kraken_db_path} && wget -c {kraken2_db_url} -O {tarball_name}",
                shell=True, check=True
            )

            # Verify MD5 checksum
            print("[Setup] Verifying download integrity (MD5 checksum)...")
            if verify_md5(tarball_file, expected_md5):
                print("[Setup] ✓ MD5 checksum verified successfully")

                # Extract the database
                print("[Setup] Extracting database...")
                subprocess.run(
                    f"cd {kraken_db_path} && tar -xzf {tarball_name}",
                    shell=True, check=True
                )

                # Clean up
                if tarball_file.exists():
                    tarball_file.unlink()
                if md5_file.exists():
                    md5_file.unlink()

                print("[Setup] ✓ Kraken2 database downloaded and extracted successfully")
                break
            else:
                print(f"[Warning] MD5 checksum verification failed (attempt {attempt})")
                # Remove corrupted file
                if tarball_file.exists():
                    tarball_file.unlink()

                if attempt == max_retries:
                    print("[Error] Failed to download valid Kraken2 database after all retries")
                    print("[Setup] You can manually download the database later.")

        except subprocess.CalledProcessError as e:
            print(f"[Warning] Download failed (attempt {attempt}): {e}")
            if attempt == max_retries:
                print("[Error] Failed to download Kraken2 database after all retries")
                print("[Setup] You can manually download the database later.")


def setup_workflow(workflow_dir: str, k2_db: str = "8GB") -> None:
    """
    Setup the XfCapture workflow:
    - Create directory structure
    - Extract reference sequences
    - Download Kraken2 database
    - Save configuration

    Args:
        workflow_dir: Directory where the workflow will be installed
        k2_db: Kraken2 database size to download ("8GB" or "16GB", default: "8GB")

    Skips components that are already installed.
    """
    workflow_path = Path(workflow_dir).resolve()

    print("\n" + "="*70)
    print("[Setup] Starting XfCapture workflow setup")
    print("="*70 + "\n")

    # Directorios principales
    conda_prefix = workflow_path / "conda_envs"
    kraken_db = workflow_path / "databases" / "kraken2"
    ref_seqs = workflow_path / "reference_seqs"

    for d in (conda_prefix, kraken_db, ref_seqs):
        d.mkdir(parents=True, exist_ok=True)
        print(f"[Setup] Created directory: {d}")

    # --------------------------------------------------
    # Referencias
    # --------------------------------------------------
    print("\n" + "="*70)
    print("[Setup] Reference sequences...")
    print("="*70)

    tarball = RESOURCES_DIR / "ref_seqs.tar.gz"
    xf_genomes = ref_seqs / "xf_genomes"
    probes = ref_seqs / "probes.fasta"

    if probes.exists() and xf_genomes.exists() and any(xf_genomes.glob("*.fna")):
        print("[Setup] ✓ Reference sequences already extracted. Skipping.")
    else:
        if not tarball.exists():
            print(f"[Error] Missing reference archive: {tarball}")
            print("[Error] Please ensure the package was installed correctly.")
            sys.exit(1)

        print(f"[Setup] Source archive: {tarball}")
        print(f"[Setup] Destination:    {ref_seqs}")
        print("[Setup] Extracting reference sequences...")

        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(ref_seqs)
        print("[Setup] ✓ Reference sequences extracted")

        # Decompress genome files
        if xf_genomes.exists():
            gz_files = list(xf_genomes.glob("*.gz"))
            if gz_files:
                print(f"[Setup] Decompressing {len(gz_files)} genome files...")
                for gz_file in gz_files:
                    output_file = gz_file.with_suffix("")
                    with gzip.open(gz_file, 'rb') as f_in:
                        with open(output_file, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    gz_file.unlink()
                print(f"[Setup] ✓ {len(gz_files)} genome files decompressed")

    # --------------------------------------------------
    # Conda environments info
    # --------------------------------------------------
    print("\n" + "="*70)
    print("[Setup] Conda environments...")
    print("="*70)
    print(f"[Setup] Conda environments will be installed to: {conda_prefix}")
    print("[Setup] Environments will be created on first 'xf_capture run'")

    # --------------------------------------------------
    # Kraken2 database
    # --------------------------------------------------
    print("\n" + "="*70)
    print(f"[Setup] Kraken2 database (PlusPF {k2_db})...")
    print("="*70)

    download_kraken2_database(kraken_db, db_key=k2_db)

    # --------------------------------------------------
    # Configuración del workflow
    # --------------------------------------------------
    config = {
        "xf_capture": {
            "workflow_dir": str(workflow_path),
            "conda_prefix": str(conda_prefix),
            "kraken2_database": str(kraken_db),
            "probes": str(probes),
            "xf_genomes": str(xf_genomes),
        }
    }

    config_file = workflow_path / "xf_capture_config.yaml"
    with open(config_file, "w") as fh:
        yaml.dump(config, fh, default_flow_style=False)

    user_cfg = save_user_config(workflow_path)

    # Print summary
    print("\n" + "="*70)
    print("[Setup] XfCapture setup completed!")
    print("="*70)
    print(f"""
Summary:
  - Reference sequences: {ref_seqs}
  - Conda environments:  {conda_prefix}
  - Kraken2 database:    {kraken_db}
  - Workflow config:     {config_file}
  - User config:         {user_cfg}

To run the pipeline, use:
  xf_capture run -i <fastq_dir> -o <results_dir>
    """)
