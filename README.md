# XfCapture (beta-version)


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
[![Python](https://img.shields.io/badge/python-≥3.11-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)](https://github.com/Luisagi/XfCapture)


<img src="misc/logo.jpeg" align="right" width="200"/></a>

<hr>

> A scalable Snakemake pipeline for end-to-end analysis of  *Xylella fastidiosa* **T**argeted **S**equence **C**apture **E**nrichment (*Xf*-TSCE) Illumina sequencing data, from raw FASTQ files to phylogenetic inference.



## Pipeline summary


<div align="center">
  <img src="misc/flowchart.png" width="500"/>
</div>


## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Input Requirements](#input-requirements)
- [Output Structure](#output-structure)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [Authors and Contributors](#authors-and-contributors)
- [License](#license)

---

## Quick Start

### Step 1: System and Software Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **OS** | Linux or macOS | - |
| **RAM** | 24 GB | 64+ GB |
| **CPU** | Multi-core | 16+ cores |
| **Disk Space** | ~8-16 GB (databases) + data storage | - |

**Required Software:**

- [Python](https://www.python.org/) ≥ 3.11
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://github.com/mamba-org/mamba)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 9.0

**Create a Conda environment with Snakemake:**

```bash
conda create -c conda-forge -c bioconda -n xfcapture snakemake
conda activate xfcapture
```

> **⚠️ IMPORTANT: Conda Channel Configuration**  
> Configure conda channels **before** installation to ensure reproducible dependency resolution:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```


### Step 2: Installation

**Three-command setup:**

```bash
# 1. Install XfCapture
git clone https://github.com/Luisagi/XfCapture.git
cd XfCapture
pip install -e .

# 2. Prepare references and databases
xf_capture setup --dir /path/to/xf_capture_db --k2-db "16GB"

# 3. Run the pipeline
xf_capture run -i test_data/ -o results/ --cores 16
```

**Available commands:**

```bash
xf_capture setup    # Prepare resources (references, databases, configs)
xf_capture run      # Run the pipeline
xf_capture --help   # Show detailed help
```

---

## Input Requirements

### File Structure

Provide paired-end reads in a single directory. Sample names are automatically inferred from filenames.

```bash
input_fastq_dir/
├── sample-A_R1.fastq.gz
├── sample-A_R2.fastq.gz
├── sample-B_R1.fastq.gz
└── sample-B_R2.fastq.gz
```

### Supported Naming Formats

The pipeline automatically recognizes multiple naming conventions:

| Format   | Example                         | Sample Name |
|----------|---------------------------------|-------------|
| Standard | `ABC123_R1.fastq.gz`            | ABC123      |
| CASAVA   | `ABC123_S1_L001_R1_001.fastq.gz`| ABC123      |
| Simple   | `ABC123_1.fastq.gz`             | ABC123      |

**Supported extensions**: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`

**Important**: Sample names cannot contain: `_*#@%^/! ?&:;|<>`

---

## Output Structure

**Pipeline workflow:** Individual sample processing (QC → taxonomic classification → gene reconstruction → MLST typing) followed by comparative analysis (alignment → phylogenetic inference) on successfully reconstructed samples.

| Directory | Description |
|-----------|-------------|
| **01.pre-processing/** | Quality control reports (fastp HTML) |
| **02.tax-classification/** | Taxonomic classification (Kraken2, Recentrifuge) |
| **03.probes_reconstruction/** | Reconstructed sequences (FASTA) + statistics (CSV) |
| **04.mlst-typing/** | MLST typing results |
| **05.phylogenetic_trees/** | Alignments, trees (Newick), and visualizations |

> **Note:** Phylogenetic analysis is optional, requires user confirmation, and takes longer to complete.

```bash
output_dir/
├── 01.pre-processing/            # fastp QC reports and trimmed reads
│   ├── qc_report.html
│   └── qc_report_data/           # fastp JSON/HTML assets, per-sample reports
├── 02.tax-classification/        # Taxonomic classification reports
│   ├── recentrifuge_report.html
│   ├── recentrifuge_report.xlsx
│   └── xf_taxid_2370/            # per-taxid outputs (FASTQ)
├── 03.probes_reconstruction/     # Reconstructed gene sequences & stats
│   ├── sample-A/                  # sample-level reconstructed FASTA and stats
│   └── sample-B/
├── 04.mlst-typing/               # MLST typing results
│   └── mlst_summary.csv
├── 05.phylogenetic_trees/        # Phylogenetic analysis per sample
│   ├── summary.txt
│   ├── sample-A/                  # per-sample alignments, trees, plots
│   └── sample-B/
└── logs/                         # Log files for each step (rule/sample.log)
```

---

## Usage

### 1. Setup Workflow Directory

**Kraken2 Database Options:**

| Database | RAM Required | Notes |
|----------|--------------|-------|
| **8 GB** (default) | ~8 GB | For limited memory systems |
| **16 GB** (recommended) | ≥24 GB | Better performance; use `--k2-mapping-memory` to reduce RAM |

> **Note:** Database loading is the most memory-intensive step. Choose according to your available RAM.

```bash
xf_capture setup --dir /path/to/xf_capture_db --k2-db "16GB"
```

**Options:**
- `--dir`: Workflow directory location
- `--k2-db`: Database size (`8GB` or `16GB`, default: `8GB`)

This command creates the following structure:

```
/path/to/xf_capture_db/
├── conda_envs/
├── databases/
│   └── kraken2/
│       ├── k2_pluspfp_16_GB
│       └── k2_pluspfp_08_GB
├── reference_seqs/
│   ├── probes.fasta
│   └── xf_genomes/
└── xf_capture_config.yaml
```

**Customization:**
- **Custom Kraken2 database**: Download from [AWS indexes](https://benlangmead.github.io/aws-indexes/k2) and place in `databases/kraken2/`
- **Additional genomes**: Add `.fna` files to `reference_seqs/xf_genomes/` ([download here](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2371))

### 2. Run the Pipeline

**Basic usage:**

```bash
xf_capture run -i test_data/ -o results/ --cores 16
```

**Key Options:**

| Option | Description | Default |
|--------|-------------|----------|
| `-i, --input-dir` | Input FASTQ directory | required |
| `-o, --output-dir` | Output directory | required |
| `--cores` | Total CPU cores | 16 |
| `--k2-mapping-memory` | Reduce Kraken2 RAM usage | False |
| `--no-auto` | Require confirmation for phylogeny | False |

**Resource Allocation:**
- `--kraken-jobs`: Parallel Kraken2 jobs (default: 1)
- `--alignment-jobs`: Parallel alignment jobs (default: 4)
- `--iqtree-jobs`: Parallel IQ-TREE jobs (default: 2)
- `--iqtree-threads`: Threads per IQ-TREE job (default: 8)
- `--kraken-threads`: Threads per Kraken2 job (default: 8)

> **Tip:** Use `test_data/` for testing the pipeline

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| **Memory errors (Kraken2)** | Use 8GB database (`--k2-db "8GB"`) or `--k2-mapping-memory` flag |
| **Missing/unrecognized files** | Verify FASTQ naming matches [supported formats](#supported-naming-formats) and R1/R2 pairs exist |
| **Rule failures** | Check logs in `output_dir/logs/` |
| **Database download fails** | Manually download from [AWS](https://benlangmead.github.io/aws-indexes/k2) and place in `databases/kraken2/` |


## References

### Core Tools

<small>

- **Snakemake**: 
  - Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. *F1000Research*, *10*, 33. https://doi.org/10.12688/f1000research.29032.2

- **fastp**: 
  - Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, *34*(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

- **MultiQC**:  
  - Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, *32*(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354

- **Kraken2**:  
  - Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, *20*(1), 257. https://doi.org/10.1186/s13059-019-1891-0
  - Lu, J., Rincon, N., Wood, D. E., Breitwieser, F. P., Pockrandt, C., Langmead, B., Salzberg, S. L., & Steinegger, M. (2022). Metagenome analysis using the Kraken software suite. *Nature Protocols*, *17*(12), 2815-2839. https://doi.org/10.1038/s41596-022-00738-y

- **Recentrifuge**:  
  - Martí, J. M. (2019). Recentrifuge: Robust comparative analysis and contamination removal for metagenomics. *PLOS Computational Biology*, *15*(4), e1006967. https://doi.org/10.1371/journal.pcbi.1006967

- **BWA**:
  - Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, *25*(14), 1754-1760. https://doi.org/10.1093/bioinformatics/btp324

- **Samtools & BCFtools**:  
  - Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, *10*(2), giab008. https://doi.org/10.1093/gigascience/giab008

- **BLAST+**:  
  - Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics*, *10*, 421. https://doi.org/10.1186/1471-2105-10-421

- **SeqKit**:  
  - Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*, *11*(10), e0163962. https://doi.org/10.1371/journal.pone.0163962

- **MLST**:  
  - Seemann, T. (2024). mlst: Scan contig files against PubMLST typing schemes. https://github.com/tseemann/mlst
  - Jolley, K. A., & Maiden, M. C. (2010). BIGSdb: Scalable analysis of bacterial genome variation at the population level. *BMC Bioinformatics*, *11*(1), 595. https://doi.org/10.1186/1471-2105-11-595

- **MAFFT**:  
  - Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. *Molecular Biology and Evolution*, *30*(4), 772–780. https://doi.org/10.1093/molbev/mst010

- **IQ-TREE**:  
  - Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, *37*(5), 1530-1534. https://doi.org/10.1093/molbev/msaa015
  - Chernomor, O., Von Haeseler, A., & Minh, B. Q. (2016). Terrace aware data structure for phylogenomic inference from supermatrices. *Systematic Biology*, *65*(6), 997-1008. https://doi.org/10.1093/sysbio/syw037
  - Minh, B. Q., Nguyen, M. A. T., & Von Haeseler, A. (2013). Ultrafast approximation for phylogenetic bootstrap. *Molecular Biology and Evolution*, *30*(5), 1188–1195. https://doi.org/10.1093/molbev/mst024

- **R Packages**:  
  - R Core Team. (2025). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. https://www.R-project.org/
  - Revell, L. J. (2024). phytools 2.0: An updated R ecosystem for phylogenetic comparative methods (and other things). *PeerJ*, *12*, e16505. https://doi.org/10.7717/peerj.16505
  - Paradis, E., & Schliep, K. (2019). ape 5.0: An environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*, *35*(3), 526-528. https://doi.org/10.1093/bioinformatics/bty633
  - Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: An R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. *Methods in Ecology and Evolution*, *8*(1), 28-36. https://doi.org/10.1111/2041-210X.12628

</small>

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Authors and Contributors
- [ ... ]

## Acknowledgments

This work was funded by the European Union's Horizon Europe research and innovation programme under **BeXyl Grant Agreement 101060593**.

**Special thanks to:**
- The *Xylella fastidiosa* research community
- Developers of the bioinformatics tools integrated in this pipeline

---

<p align="center">
  <sub>Made with ❤️ for the <i>Xylella fastidiosa</i> research community</sub>

