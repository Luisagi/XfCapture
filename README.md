# XfCapture


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
[![Python](https://img.shields.io/badge/python-≥3.11-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)](https://github.com/Luisagi/XfCapture)


<img src="misc/logo.jpeg" align="right" width="200"/></a>

<hr>

> A scalable Snakemake pipeline for end-to-end analysis of  *Xylella fastidiosa* **T**argeted **S**equence **C**apture **E**nrichment (*Xf*-TSCE) Illumina sequencing data, from raw FASTQ files to phylogenetic inference.



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

### Step 1: System and software requirements

Before installing, ensure your system meets these requirements:

- **Operating System**: Linux or macOS
- **RAM**: Minimum 32 GB (64+ GB recommended for large datasets)
- **CPU**: Multi-core system (16+ cores recommended)
- **Disk Space**: ~8-16 GB for databases + storage for your data

The following software must be installed:

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 9.0
- [Python](https://www.python.org/) ≥ 3.11
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://github.com/mamba-org/mamba)

```bash
# optional, if not installed already
conda create -c conda-forge -c bioconda -n xfcapture snakemake

# activate environment
conda activate xfcapture
```
### Step 2: Installation

If you already have Conda/Mamba and Snakemake installed, you can run the pipeline in three commands:

```bash
# 1. Install XfCapture
git clone https://github.com/Luisagi/XfCapture.git
cd XfCapture
pip install -e .

# 2. Prepare references and databases
xf_capture setup --dir /path/to/xf_capture_db --k2-db "16Gb"

# 3. Run the pipeline
xf_capture run -i test_data/ -o results/ --cores 16
```

run `xf_capture` to see full help:

```bash
xf_capture
usage: xf_capture [-h] <command> ...

Available commands:
  setup     Prepare workflow resources (references, databases, configs)
  run       Run the XfCapture Snakemake pipeline

Typical usage:
  xf_capture setup --dir /path/to/xf_workflow
  xf_capture run -i reads/ -o results/

For help on a specific command:
  xf_capture <command> --help
```

---

## Input Requirements

### File Structure

All samples must be provided as paired-end reads located in the same directory. The pipeline automatically infers sample identities from filenames.

```bash
input_fastq_dir/
├── sampleA_R1.fastq.gz
├── sampleA_R2.fastq.gz
├── sampleB_R1.fastq.gz
└── sampleB_R2.fastq.gz
```

### Supported Naming Formats

The pipeline automatically recognizes multiple naming conventions:

| Format | Example | Sample Name |
|--------|---------|-------------|
| Standard | `ABC123_R1.fastq.gz` | ABC123 |
| CASAVA | `ABC123_S1_L001_R1_001.fastq.gz` | ABC123 |
| Simple | `ABC123_1.fastq.gz` | ABC123 |

**Supported extensions**: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`

**Important**: Sample names cannot contain: `_*#@%^/! ?&:;|<>`

---

## Output Structure

The XfCapture pipeline is organized into two main phases. First, each sample is processed independently (quality control, taxonomic classification, gene reconstruction, and MLST typing). Second, comparative analyses are performed only on samples that successfully pass the reconstruction step, including multiple sequence alignment and phylogenetic inference.

Each output directory corresponds to a key analysis stage:

- **01.pre-processing/**: Contains an HTML quality-control report of FASTQ files processed by `fastp`.
- **02.tax-classification/**: Includes taxonomic classification results from `Kraken2` and `Recentrifuge`.
- **03.probes_reconstruction/**: Stores reconstructed gene/probe sequences per sample (FASTA) and reconstruction statistics (CSV).
- **04.mlst-typing/**: Presents `MLST` typing results (`mlst_summary.csv`).
- **05.phylogenetic_trees/**: Holds alignments (FASTA), phylogenetic trees (Newick), and visualizations produced from successfully reconstructed samples.

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
│   ├── sampleA/                  # sample-level reconstructed FASTA and stats
│   └── sampleB/
├── 04.mlst-typing/               # MLST typing results
│   └── mlst_summary.csv
├── 05.phylogenetic_trees/        # Phylogenetic analysis per sample
│   ├── summary.txt
│   ├── sampleA/                  # per-sample alignments, trees, plots
│   └── sampleB/
└── logs/                         # Log files for each step (rule/sample.log)
```

---

## Usage

### 1. Set up the working directory

Create a working directory where the pipeline will store all required references, probe files, and the Kraken2 database. This directory can be located anywhere on your system.

XfCapture uses two pre-built Kraken2 databases:
- **8 GB** (default): suitable for machines with limited memory.
- **16 GB** (recommended): requires ≥64 GB RAM and provides improved classification performance.

Loading the Kraken2 database is the most memory-intensive step of the pipeline, so choose the database size according to your available RAM.

```bash
# Set up directories and config (choose 16Gb or 8Gb for the database)
xf_capture setup --dir /path/to/xf_capture_db --k2-db "16Gb"
```

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

**Notes:**

- You can manually download any Kraken2 database and place it in `databases/kraken2/`. 
See [kraken2 AWS indexes](https://benlangmead.github.io/aws-indexes/k2). 
- To include additional *Xf* [genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2371) in the phylogenetic analysis, add their FASTA (`.fna`) files to `reference_seqs/xf_genomes/`. No special naming is required, but avoid special characters in filenames.

### 2. Run the pipeline

There is a mock dataset available in the clone repository `test_data/` folder for testing purposes.

```bash
xf_capture run -i test_data/ -o output_dir/ --cores 8
```

---

## Troubleshooting

- Memory errors loading Kraken2 DB: use the 8GB database (`--k2-db "8Gb"`) or run on a machine with >=64 GB RAM. If you have limited RAM, consider running only sample-level steps first.
- Missing read pairs or unrecognized file names: ensure FASTQ names follow one of the supported patterns and that R1/R2 pairs are present in the input directory.
- Rule failures: inspect `output_dir/logs/` for per-rule logs.


If a database download fails during `xf_capture setup`, you can manually download the Kraken2 index and place it under `databases/kraken2/` as described above.


## References

### Core Tools

<small>

- **Snakemake**: 
  - Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33.

- **fastp**: 
  - Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp:  - An ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.

- **multiQC**:  
  - Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047-3048.

- **Kraken2**:  
  - Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.  
  - Lu, J., Rincon, N., Wood, D. E., Breitwieser, F. P., Pockrandt, C., Langmead, B., ... & Steinegger, M. (2022). - Metagenome analysis using the Kraken software suite. *Nature Protocols*, 17(12), 2815-2839.

- **Recentrifuge**:  
  - Martí, J. M. (2019). Recentrifuge: Robust comparative analysis and contamination removal for metagenomics. *PLOS Computational Biology*, 15(4), e1006967.

- **BWA**:
  - Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, 25(14), 1754-1760.

- **Samtools & BCFtools**:  
  - Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.

- **BLAST+**:  
  - Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics*, 10, 421.

- **SeqKit**:  
  - Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PloS one*, 11(10), e0163962.

- **MLST**:  
  - Seemann T., mlst  Github https://github.com/tseemann/mlst  
  - Jolley, K. A., & Maiden, M. C. (2010). BIGSdb: scalable analysis of bacterial genome variation at the population level. *BMC bioinformatics*, 11(1), 595.

- **MAFFT**:  
  - Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.

- **IQ-TREE**:  
  - Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, 37(5), 1530-1534.  
  - Chernomor, O., Von Haeseler, A., & Minh, B. Q. (2016). Terrace aware data structure for phylogenomic inference from supermatrices. *Systematic Biology*, 65, 997-1008.  
  - Minh, B. Q., Nguyen, M. A. T., & Von Haeseler, A. (2013). Ultrafast approximation for phylogenetic bootstrap. *Molecular Biology and Evolution*, 30(5), 1188–1195.

- **R packages**:  
  - R Core Team (2025). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/  
  - Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). *PeerJ*, 12, e16505.  
  - Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35(3), 526-528.  
  - Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. *Methods in ecology and evolution*, 8(1), 28-36.

</small>

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Authors and Contributors

- Luis F. Arias-Giraldo (ORCID: 0000-0003-4861-8064)[https://orcid.org/0000-0003-4861-8064]
- Maria P. Velasco-Amo  (ORCID: 0000-0001-7176-0435)[https://orcid.org/0000-0001-7176-0435]
- Blanca B. Landa       (ORCID: 0000-0002-9511-3731)[https://orcid.org/0000-0002-9511-3731]

## Acknowledgments

This work was funded by the European Union's Horizon Europe research and innovation programme under **BeXyl Grant Agreement 101060593**.

> We thank:
The *Xylella fastidiosa* research community and all developers of the bioinformatics tools integrated in this pipeline.

---

<p align="center">
  <sub>Made with ❤️ for the <i>Xylella fastidiosa</i> research community</sub>

