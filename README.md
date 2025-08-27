# XfCapture

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.4.7-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)

A scalable Snakemake pipeline for analyzing *Xylella fastidiosa* targeted sequence capture enrichment (*Xf*-TSCE) sequencing data, featuring automated quality control, targeted gene reconstruction, MLST typing, and phylogenetic analysis.

## Table of Contents

- [XfCapture](#xfcapture)
  - [Table of Contents](#table-of-contents)
  - [Authors and Contributors](#authors-and-contributors)
  - [Overview](#overview)
  - [Quick Start](#quick-start)
  - [Input Requirements](#input-requirements)
    - [File Structure](#file-structure)
    - [Naming Conventions](#naming-conventions)
  - [Output Structure](#output-structure)
  - [Prerequisites](#prerequisites)
    - [Required Software](#required-software)
    - [System Requirements](#system-requirements)
    - [Required Databases (see Database Setup)](#required-databases-see-database-setup)
  - [Database Setup](#database-setup)
    - [Kraken2 Database](#kraken2-database)
    - [Reference Genomes](#reference-genomes)
    - [Update config.yaml](#update-configyaml)
  - [Usage Examples](#usage-examples)
    - [Basic Usage](#basic-usage)
    - [Advanced Options](#advanced-options)
  - [Performance Optimization](#performance-optimization)
    - [Resource Configuration](#resource-configuration)
  - [References](#references)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)

<!-- Add a blank line before the next heading -->
  
## Authors and Contributors

- Luis F. Arias-Giraldo (ORCID: 0000-0003-4861-8064)
- Maria P. Velasco-Amo  (ORCID: 0000-0001-7176-0435)
- Blanca B. Landa       (ORCID: 0000-0002-9511-3731)

## Overview

This pipeline processes paired-end sequencing data from ***X. fastidiosa*** targeted sequence capture enrichment (*Xf*-TSCE) through multiple analysis steps:

1. **Quality Control & Trimming** (fastp) - Remove adapters and low-quality reads
2. **Taxonomic Classification** (Kraken2 + Recentrifuge) - Identify and extract *Xylella* reads
3. **Target Gene Reconstruction** (BWA + BCFtools) - Map reads to probe sequences and generate consensus
4. **MLST Typing** (MLST) - Multi-locus sequence typing for strain identification
5. **Phylogenetic Analysis** (MAFFT + IQ-TREE) - Multiple alignment and tree reconstruction

The workflow uses a **checkpoint system** to ensure only successfully reconstructed samples proceed to computationally expensive phylogenetic analysis.

## Quick Start

```bash
# 1. Install XfCapture
git clone https://github.com/Luisagi/XfCapture.git
cd XfCapture/

# 2. Configure your data paths
nano config.yaml

# 3. Run the pipeline
./run_pipeline.sh --auto
```

## Input Requirements

### File Structure

```bash
input_fastq_dir/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

### Naming Conventions

- **Standard format**: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- **CASAVA format**: `{sample}_L001_R1_001.fastq.gz` and `{sample}_L001_R2_001.fastq.gz`
- Sample names cannot contain: `_*#@%^/! ?&:;|<>`
- Supported extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`

> **Note**: The pipeline automatically extracts sample names from any supported format. For CASAVA format, sample names are extracted from the portion before the first underscore (e.g., `ABC123_L001_R1_001.fastq.gz` → sample name: `ABC123`).

## Output Structure

The pipeline generates the following main output directories:

- **01.pre-processing/**: Contains quality control reports and trimmed FASTQ files produced by fastp.
- **02.tax-classification/**: Includes Kraken2 classification reports and Recentrifuge visualizations for taxonomic profiling.
- **03.probes_reconstruction/**: Stores reconstructed target gene sequences, mapping statistics, and consensus FASTA files for each sample.
- **04.mlst-typing/**: Provides MLST typing results, including allele assignments and sequence type predictions.
- **05.phylogenetic_trees/**: Contains multiple sequence alignments, phylogenetic trees (Newick format), and tree visualizations for successful samples.

```bash
output_dir/
├── 01.pre-processing/           # Quality control results
├── 02.tax-classification/       # Taxonomic classification
├── 03.probes_reconstruction/    # Targeted genes reconstruction
├── 04.mlst-typing/              # MLST results
└── 05.phylogenetic_trees/       # Phylogenetic analysis
```

Each directory contains subfolders and files organized by sample name and analysis step, allowing easy tracking and downstream use of results.

## Prerequisites

### Required Software

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 9.5.1
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://github.com/mamba-org/mamba) (recommended)

```bash
# optional, if not installed already
conda create -c conda-forge -c bioconda -n snakemake snakemake
# activate Snakemake environment
conda activate snakemake
```

### System Requirements

- **RAM**: Minimum 32 GB (64+ GB recommended for large datasets)
- **CPU**: Multi-core system (16+ cores recommended)

### Required Databases (see [Database Setup](#database-setup))

- Kraken2 PlusPFP database (8 or 16 GB)
- *Xylella fastidiosa* reference genomes
- Probe sequences (included: `reference_seqs/probes.fasta`)

## Database Setup

### Kraken2 Database

```bash
# Download Kraken2 PlusPFP 16GB database (recommended)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16_GB_20250714.tar.gz
tar -xzf k2_pluspfp_16_GB_20250714.tar.gz

# For smaller systems, use 8GB version
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08_GB_20250714.tar.gz
```

### Reference Genomes

```bash
# Download Xylella fastidiosa genomes from NCBI
mkdir -p reference_genomes
# Use NCBI Datasets or manual download
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2371
```

### Update config.yaml

```yaml
kraken2:
  # Path to Kraken2 PlusPFP database (choose one)
  database: "/absolute/path/to/k2_pluspfp_16_GB_20250714"
  #database: "/absolute/path/to/k2_pluspfp_08_GB_20250714"

references:
  xf_genomes: "/absolute/path/to/reference_genomes" # Directory containing X. fastidiosa ref genomes
  probes: "reference_seqs/probes.fasta" # Path to probe sequences (included by default)
```

## Usage Examples

### Basic Usage

```bash
# Check configuration and dependencies
./run_pipeline.sh --dry-run

# Run full pipeline with user confirmation
./run_pipeline.sh

# Automatic execution (no prompts)
./run_pipeline.sh --auto
```

### Advanced Options

```bash
# Custom resource allocation
./run_pipeline.sh --cores 10 --auto

# Force re-run specific steps
./run_pipeline.sh --forcerun verify_reconstruction

# Run only Phase 1 (up to MLST)
snakemake --cores 10 --use-conda

# Run only Phase 2 (phylogenetic analysis)
snakemake phylogeny_phase --cores 10 --use-conda --resources iqtree_jobs=3 alignment_jobs=3
```

## Performance Optimization

### Resource Configuration

```bash
# Edit run_pipeline.sh with your preferred editor to customize these values:
# nano run_pipeline.sh
# vim run_pipeline.sh
# code run_pipeline.sh

# For HPC systems (128+ GB RAM, 32+ cores)
CORES=10
KRAKEN_MAX_JOBS=2
ALIGNMENT_MAX_JOBS=4
IQTREE_MAX_JOBS=4

# For high-memory systems (64+ GB RAM, 16-24 cores)
CORES=5
KRAKEN_MAX_JOBS=1
ALIGNMENT_MAX_JOBS=3
IQTREE_MAX_JOBS=3
```

<!-- ## Citation

If you use this pipeline, please cite:

```
[Your pipeline citation]

``` -->

## References

- **Snakemake**: 
  - Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33.
- **fastp**: 
  - Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
- **Kraken2**: 
  - Lu, J., Rincon, N., Wood, D. E., Breitwieser, F. P., Pockrandt, C., Langmead, B., ... & Steinegger, M. (2022). Metagenome analysis using the Kraken software suite. *Nature Protocols*, 17(12), 2815-2839.  
  - Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.
- **Recentrifuge**: 
  - Martí, J. M. (2019). Recentrifuge: Robust comparative analysis and contamination removal for metagenomics. *PLOS Computational Biology*, 15(4), e1006967.
- **BWA**: 
  - Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, 25(14), 1754-1760.
- **Samtools & BCFtools**: 
  - Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.
- **BLAST+**: 
  - Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics*, 10, 421.
- **SeqKit**:
  - Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), e0163962.
- **MLST-typing**:
  - https://github.com/tseemann/mlst
  - Jolley, K. A., & Maiden, M. C. (2010). BIGSdb: scalable analysis of bacterial genome variation at the population level. BMC bioinformatics, 11(1), 595.
- **MAFFT**: 
  - Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780.
- **IQ-TREE**: 
  - Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, 37(5), 1530-1534.  
  - Chernomor, O., Von Haeseler, A., & Minh, B. Q. (2016). Terrace aware data structure for phylogenomic inference from supermatrices. *Systematic Biology*, 65, 997-1008.  
  - Minh, B. Q., Nguyen, M. A. T., & Von Haeseler, A. (2013). Ultrafast approximation for phylogenetic bootstrap. *Molecular Biology and Evolution*, 30(5), 1188–1195.
- **R tools**: 
  - R Core Team (2025). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. [https://www.R-project.org/](https://www.R-project.org/)
  - Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ, 12, e16505.
  - Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526-528.
  - Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in ecology and evolution, 8(1), 28-36.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The *Xylella fastidiosa* research community (**BeXyl  Grant Agreement 101060593**)
- Snakemake workflow management system
- All the bioinformatics tools integrated in this pipeline
