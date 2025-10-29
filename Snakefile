# Modules:
import os
import sys  
import re

# Import custom functions
from utils.auxiliar_functions import (
    get_samples, get_r1, get_r2, get_fasta_files, get_gene_names,
    get_alignment_outputs, get_phylo_outputs, print_fasta_headers,
    get_successful_samples 
)

# Configuration file:
configfile: 'config.yaml'

# Path to DBs and input/output directories in config file:
workdir: config['directories']['output_dir']
INPUT_DIR = config['directories']['fastq_dir']
XF_REFS = config['references']['xf_genomes']

# Define not allowed characters
not_allowed_chars = set("_*#@%^/! ?&:;|<>")

# Valid extensions
valid_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

# --------------------------------------------------------------------------------------
# Get the samples from the input directory

SAMPLES_DICT = get_samples(INPUT_DIR, valid_extensions, not_allowed_chars)
SAMPLES = SAMPLES_DICT.keys()

# Print the samples found
print(f"\nTotal samples found: {len(SAMPLES)}")

for sample, (r1, r2) in SAMPLES_DICT.items():
    print(f"  - {sample}:\n      R1 → {r1}\n      R2 → {r2}")

# Get the list of FASTA files in the reference directory
XF_DICT = get_fasta_files(XF_REFS)
XF_GENOMES = XF_DICT.keys()

# print(f"\nTotal Xylella FASTA files found: {len(XF_DICT)}")
# print(f"Reference genomes: {list(XF_DICT.keys())}")

# --------------------------------------------------------------------------------------
# Helper functions for rules

def get_r1_wrapper(wildcards):
    """Wrapper function for get_r1 to pass SAMPLES_DICT"""
    return get_r1(wildcards, SAMPLES_DICT)

def get_r2_wrapper(wildcards):
    """Wrapper function for get_r2 to pass SAMPLES_DICT"""
    return get_r2(wildcards, SAMPLES_DICT)

# --------------------------------------------------------------------------------------
# Rules:

# Fase 1: Pipeline principal hasta MLST (incluye checkpoints)
rule all:
  input:
    fastp_json = expand("01.pre-processing/{sample}_fastp.json", sample = SAMPLES),
    fastp_html = expand("01.pre-processing/{sample}_fastp.html", sample = SAMPLES),
    rcf_report = "02.tax-classification/recentrifuge_report.html",
    xf_reads_R1 = expand("02.tax-classification/xf_taxid_2370/{sample}_R1.fastq.gz", sample = SAMPLES),
    xf_reads_R2 = expand("02.tax-classification/xf_taxid_2370/{sample}_R2.fastq.gz", sample = SAMPLES),
    bam_stats = expand("03.probes_reconstruction/stats/{sample}_mapped_stats.tsv", sample = SAMPLES),
    consensus_filtered = expand("03.probes_reconstruction/consensus/{sample}_cns_filtered.fasta", sample = SAMPLES),
    recon_stats = expand("03.probes_reconstruction/stats/{sample}_recon_stats.tsv", sample = SAMPLES),
    summary_stats =  "03.probes_reconstruction/reconstruction_summary.xlsx",
    unmapped_R1 = expand("03.probes_reconstruction/unmapped/{sample}_unmapped_R1.fastq.gz", sample = SAMPLES),
    unmapped_R2 = expand("03.probes_reconstruction/unmapped/{sample}_unmapped_R2.fastq.gz", sample = SAMPLES),
    mlst = "04.mlst-typing/mlst_identification.csv",
    checkpoint = expand("03.probes_reconstruction/checkpoint/{sample}_reconstruction_check.txt", sample = SAMPLES)
    
    
    #pre_alignment =  lambda wildcards: expand("05.phylogenetic_trees/{sample}/genes_ids.txt", sample = get_successful_samples()),
    #alignment_summary = lambda wildcards: expand("05.phylogenetic_trees/{sample}/alignment_summary.txt", sample = get_successful_samples()),
    #iqtree_tree = lambda wildcards: expand("05.phylogenetic_trees/{sample}/consensus_tree/{sample}.contree", sample = get_successful_samples()),
    #iqtree_log = lambda wildcards: expand("05.phylogenetic_trees/{sample}/consensus_tree/{sample}.log", sample = get_successful_samples()),

# --------------------------------------------------------------------------------------
# Quality control and trimming of FASTQ files using fastp
# This rule processes paired-end reads, trimming adapters and low-quality bases.
# It generates trimmed FASTQ files, HTML, and a JSON report.

rule fastq:
    input:
      r1 = get_r1_wrapper,
      r2 = get_r2_wrapper
    output:
      r1 = temp("01.pre-processing/{sample}_trim_R1.fastq"),
      r2 = temp("01.pre-processing/{sample}_trim_R2.fastq"),
      json = "01.pre-processing/{sample}_fastp.json",
      html = "01.pre-processing/{sample}_fastp.html"
    resources:
        cpus = 8
    conda:
        "envs/fastp.yaml"
    message:
        "--- Fastp: Filter and trim reads. ---"
    log:
        "logs/fastp_{sample}.log"
    shell:
        """
        fastp --detect_adapter_for_pe --length_required 50 \
        --cut_front --cut_right \
        -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        --thread {resources.cpus} \
        -j {output.json} -h {output.html} > {log} 2>&1
        """

# --------------------------------------------------------------------------------------
# Taxonomic classification of reads using Kraken2
# The rule uses paired-end reads and requires a pre-built Kraken2 database.
# The database path is specified in the configuration file.

rule kraken2:
    input:
        r1 = "01.pre-processing/{sample}_trim_R1.fastq",
        r2 = "01.pre-processing/{sample}_trim_R2.fastq",
    output:
        k2_krk = temp("02.tax-classification/{sample}.krk"),
        k2_rep = temp("02.tax-classification/{sample}.txt")
    params:
        db = config['kraken2']['database'],  # Use the 16GB database
        memory_mapping = "--memory-mapping" if config['kraken2'].get('memory_mapping', False) else ""
    resources:
        cpus = config['kraken2']['threads'],
        kraken_jobs = 1  # Limit to 1 job for Kraken2
    conda:
        "envs/kraken2.yaml"
    message:
        "--- Kraken2: Taxonomic classification of reads. ---"
    log:
        "logs/kraken2_{sample}.log"
    shell:
        """
        kraken2 --db {params.db} --threads {resources.cpus} {params.memory_mapping} \
        --output {output.k2_krk} --report {output.k2_rep} \
        --paired {input.r1} {input.r2} > {log} 2>&1
        """

# Download taxonomic classification files for Recentrifuge
# This rule downloads the taxdump files required for Recentrifuge.

rule download_taxdump:
    output:
        tax_dump = temp(directory("02.tax-classification/tax_dump"))
    conda:
        "envs/recentrifuge.yaml"
    message:
        "--- Download taxdump files for Recentrifuge. ---"
    log:
        "logs/download_taxdump.log"
    shell:
        """
        retaxdump --nodespath {output.tax_dump} > {log} 2>&1
        rm taxdmp.zip
        """

# This rule runs Recentrifuge to generate a taxonomic classification report
# using the Kraken2 output files. It requires the taxdump files downloaded in the previous rule

rule recentrifuge:
    input:
        k2_rep = expand("02.tax-classification/{sample}.krk", sample = SAMPLES),
        tax_dump = "02.tax-classification/tax_dump"
    output:
        rcf_rep = "02.tax-classification/recentrifuge_report.html",
    params:
        summary = "AVOID"
    conda:
        "envs/recentrifuge.yaml"
    message:
        "--- Recentrifuge: Generate taxonomic classification report. ---"
    log:
        "logs/recentrifuge.log"
    shell:
        """
        rcf --nodespath {input.tax_dump} --kraken 02.tax-classification/ \
        --outprefix {output.rcf_rep} --summary {params.summary} --avoidcross > {log} 2>&1
        """

# --------------------------------------------------------------------------------------
# Extracting Xylella fastidiosa reads from Kraken2 output
# This rule extracts reads classified as Xylella fastidiosa from the Kraken2 report.
# It uses the `extract_kraken_reads.py` script to filter reads based on taxonomic ID.
# The taxonomic ID for Xylella fastidiosa is specified in the parameters

rule extract_xf_reads:
    input:
        r1 = "01.pre-processing/{sample}_trim_R1.fastq",
        r2 = "01.pre-processing/{sample}_trim_R2.fastq",
        k2_krk = "02.tax-classification/{sample}.krk",
        k2_rep = "02.tax-classification/{sample}.txt"
    output:
        r1 = "02.tax-classification/xf_taxid_2370/{sample}_R1.fastq.gz",
        r2 = "02.tax-classification/xf_taxid_2370/{sample}_R2.fastq.gz",
    params:
        taxID = 2370, # Xf taxID
        tmp_r1 = "02.tax-classification/xf_taxid_2370/{sample}_R1.fastq",
        tmp_r2 = "02.tax-classification/xf_taxid_2370/{sample}_R2.fastq"
    conda:
        "envs/kraken2.yaml"
    message:
        "--- Krakentools: Extract Xylella fastidiosa reads. ---"
    log:
        "logs/extract_reads_{sample}.log"
    shell:
        """
        extract_kraken_reads.py -k {input.k2_krk} -r {input.k2_rep} \
        -1 {input.r1} -2 {input.r2} -t {params.taxID} \
        -o {params.tmp_r1} -o2 {params.tmp_r2} --fastq-output --include-children > {log} 2>&1

        # Compress the output files
        pigz -c {params.tmp_r1} > {output.r1}
        pigz -c {params.tmp_r2} > {output.r2}
        
        # Remove temporary files
        rm {params.tmp_r1} {params.tmp_r2}
        """

# --------------------------------------------------------------------------------------
# Probes reconstruction mapping a reference fasta

# Build BWA index for probes
# This rule builds a BWA index for the probes sequences specified in the configuration file.

rule build_probes:
    input: 
        fasta = config['references']['probes']
    output:
        index1 = temp("03.probes_reconstruction/probes.amb"),
        index2 = temp("03.probes_reconstruction/probes.ann"),
        index3 = temp("03.probes_reconstruction/probes.bwt"),
        index4 = temp("03.probes_reconstruction/probes.pac"),
        index5 = temp("03.probes_reconstruction/probes.sa")
    params:
        basename = "03.probes_reconstruction/probes",
    conda:
        "envs/mapping.yaml"
    message:
        "--- BWA: Build probes index for mapping. ---"
    log:
        "logs/build_bwa_index_probes.log"
    shell:
        """
        bwa index -p {params.basename} {input.fasta} > {log} 2>&1
        """

# Map reads to probes using BWA
# This rule maps the Xylella fastidiosa reads to the probes using BWA MEM
# It requires the BWA index files generated in the previous rule.

rule bwa_mapping:
    input:
        index1 = "03.probes_reconstruction/probes.amb",
        index2 = "03.probes_reconstruction/probes.ann",
        index3 = "03.probes_reconstruction/probes.bwt",
        index4 = "03.probes_reconstruction/probes.pac",
        index5 = "03.probes_reconstruction/probes.sa",
        r1 = "02.tax-classification/xf_taxid_2370/{sample}_R1.fastq.gz",
        r2 = "02.tax-classification/xf_taxid_2370/{sample}_R2.fastq.gz"
    output:
        bam_sorted = temp("03.probes_reconstruction/{sample}_mapped.sorted.bam")
    params:
        db = "03.probes_reconstruction/probes",
        cpus = 4
    conda:
        "envs/mapping.yaml"
    message:
        "--- BWA: Map reads to probes. ---"
    log:
        "logs/bwa_mapping_{sample}.log"
    shell:
        """
        bwa mem -t {params.cpus} {params.db} {input.r1} {input.r2} 2>> {log} | \
        samtools sort -@ {params.cpus} -o {output.bam_sorted} - 2>> {log}
        """

# Extract BAM statistics
# This rule generates statistics from the BAM file produced by the BWA mapping.

rule bam_stats:
    input:
        bam = "03.probes_reconstruction/{sample}_mapped.sorted.bam"
    output:
        stats = "03.probes_reconstruction/stats/{sample}_mapped_stats.tsv",
        bai = temp("03.probes_reconstruction/{sample}_mapped.sorted.bam.bai")
    conda:
        "envs/mapping.yaml"
    message:
        "--- Samtools: Generate BAM statistics. ---"
    shell:
        """
        samtools index {input.bam}
        samtools idxstats {input.bam} | \
        sed '1i #Refsequence\tSequence_length\tMapped_read_segments\tUnmapped_read_segments' \
        > {output.stats} 
        """

# Extract Xf unmapped reads with the probes from the BAM file
rule extract_non_mapped_reads:
    input:
        bam = "03.probes_reconstruction/{sample}_mapped.sorted.bam"
    output:
        unmapped_R1 = "03.probes_reconstruction/unmapped/{sample}_unmapped_R1.fastq.gz",
        unmapped_R2 = "03.probes_reconstruction/unmapped/{sample}_unmapped_R2.fastq.gz"
    params:
        cpus = 4
    conda:
        "envs/mapping.yaml"
    message:
        "--- Samtools: Extract Xf unmapped reads. ---"
    log:
        "logs/extract_unmapped_reads_{sample}.log"
    shell:
        """
        samtools sort -n {input.bam} 2>> {log} | \
        samtools fastq -@ {params.cpus} -f 4 -F 264 -c 1 \
        -1 {output.unmapped_R1} -2 {output.unmapped_R2} 2>> {log}
        """



# Call variants and generate consensus sequence
# This rule uses BCFtools to call variants from the BAM file and generate a consensus sequence
# It requires the BAM file produced by the BWA mapping and the reference FASTA file.

rule bcftools_call:
    input:
        bam = "03.probes_reconstruction/{sample}_mapped.sorted.bam"
    output:
        vcf = temp("03.probes_reconstruction/consensus/{sample}.vcf.gz"),
        csi = temp("03.probes_reconstruction/consensus/{sample}.vcf.gz.csi"),
        bed = temp("03.probes_reconstruction/consensus/{sample}_lowcov.bed"),
        consensus = "03.probes_reconstruction/consensus/{sample}_cns.fasta",
        consensus_filtered = "03.probes_reconstruction/consensus/{sample}_cns_filtered.fasta"
    params:
        fasta = config['references']['probes'],
        ploidy = 1,  # Assuming Xylella fastidiosa is haploid
        min_depth = 3  # Minimum number of reads for variant calling
    conda:
        "envs/mapping.yaml"
    message:
        "--- BCFtools: Call variants and generate consensus sequence. ---"
    log:
        "logs/call_variants_{sample}.log"
    shell:
        """
        # Generate VCF
        bcftools mpileup -Ou -f {params.fasta} {input.bam} 2>> {log} | \
        bcftools call -c --ploidy {params.ploidy} -Oz -o {output.vcf} 2>> {log}
        bcftools index {output.vcf}

        # BED mask for low coverage regions
        samtools depth -aa {input.bam} | \
        awk -v min_depth={params.min_depth} '$3 < min_depth {{print $1"\t"$2-1"\t"$2}}' > {output.bed} 2>> {log}

        # Generate FASTA consensus sequence with Ns for low coverage regions
        bcftools consensus -f {params.fasta} -m {output.bed} {output.vcf} > {output.consensus} 2>> {log}

        # Filter out sequences with Ns
        seqkit grep -s -v -r -p "N" {output.consensus} > {output.consensus_filtered}
        """

# Statistics for consensus sequence reconstruction
# This rule calculates the reconstruction statistics for the consensus sequence generated in the previous step.
# It computes the length of the sequence, the number of N bases, and the reconstruction percentage
# It also checks if the reconstruction percentage is below 90% and logs a warning if so

rule reconstruction_stats:
    input:
        consensus = "03.probes_reconstruction/consensus/{sample}_cns.fasta"
    output:
        recon_stats = "03.probes_reconstruction/stats/{sample}_recon_stats.tsv"
    message:
        "--- Reconstruction statistics for consensus sequence. ---"
    log:
        "logs/reconstruction_stats_{sample}.log"
    shell:
        """
        echo -e "#Refsequence\tSequence_length\tReconstructed_bases\tReconstruction_percent" > {output.recon_stats}
        awk '
            BEGIN {{
                OFS = "\\t"
            }}
            /^>/ {{
                if (seq != "") {{
                    len = length(seq)
                    n = gsub(/[Nn]/, "", seq)
                    recon = ((len - n) / len) * 100

                    print substr(header,2), len, len - n, recon >> "{output.recon_stats}"


                    if (recon < 100) {{
                        print "WARNING: Low reconstruction quality for gene " \
                              substr(header,2) \
                              ". Reconstruction percentage: " recon "%" >> "{log}"
                    }}

                    seq = ""
                }}
                header = $0
                next
            }}
            {{
                seq = seq $0
            }}
            END {{
                if (seq != "") {{
                    len = length(seq)
                    n = gsub(/[Nn]/, "", seq)
                    recon = ((len - n) / len) * 100

                    print substr(header,2), len, len - n, recon >> "{output.recon_stats}"

                }}
            }}
            ' {input.consensus} > {log} 2>&1 
        """

# Summary statistics for reconstruction
rule summary_reconstruction:
    input:
        recon_stats = expand("03.probes_reconstruction/stats/{sample}_recon_stats.tsv", sample = SAMPLES),
        mapped_stats = expand("03.probes_reconstruction/stats/{sample}_mapped_stats.tsv", sample = SAMPLES)
    output:
        summary_out = "03.probes_reconstruction/reconstruction_summary.xlsx"
    params:
        stats_dir = directory("03.probes_reconstruction/stats/"),
        summary_script = workflow.basedir + "/utils/reconstruction_summary.R"
    conda:
        "envs/R_tools.yaml"
    message:
        "--- Summary statistics for reconstruction ---"
    log:
        "logs/summary_reconstruction.log"
    shell:
        """
        unset R_HOME
        export R_LIBS_USER=""
        export R_LIBS=""
        Rscript {params.summary_script} {params.stats_dir} {output.summary_out} >> {log} 2>&1
        """




# --------------------------------------------------------------------------------------
# Checkpoint: Verify if consensus sequences were successfully reconstructed
# This rule checks if any consensus sequences were generated for the sample
# If no sequences are found, downstream analyses are skipped

rule verify_reconstruction:
    input:
        consensus_filtered = "03.probes_reconstruction/consensus/{sample}_cns_filtered.fasta"
    output:
        check_file = "03.probes_reconstruction/checkpoint/{sample}_reconstruction_check.txt"
    message:
        "--- Checkpoint: Verifying reconstruction success for sample {sample} ---"
    log:
        "logs/verify_reconstruction_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.check_file})
        
        # Check if consensus_filtered file has sequences
        if [ -s {input.consensus_filtered} ] && grep -q "^>" {input.consensus_filtered}; then
            echo "SUCCESS: Sample {wildcards.sample} has reconstructed sequences" > {output.check_file}
            echo "Reconstruction successful for sample {wildcards.sample}" >> {log}
            
            # Count sequences
            seq_count=$(grep -c "^>" {input.consensus_filtered})
            echo "Total sequences reconstructed: $seq_count" >> {output.check_file}
            echo "Total sequences reconstructed: $seq_count" >> {log}
            
        else
            echo "FAILED: Sample {wildcards.sample} has no reconstructed sequences" > {output.check_file}
            echo "No sequences found in {input.consensus_filtered}" >> {log}
            echo "This sample will be excluded from downstream analyses" >> {output.check_file}
            echo "Sample {wildcards.sample} excluded from downstream analysis" >> {log}
        fi
        """

# ----------------------------------------------------------------------------------------
# MLST Typing of Xylella fastidiosa strains
# This rule uses the consensus sequences generated in the previous step to perform MLST typing.

rule mlst_typing:
    input:
        consensus = expand("03.probes_reconstruction/consensus/{sample}_cns.fasta", sample = SAMPLES)
    output:
        mlst = "04.mlst-typing/mlst_identification.csv"
    params:
        cpus = 8
    conda:
        "envs/mlst.yaml"
    message:
        "--- MLST: Typing of Xylella fastidiosa strains. ---"
    log:
        "logs/mlst_typing.log"
    shell:
        """
        mlst --threads {params.cpus} --csv --nopath {input.consensus} >> {output.mlst} 2>> {log}
        """

# ----------------------------------------------------------------------------------------
# Extract genes from all Xylella fastidiosa reference genomes using BLAST

rule extract_genes_all:
    input:
        xf_genomes = [XF_DICT[fasta][0] for fasta in XF_GENOMES]
    output:
        refs = expand("05.phylogenetic_trees/refs/{fasta}/{fasta}.bed.fasta", fasta=XF_GENOMES),
        ref_dir = temp(directory("05.phylogenetic_trees/refs/")),
    params:
        query = config['references']['probes'],
        basenames = list(XF_GENOMES)  # Move the list creation here
    conda:
        "envs/extract_genes.yaml"
    message:
        "--- Extracting genes from all reference genomes ---"
    log:
        "logs/extract_xf-ref_genes_all.log"
    shell:
        """
        mkdir -p {output.ref_dir}

        # Arrays of genomes and basenames
        xf_genomes=({input.xf_genomes})
        basenames=({params.basenames})

        for i in "${{!xf_genomes[@]}}"; do
            genome="${{xf_genomes[$i]}}"
            basename="${{basenames[$i]}}"

            echo "Processing genome: $basename" >> {log}

            # Create output directory
            mkdir -p "{output.ref_dir}/$basename"

            db_name="{output.ref_dir}/$basename/$basename"
            blastn_out="{output.ref_dir}/$basename/$basename.blast.out"
            bed_file="{output.ref_dir}/$basename/$basename.bed"
            output_fasta="{output.ref_dir}/$basename/$basename.bed.fasta"

            # Create BLAST database from the Xylella genome
            makeblastdb -in "$genome" -dbtype nucl -out "$db_name" >> {log} 2>&1

            # Run BLASTn
            blastn -db "$db_name" -query {params.query} \
            -outfmt 6 -max_target_seqs 1 -max_hsps 1 \
            -qcov_hsp_perc 90 -out "$blastn_out" >> {log} 2>&1

            # Convert BLAST output to BED file
            grep -v '^#' "$blastn_out" | \
            awk 'BEGIN{{OFS="\\t"}} {{if($9<=$10) print $2,$9-1,$10,$1,"0","+"; else print $2,$10-1,$9,$1,"0","-"}}' | \
            sort > "$bed_file"

            # Extract sequences from the Xylella genome based on the BED file
            bedtools getfasta -fi "$genome" -bed "$bed_file" -name -s  2>> {log} | \
            sed "s/::.*/__$basename/g" > "$output_fasta" 2>> {log}

            # Clean up
            rm -f "$db_name".n* "$blastn_out" "$bed_file"

            echo "Completed: $basename" >> {log}
        done

        # Clean auxiliary files
        find {config[references][xf_genomes]} -type f -name "*.fai" -delete
        find $(dirname {config[references][probes]}) -type f -name "*.fai" -delete

        echo "All genomes processed successfully" >> {log}
        """

# Prepare FASTA files for alignment 
rule sort_fasta_to_alignment:
    input:
        consensus_filtered = "03.probes_reconstruction/consensus/{sample}_cns_filtered.fasta",
        ref_dir = "05.phylogenetic_trees/refs/",
        xf_refs = expand("05.phylogenetic_trees/refs/{fasta}/{fasta}.bed.fasta", fasta = XF_GENOMES),
        checkpoint = "03.probes_reconstruction/checkpoint/{sample}_reconstruction_check.txt"
    output:
        gene_ids = "05.phylogenetic_trees/{sample}/genes_ids.txt",
        gene_outdir = temp(directory("05.phylogenetic_trees/{sample}/genes"))
    conda:
        "envs/mapping.yaml"
    message:
        "--- Preparing per-gene FASTA files for alignment for sample {sample} ---"
    log:
        "logs/sort_fasta_to_alignment_{sample}.log"
    shell:
        """
        # Create output directory
        mkdir -p {output.gene_outdir}
        
        # Extract gene names that don't contain N's
        grep "^>" {input.consensus_filtered} | sed 's/^>//' > {output.gene_ids}

        # Process each gene
        for gene_id in $(cat {output.gene_ids}); do
            echo "Processing gene: ${{gene_id}}" >> {log}
            
            # Extract sequence from referece genomes
            seqkit grep -n -r -p "^${{gene_id}}__" 05.phylogenetic_trees/refs/*/*.fasta > {output.gene_outdir}/refs_${{gene_id}}.fasta

            # Extract sequence from consensus
            seqkit grep -n -r -p "^${{gene_id}}$" {input.consensus_filtered} > {output.gene_outdir}/cns_${{gene_id}}.fasta

            # Clean up header to just sample name
            sed -i "s/^>.*__/>/g" {output.gene_outdir}/refs_${{gene_id}}.fasta
            sed -i "s/^>.*/>sample_{wildcards.sample}/g" {output.gene_outdir}/cns_${{gene_id}}.fasta

            # Combine reference and consensus sequences
            cat {output.gene_outdir}/refs_${{gene_id}}.fasta {output.gene_outdir}/cns_${{gene_id}}.fasta > {output.gene_outdir}/${{gene_id}}.fasta
            
            # Remove temporary files
            rm {output.gene_outdir}/refs_${{gene_id}}.fasta {output.gene_outdir}/cns_${{gene_id}}.fasta
        done

        """

# Align and trim genes (**time consuming step**)
# This rule aligns the extracted gene sequences using MAFFT and trims them using ClipKit.

rule align_and_trim_genes:
    input:
        gene_ids = "05.phylogenetic_trees/{sample}/genes_ids.txt",
        genes_dir = "05.phylogenetic_trees/{sample}/genes",
        checkpoint = "03.probes_reconstruction/checkpoint/{sample}_reconstruction_check.txt"
    output:
        summary = "05.phylogenetic_trees/{sample}/alignment_summary.txt",
        align_dir = directory("05.phylogenetic_trees/{sample}/alignments")
    params:
        cpus = 4
    resources:
        alignment_jobs = 1  # Limit to 1 job for alignment
    conda:
        "envs/phylogeny.yaml"
    message:
        "--- Aligning genes for sample {sample} ---"
    log:
        "logs/align_genes_{sample}.log"
    shell:
        """
        # Create output directories
        mkdir -p {output.align_dir}
        
        # Initialize summary
        echo "Alignment and trimming summary for sample {wildcards.sample}" > {output.summary}
        echo "Analysis started: $(date)" >> {output.summary}
        echo "" >> {output.summary}
        
        # Process each gene
        for gene_id in $(cat {input.gene_ids}); do
            if [ -z "$gene_id" ]; then
            continue
            fi

            echo "Processing gene: $gene_id" >> {log}
            echo "Processing gene: $gene_id" >> {output.summary}

            gene_file="{input.genes_dir}/${{gene_id}}.fasta"
            aligned_file="{output.align_dir}/${{gene_id}}_aligned.fasta"
            trimmed_file="{output.align_dir}/${{gene_id}}_alignment_clipkitted.fasta"

            if [ ! -f "$gene_file" ]; then
            echo "  ERROR: Gene file $gene_file not found" >> {log}
            echo "  - Status: FAILED (file not found)" >> {output.summary}
            continue
            fi

            mafft --auto --thread {params.cpus} "$gene_file" > "$aligned_file" 2>> {log}
            clipkit "$aligned_file" -o "$trimmed_file" >> {log} 2>&1
            rm -f "$aligned_file"

            if [ -s "$trimmed_file" ]; then
            echo "  - Status: SUCCESS" >> {output.summary}
            else
            echo "  - Status: FAILED" >> {output.summary}
            fi
            
        done

        echo "Analysis completed: $(date)" >> {output.summary}
        """

# ----------------------------------------------------------------------------------------
# Phylogenetic reconstruction using IQ-TREE
# This rule constructs phylogenetic trees for each gene using IQ-TREE with bootstrap support.

rule phylo_analysis:
    input:
        gene_ids = "05.phylogenetic_trees/{sample}/genes_ids.txt",
        alignment_summary = "05.phylogenetic_trees/{sample}/alignment_summary.txt",
        checkpoint = "03.probes_reconstruction/checkpoint/{sample}_reconstruction_check.txt"
    output:
        concated_alignment = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}_concatenated_alignment.phy",
        partition_file = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}_partitions.txt",
        iqtree_log = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}.log",
        iqtree_tree = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}.contree"
    params:
        cpus = 8, # Determine the best number of cores
        model = "MFP+MERGE", # ModelFinderPlus with Merge
        tree_dir = "05.phylogenetic_trees/{sample}/consensus_tree",
        align_dir = "05.phylogenetic_trees/{sample}/alignments"
    resources:
        iqtree_jobs = 1  
    conda:
        "envs/phylogeny.yaml"
    message:
        "--- Building phylogenetic trees for sample {sample} ---"
    log:
        "logs/phylogeny_{sample}.log"
    shell:
        """
        mkdir -p {params.tree_dir}

        # create alignment file concatenated
        AMAS.py concat \
        --in-files {params.align_dir}/*_alignment_clipkitted.fasta --in-format fasta --data-type dna \
        --concat-out {output.concated_alignment} --out-format phylip \
        --concat-part {output.partition_file} --part-format nexus >> {log} 2>&1

        # IQ-TREE: Partitioned analysis for multi-gene alignments
        iqtree -s {output.concated_alignment} -p {output.partition_file} \
        -T {params.cpus} \
        -pre {params.tree_dir}/{wildcards.sample} \
        -m {params.model} -B 1000 \
        --quiet >> {log} 2>&1

        # Clean aux files
        find {params.tree_dir} -type f ! \
        \\( -name "*.contree" -o -name "*.log" -o -name "*_concatenated_alignment.phy" -o -name "*_partitions.txt" \\) \
        -delete 2>/dev/null || true
        """

# Plot phylogenetic tree
rule plot_tree:
    input:
        iqtree_tree = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}.contree"
    output:
        tree_png = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}_tree.png"
    params:
        plot_script = workflow.basedir + "/utils/plot_tree.R",
        plot_prefix = "05.phylogenetic_trees/{sample}/consensus_tree/{sample}"
    conda:
        "envs/R_tools.yaml"
    message:
        "--- Plotting phylogenetic tree for sample {sample} ---"
    log:
        "logs/plot_tree_{sample}.log"
    shell:
        """
        # Plot tree
        unset R_HOME
        export R_LIBS_USER=""
        export R_LIBS=""
        Rscript {params.plot_script} {input.iqtree_tree} {params.plot_prefix} >> {log} 2>&1
        """

# Phase 2: Phylogenetic analysis only for successful samples

rule phylogeny_phase:
    input:
        # Ensure phase 1 is complete
        mlst = "04.mlst-typing/mlst_identification.csv",
        reconstruction_check = expand("03.probes_reconstruction/checkpoint/{sample}_reconstruction_check.txt", sample = SAMPLES),
        # Global references required for phylogeny
        xf_ref_dir = "05.phylogenetic_trees/refs/",
        xf_ref_fasta = expand("05.phylogenetic_trees/refs/{fasta}/{fasta}.bed.fasta", fasta = XF_GENOMES),
        # Phylogenetic targets dynamically expanded only for successful samples
        gene_ids = lambda wildcards: expand("05.phylogenetic_trees/{sample}/genes_ids.txt", sample = get_successful_samples()),
        alignment_summary = lambda wildcards: expand("05.phylogenetic_trees/{sample}/alignment_summary.txt", sample = get_successful_samples()),
        iqtree_tree = lambda wildcards: expand("05.phylogenetic_trees/{sample}/consensus_tree/{sample}.contree", sample = get_successful_samples()),
        tree_plots = lambda wildcards: expand("05.phylogenetic_trees/{sample}/consensus_tree/{sample}_tree.png", sample = get_successful_samples())

    output:
        completion = "05.phylogenetic_trees/phylogeny_analysis_complete.txt"
    message:
        "--- Phase 2: Phylogenetic analysis for successful samples ---"
    log:
        "logs/phylogeny_phase.log"
    message:
        "--- Phase 2: Phylogenetic analysis for successful samples ---"
    shell:
        """
        echo "Phylogenetic Analysis Summary" > {output.completion}
        echo "============================" >> {output.completion}
        echo "Analysis started: $(date)" >> {output.completion}
        echo "" >> {output.completion}

        # Count samples processed
        total_samples=$(ls 03.probes_reconstruction/checkpoint/*_reconstruction_check.txt 2>/dev/null | wc -l)
        successful_count=0
        failed_count=0

        echo "Processing checkpoint results:" >> {output.completion}
        for checkpoint_file in 03.probes_reconstruction/checkpoint/*_reconstruction_check.txt; do
            [ -e "$checkpoint_file" ] || continue
            sample=$(basename "$checkpoint_file" _reconstruction_check.txt)
            
            if grep -q "SUCCESS:" "$checkpoint_file"; then
                # Count genes used for this sample
                genes_count=0
                partition_file="05.phylogenetic_trees/$sample/consensus_tree/${{sample}}_partitions.txt"

                if [ -f "$partition_file" ]; then
                    # Count partitions (genes) from the partitions file
                    genes_count=$(grep -c "charset" "$partition_file" 2>/dev/null || echo 0)
                fi

                echo "✓ Sample $sample: phylogenetic analysis completed ($genes_count targeted genes used)" >> {output.completion}
                successful_count=$((successful_count + 1))
            else
                echo "✗ Sample $sample: excluded from phylogenetic analysis" >> {output.completion}
                failed_count=$((failed_count + 1))
            fi
        done

        echo "" >> {output.completion}
        echo "Summary:" >> {output.completion}
        echo "- Total samples processed: $total_samples" >> {output.completion}
        echo "- Successful samples (phylogeny): $successful_count" >> {output.completion}
        echo "- Failed samples (excluded): $failed_count" >> {output.completion}
        echo "" >> {output.completion}
        echo "Analysis completed: $(date)" >> {output.completion}
        
        if [ $successful_count -eq 0 ]; then
            echo "WARNING: No samples were successful for phylogenetic analysis!" >> {output.completion}
            echo "WARNING: No samples were successful for phylogenetic analysis!" >> {log}
        else
            echo "SUCCESS: Phylogenetic analysis completed for $successful_count samples" >> {log}
        fi
        """
