import os
import sys
import re

def get_samples(input_dir, valid_extensions, not_allowed_chars):
    """
    Get samples from the input directory and validate them.
    Supports multiple naming conventions:
    - CASAVA format: sample_S1_L001_R1_001.fastq.gz, sample_S1_L001_R2_001.fastq.gz
    - Standard format: sample_R1.fastq.gz, sample_R2.fastq.gz
    - Simple format: sample_1.fastq.gz, sample_2.fastq.gz
    
    Args:
        input_dir (str): Path to the input directory
        valid_extensions (tuple): Valid file extensions
        not_allowed_chars (set): Set of not allowed characters in sample names
    
    Returns:
        dict: Dictionary with sample names as keys and [R1, R2] paths as values
    """
    print(f"Looking for fastq files in directory: {input_dir}")
    
    if not os.path.exists(input_dir):
        print(f"Directory {input_dir} does not exist!")
        sys.exit(0)
    
    files = os.listdir(input_dir)
    
    # Filter files with valid extensions
    fastq_files = [f for f in files if f.endswith(valid_extensions)]
    
    if not fastq_files:
        print(f"No FASTQ files found with valid extensions: {valid_extensions}")
        sys.exit(1)
    
    print(f"Found {len(fastq_files)} FASTQ files")
    
    # Create a dictionary to hold sample names and their corresponding paths
    sample_dict = {}

    for f in fastq_files:
        sample_name = None
        read_type = None
        
        # Pattern 1: CASAVA format with sample number (sample_S1_L001_R1_001.fastq.gz)
        casava_match = re.match(r"(.+)_S\d+_L\d+_R([12])_\d+", f)
        if casava_match:
            sample_name = casava_match.group(1)
            read_type = casava_match.group(2)
        
        # Pattern 2: CASAVA format without sample number (sample_L001_R1_001.fastq.gz)
        elif re.match(r".+_L\d+_R[12]_\d+", f):
            casava_simple_match = re.match(r"(.+)_L\d+_R([12])_\d+", f)
            if casava_simple_match:
                sample_name = casava_simple_match.group(1)
                read_type = casava_simple_match.group(2)
        
        # Pattern 3: Standard format (sample_R1.fastq.gz, sample_R2.fastq.gz)
        elif "_R1" in f or "_R2" in f:
            standard_match = re.match(r"(.+)_R([12])", f)
            if standard_match:
                sample_name = standard_match.group(1)
                read_type = standard_match.group(2)
        
        # Pattern 4: Simple format (sample_1.fastq.gz, sample_2.fastq.gz)
        elif "_1." in f or "_2." in f:
            simple_match = re.match(r"(.+)_([12])\.", f)
            if simple_match:
                sample_name = simple_match.group(1)
                read_type = simple_match.group(2)
        
        # If we found a valid pattern, add to dictionary
        if sample_name and read_type:
            path = os.path.join(input_dir, f)
            
            if sample_name not in sample_dict:
                sample_dict[sample_name] = ["", ""]
            
            # Assign to R1 (index 0) or R2 (index 1)
            if read_type == "1":
                sample_dict[sample_name][0] = path
            elif read_type == "2":
                sample_dict[sample_name][1] = path

    if not sample_dict:
        print("ERROR: No valid sample pairs found!")
        print("Supported formats:")
        print("  - CASAVA with sample: sample_S1_L001_R1_001.fastq.gz")
        print("  - CASAVA simple: sample_L001_R1_001.fastq.gz")
        print("  - Standard: sample_R1.fastq.gz")
        print("  - Simple: sample_1.fastq.gz")
        sys.exit(1)

    # Check if all samples have both R1 and R2 files
    incomplete = {s: p for s, p in sample_dict.items() if not all(p)}
    if incomplete:
        print(f"WARNING: Found samples with missing pairs:")
        for sample, paths in incomplete.items():
            print(f"  - {sample}: R1: {paths[0] or 'MISSING'}, R2: {paths[1] or 'MISSING'}")
        
        # Remove incomplete samples
        for sample in incomplete:
            del sample_dict[sample]
        
        if not sample_dict:
            print("ERROR: No complete sample pairs found!")
            sys.exit(1)

    # Extract sample names from the dictionary keys
    samples = sorted(sample_dict.keys())
    print(f"Found {len(samples)} complete sample pairs:")
    for sample in samples:
        print(f"  - {sample}")

    # Check for non-allowed characters in sample names
    invalid_samples = []
    for sample in samples:
        invalid_chars = [char for char in sample if char in not_allowed_chars]
        if invalid_chars:
            invalid_samples.append((sample, invalid_chars))
            print(f"WARNING: Sample '{sample}' contains non-allowed characters: {invalid_chars}")
    
    if invalid_samples:
        print(f"\nERROR: Found {len(invalid_samples)} sample(s) with non-allowed characters:")
        for sample, chars in invalid_samples:
            print(f"  - Sample '{sample}' has characters: {chars}")
        print(f"Non-allowed characters are: {sorted(not_allowed_chars)}")
        sys.exit(1)
    else:
        print(f"✓ All sample names are valid (no forbidden characters found)")
    
    return sample_dict

# def get_samples(input_dir, valid_extensions, not_allowed_chars):
#     """
#     Get samples from the input directory and validate them.
    
#     Args:
#         input_dir (str): Path to the input directory
#         valid_extensions (tuple): Valid file extensions
#         not_allowed_chars (set): Set of not allowed characters in sample names
    
#     Returns:
#         dict: Dictionary with sample names as keys and [R1, R2] paths as values
#     """
#     print(f"Looking for fastq files in directory: {input_dir}")
    
#     if not os.path.exists(input_dir):
#         print(f"Directory {input_dir} does not exist!")
#         sys.exit(0)
    
#     files = os.listdir(input_dir)
    
#     # Filter files with valid extensions (both R1 and R2)
#     fastq_files = [f for f in files if f.endswith(valid_extensions)]
    
#     # Create a dictionary to hold sample names and their corresponding paths
#     sample_dict = {}

#     for f in fastq_files:
#         if "_R1_001" in f or "_R2_001" in f:
#             match = re.match(r"([^_]+)_", f)
#             if match:
#                 sample_name = match.group(1)
#                 path = os.path.join(input_dir, f)
#                 if sample_name not in sample_dict:
#                     sample_dict[sample_name] = ["", ""]
#                 if "_R1_001" in f:
#                     sample_dict[sample_name][0] = path
#                 elif "_R2_001" in f:
#                     sample_dict[sample_name][1] = path

#     # Check if all samples have both R1 and R2 files
#     incomplete = {s: p for s, p in sample_dict.items() if not all(p)}
#     if incomplete:
#         print(f"WARNING: Found samples with missing pairs:")
#         for sample, paths in incomplete.items():
#             print(f"  - {sample}: R1: {paths[0] or 'MISSING'}, R2: {paths[1] or 'MISSING'}")
#         sys.exit(1)

#     # Extract sample names from the dictionary keys
#     samples = sorted(sample_dict.keys())

#     # Check for non-allowed characters in sample names
#     invalid_samples = []
#     for sample in samples:
#         invalid_chars = [char for char in sample if char in not_allowed_chars]
#         if invalid_chars:
#             invalid_samples.append((sample, invalid_chars))
#             print(f"WARNING: Sample '{sample}' contains non-allowed characters: {invalid_chars}")
    
#     if invalid_samples:
#         print(f"\nERROR: Found {len(invalid_samples)} sample(s) with non-allowed characters:")
#         for sample, chars in invalid_samples:
#             print(f"  - Sample '{sample}' has characters: {chars}")
#         print(f"Non-allowed characters are: {sorted(not_allowed_chars)}")
#         sys.exit(1)
#     else:
#         print(f"✓ All sample names are valid (no forbidden characters found)")
    
#     return sample_dict


def get_r1(wildcards, samples_dict):
    """Get R1 file path for a given sample"""
    return samples_dict[wildcards.sample][0]


def get_r2(wildcards, samples_dict):
    """Get R2 file path for a given sample"""
    return samples_dict[wildcards.sample][1]


def get_fasta_files(directory):
    """
    Detect FASTA files in the input directory.
    
    Args:
        directory (str): Path to the directory to search
    
    Returns:
        dict: Dictionary with basename as keys and [full_path] as values
    """
    fasta_extensions = ('.fa', '.fasta', '.fna')
    fasta_dict = {}

    if not os.path.exists(directory):
        print(f"Directory {directory} does not exist!")
        return {}
    
    for file in os.listdir(directory):
        if file.lower().endswith(fasta_extensions):
            full_path = os.path.join(directory, file)
            basename = os.path.splitext(os.path.basename(file))[0]
            fasta_dict[basename] = [full_path]
    
    return fasta_dict

def print_fasta_headers(fasta_dict):
    """
    Print the headers of the FASTA files (placeholder function)
    
    Args:
        fasta_dict (dict): Dictionary of FASTA files
    """
    for fasta, paths in fasta_dict.items():
        print(f"Headers in {fasta}:")
        # This would need the get_fasta_headers function implementation
        # headers = get_fasta_headers(paths[0])
        # for header in headers:
        #     print(f"  - {header}")
        print()  # New line for better readability

def get_gene_names(sample):
    """
    Get list of gene names for a given sample from headers_ids.txt
    
    Args:
        sample (str): Sample name
    
    Returns:
        list: List of gene names
    """
    headers_file = f"05.phylogenetic_trees/{sample}/headers_ids.txt"
    gene_names = []
    
    if os.path.exists(headers_file):
        try:
            with open(headers_file, 'r') as f:
                gene_names = [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"Error reading {headers_file}: {e}")
    
    return gene_names

def get_alignment_outputs(wildcards):
    """
    Get alignment output files for a sample based on available genes
    
    Args:
        wildcards: Snakemake wildcards object
    
    Returns:
        list: List of alignment output file paths
    """
    genes = get_gene_names(wildcards.sample)
    if not genes:
        return []
    return [f"05.phylogenetic_trees/{wildcards.sample}/alignments/{gene}_trimmed.fasta" for gene in genes]


def get_phylo_outputs(wildcards):
    """
    Get phylogeny output files for a sample based on available genes
    
    Args:
        wildcards: Snakemake wildcards object
    
    Returns:
        list: List of phylogeny output file paths
    """
    genes = get_gene_names(wildcards.sample)
    if not genes:
        return []
    return [f"05.phylogenetic_trees/{wildcards.sample}/trees/{gene}.treefile" for gene in genes]

def get_successful_samples():
    """
    Get list of samples that successfully reconstructed sequences
    by reading checkpoint files.
    
    Returns:
        list: List of sample names with successful reconstruction
    """
    import os
    successful_samples = []
    
    # Directory where checkpoint files are stored
    checkpoint_dir = "03.probes_reconstruction/checkpoint"
    
    if not os.path.exists(checkpoint_dir):
        return successful_samples
    
    # Look for all checkpoint files
    for filename in os.listdir(checkpoint_dir):
        if filename.endswith("_reconstruction_check.txt"):
            # Extract sample name from filename
            sample_name = filename.replace("_reconstruction_check.txt", "")
            check_file = os.path.join(checkpoint_dir, filename)
            
            try:
                with open(check_file, 'r') as f:
                    content = f.read()
                    if "SUCCESS:" in content:
                        successful_samples.append(sample_name)
            except Exception as e:
                print(f"Error reading {check_file}: {e}")
    
    return sorted(successful_samples)