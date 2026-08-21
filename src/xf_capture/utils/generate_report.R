################################################################################
##
##  Script to create the summary table of counts for all analyzed samples.
##
################################################################################

##
##  Uso: Rscript generate_pipeline_report.R qc_file gene_counts mlst_file stats_dir [output_file]
##

# Load required libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(stringr)
})

# Function to process mapped reads stats files
# This function reads all *_mapped_stats.tsv files and creates summary tables
# for mapped reads and coverage percentages across samples
process_stats <- function(stats_dir) {
  
  # Find all *_mapped_stats.tsv files in the stats directory
  mapped_files <- list.files(stats_dir, pattern = "*_mapped_stats.tsv$", full.names = TRUE)
  
  if (length(mapped_files) == 0) {
    stop("No mapped stats files found in ", stats_dir)
  }
  
  # Initialize lists to store data from all files
  all_mapped_data <- list()
  all_coverage_data <- list()
  all_length_data <- list()
  
  for (file in mapped_files) {
    # Extract sample name from filename by removing the "_mapped_stats.tsv" suffix
    sample_name <- basename(file) %>% str_remove("_mapped_stats.tsv$")
    
    # Read the file with error handling
    tryCatch({
      # Read the TSV file containing mapping statistics
      data <- read_tsv(file, show_col_types = FALSE)
      
      # Rename columns to standard names for easier processing
      data <- 
        data %>%
        rename(
          target_gene = `#rname`,     # Reference gene name
          sequence_length  = endpos,   # Sequence length
          mapped_reads = numreads,     # Number of mapped reads
          coverage = coverage          # Coverage percentage
        )

      # Extract mapped reads data for this sample
      mapped_data <- data %>%
        select(target_gene, !!sample_name := mapped_reads)
      
      # Extract coverage data for this sample
      cover_data <- data %>%
        select(target_gene, !!sample_name := coverage)
        
      # Extract sequence lengths for this sample
      length_data <- data %>%
          select(target_gene, sequence_length)
        
      # Store the processed data in lists
      all_mapped_data[[sample_name]] <- mapped_data
      all_coverage_data[[sample_name]] <- cover_data
      all_length_data[[sample_name]] <- length_data

      
    }, error = function(e) {
      cat("Warning: Could not read file", file, ":", e$message, "\n")
    })
  }
  
  # If no files were successfully processed, return NULL
  if (length(all_mapped_data) == 0) {
    return(NULL)
  }

  # Merge all mapped data by target_gene using full_join to preserve all genes
  mapped_summary <- Reduce(function(x, y) {
    full_join(x, y, by = "target_gene")
  }, all_mapped_data)
  
  # Merge all coverage data by target_gene using full_join to preserve all genes
  cover_summary <- Reduce(function(x, y) {
    full_join(x, y, by = "target_gene")
  }, all_coverage_data)
  
  # Get sequence lengths from the first file (assuming they're consistent across files)
  lengths <- all_length_data[[1]] %>%
    select(target_gene, sequence_length)
  
  # Add sequence lengths to both summary tables
  mapped_summary <- left_join(lengths, mapped_summary, by = "target_gene")
  cover_summary <- left_join(lengths, cover_summary, by = "target_gene")
  
  
  return(list(mapped_summary = mapped_summary, cover_summary = cover_summary))
}


# Input parsing ----------------------------------------------------------------
# Parse command line arguments for input files and output file
args <- commandArgs(trailingOnly = TRUE)

# Expected arguments:
# 1. Preprocessing summary file (qc_data)
# 2. Gene counts file
# 3. MLST results file
# 4. Stats directory, OR the first mapped-stats file when the caller passes many files
# 5. Output Excel file name (optional, defaults to "reconstruction_summary.xlsx")

preprocessing_data <- args[1]
gene_counts <- args[2]
mlst_res <- args[3]

# Accept either:
#   Rscript script.R qc.tsv genes.tsv mlst.csv stats_dir output.xlsx
# or:
#   Rscript script.R qc.tsv genes.tsv mlst.csv file1.tsv file2.tsv ... output.xlsx
if (length(args) >= 4 && dir.exists(args[4])) {
  stats_dir <- args[4]
  output_file <- if (length(args) >= 5) args[5] else "pipeline_report.xlsx"
} else {
  stats_dir <- dirname(args[4])
  output_file <- args[length(args)]
}

# Read input data files
qc <- read.delim(preprocessing_data)
gcounts <- read.delim(gene_counts)
mlst <- read.delim(mlst_res, sep = ",")

# Combine all data into a single summary table
# Join preprocessing data with gene counts and MLST results
summary <- 
  qc %>%
  left_join(gcounts, by = "sample_id") %>%
  left_join(mlst, by = c("sample_id" = "FILE")) %>%
  rename(
    sequence_type = ST,   # Rename ST column to sequence_type for clarity
    alleles = ALLELES     # Rename ALLELES column to alleles for clarity
  ) %>%
  select(-STATUS, -SCORE, -SCHEME)  # Remove unnecessary columns

# Validate that stats directory exists
if (!dir.exists(stats_dir)) {
  stop("Directory does not exist: ", stats_dir)
}

# Process reconstruction stats to get mapped reads and coverage data
sorted_data <- process_stats(stats_dir)

# Create Excel workbook for output
wb <- createWorkbook()

# Add the main summary sheet
addWorksheet(wb, "General_summary")
writeData(wb, "General_summary", summary)

# Add sheets for detailed mapping statistics
addWorksheet(wb, "Coverage_Percentage")
writeData(wb, "Coverage_Percentage", sorted_data$cover_summary)

addWorksheet(wb, "Mapped_Reads")
writeData(wb, "Mapped_Reads", sorted_data$mapped_summary)

# Save the workbook to the specified output file
saveWorkbook(wb, output_file, overwrite = TRUE)
