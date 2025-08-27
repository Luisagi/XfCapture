################################################################################
##
##  Script to create the summary table of counts for all analyzed samples.
##
################################################################################

##
##  Uso: Rscript reconstruction_summary.R input_dir
##

# Load required libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(stringr)
})

# Function to process reconstruction stats files
process_recon_stats <- function(stats_dir) {
  
  # Find all *_recon_stats.tsv files
  recon_files <- list.files(stats_dir, pattern = "*_recon_stats.tsv$", full.names = TRUE)
  
  if (length(recon_files) == 0) {
    stop("No reconstruction stats files found in ", stats_dir)
  }
  
  # Initialize lists to store data
  all_recon_data <- list()
  all_length_data <- list()
  
  for (file in recon_files) {
    # Extract sample name from filename
    sample_name <- basename(file) %>% str_remove("_recon_stats.tsv$")
    
    # Read the file
    tryCatch({
      data <- read_tsv(file, col_types = cols(
        `#Refsequence` = col_character(),
        Sequence_length = col_double(),
        Reconstructed_bases = col_double(),
        Reconstruction_percent = col_double()
      ), show_col_types = FALSE)
      
      # Extract reconstruction percentages
      recon_data <- data %>%
        select(target_gene = `#Refsequence`, !!sample_name := Reconstruction_percent)
      
      # Extract sequence lengths (only need once)
      length_data <- data %>%
        select(target_gene = `#Refsequence`, sequence_length = Sequence_length)
      
      all_recon_data[[sample_name]] <- recon_data
      all_length_data[[sample_name]] <- length_data
      
    }, error = function(e) {
      cat("Warning: Could not read file", file, ":", e$message, "\n")
    })
  }
  
  # Merge all reconstruction data
  recon_summary <- Reduce(function(x, y) {
    full_join(x, y, by = "target_gene")
  }, all_recon_data)
  
  # Get sequence lengths (from first file)
  lengths <- all_length_data[[1]] %>%
    select(target_gene, sequence_length)
  
  # Add sequence lengths to reconstruction summary
  recon_summary <- left_join(lengths, recon_summary, by = "target_gene")
  
  return(recon_summary)
  
}

# Function to process mapped reads stats files
process_mapped_stats <- function(stats_dir) {
  
  # Find all *_mapped_stats.tsv files
  mapped_files <- list.files(stats_dir, pattern = "*_mapped_stats.tsv$", full.names = TRUE)
  
  if (length(mapped_files) == 0) {
    stop("No mapped stats files found in ", stats_dir)
  }
  
  # Initialize list to store data
  all_mapped_data <- list()
  all_length_data <- list()
  
  for (file in mapped_files) {
    # Extract sample name from filename
    sample_name <- basename(file) %>% str_remove("_mapped_stats.tsv$")
    
    # Read the file
    tryCatch({
      # Assuming mapped stats has columns: #Refsequence, Sequence_length, Mapped_reads
      # Adjust column names based on actual file structure
      data <- read_tsv(file, show_col_types = FALSE)
      
      # remove last row 
      data <- data[-nrow(data), ]
      
      # Check column names and adjust accordingly
      col_names <- names(data)
      
      if (ncol(data) >= 3) {
        names(data)[1:3] <- c("target_gene", "sequence_length", "mapped_reads")
        
        # Extract mapped reads
        mapped_data <- data %>%
          select(target_gene, !!sample_name := mapped_reads)
        
        # Extract sequence lengths
        length_data <- data %>%
          select(target_gene, sequence_length)
        
        all_mapped_data[[sample_name]] <- mapped_data
        all_length_data[[sample_name]] <- length_data
      }
      
    }, error = function(e) {
      cat("Warning: Could not read file", file, ":", e$message, "\n")
    })
  }
  
  if (length(all_mapped_data) == 0) {
    return(NULL)
  }
  
  # Merge all mapped data
  mapped_summary <- Reduce(function(x, y) {
    full_join(x, y, by = "target_gene")
  }, all_mapped_data)
  
  # Get sequence lengths (from first file)
  lengths <- all_length_data[[1]] %>%
    select(target_gene, sequence_length)
  
  # Add sequence lengths to mapped summary
  mapped_summary <- left_join(lengths, mapped_summary, by = "target_gene")
  
  return(mapped_summary)
}


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

stats_dir <- args[1]
output_file <- ifelse(length(args) >= 2, args[2], "reconstruction_summary.xlsx")

# Check if directory exists
if (!dir.exists(stats_dir)) {
  stop("Directory does not exist: ", stats_dir)
}

# Process reconstruction stats
recon_data <- process_recon_stats(stats_dir)

# Process mapped reads stats
mapped_data <- process_mapped_stats(stats_dir)

# Create Excel workbook
wb <- createWorkbook()

# Add reconstruction percentages and mapped sheets
addWorksheet(wb, "Reconstruction_Percentage")
writeData(wb, "Reconstruction_Percentage", recon_data)

addWorksheet(wb, "Mapped_Reads")
writeData(wb, "Mapped_Reads", mapped_data)

# Save workbook
saveWorkbook(wb, output_file, overwrite = TRUE)
