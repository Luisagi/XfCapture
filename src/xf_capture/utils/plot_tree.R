# Load required libraries
suppressMessages({
  library(ape)       # For phylogenetic tree handling
  library(phytools)  # For midpoint rooting
  library(ggtree)    # For advanced tree visualization
})

# ------------------------------
# 1. Read command-line arguments
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_tree.R <tree_file> <tree_name> <sample_name>")
}

tree_file   <- args[1]   # Path to the input tree file
tree_name   <- args[2]   # Name identifier for output file
sample_name <- args[3]   # Sample name to highlight in the tree


# ------------------------------
# 2. Load the phylogenetic tree
# ------------------------------
# The script expects a Newick-formatted file,
# but it also works with Nexus if you change to read.nexus()
tree <- read.tree(tree_file)

# ------------------------------
# 3. Root tree at the midpoint
# ------------------------------
tree_mid <- midpoint.root(tree)

# ------------------------------
# 4. Ladderize tree (reorder clades)
# ------------------------------
# right = TRUE  -> ascending order
# right = FALSE -> descending order
tree_mid <- ladderize(tree_mid, right = FALSE)

# ------------------------------
# 5. Identify tip labels that contain "sample"
# ------------------------------
tip_data <- data.frame(
  label = tree_mid$tip.label,
  highlight = grepl(sample_name, tree_mid$tip.label)  # TRUE if "sample" is in the name
)

# ------------------------------
# 6. Plot tree with ggtree
# ------------------------------
p <- ggtree(tree_mid) %<+% tip_data +
  geom_tiplab(aes(color = highlight), size = 2) +   # tip labels, highlight matches in red
  scale_color_manual(values = c("black", "red")) +  # black for normal tips, red for matches
  theme(legend.position = "none")                   # remove legend

# ------------------------------
# 7. Save plot as PNG
# ------------------------------
out_file <- paste0(tree_name, "_tree.png")

png(out_file, width = 13, height = 14, units = "in", res = 600)
print(p)
dev.off()

message("Tree plot saved as: ", out_file)
