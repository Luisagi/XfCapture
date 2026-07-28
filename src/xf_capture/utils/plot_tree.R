# ==============================================================
# plot_tree.R
# Generates a phylogenetic tree (ggtree) colored by subspecies
# of Xylella fastidiosa, with automatic subspecies inference
# for samples without metadata (e.g., sample-A) via phylogenetic proximity.
# ==============================================================

suppressMessages({
  library(ape)       # Phylogenetic tree handling
  library(phytools)  # Midpoint rooting
  library(ggtree)    # Advanced tree visualization
  library(treeio)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------
# 1. Parse command-line arguments
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript plot_tree.R <tree_file> <tree_name> <sample_name> <meta_data> <n_genes>")
}

tree_file   <- args[1]   # Input tree file path
tree_name   <- args[2]   # Output file identifier
sample_name <- args[3]   # Sample(s) to highlight in the tree ("|"-joined regex for multiple)
meta_data   <- args[4]   # Reference metadata CSV path
n_genes     <- args[5]   # Number of genes in the concatenated alignment

# ------------------------------
# 2. Load tree and metadata
# ------------------------------
tree  <- read.tree(tree_file)
mdata <- read.table(meta_data, sep = "\t", header = TRUE)

# ------------------------------
# 3. Midpoint root the tree
# ------------------------------
tree_mid <- midpoint.root(tree)

# ------------------------------
# 4. Ladderize (reorder clades for visual readability)
# ------------------------------
tree_mid <- ladderize(tree_mid, right = FALSE)

# ------------------------------
# 5. Build full labels and join metadata to tree
# ------------------------------
mdata$full_label <- paste0("ST-", mdata$ST, "_", mdata$Country_code, "_", mdata$Strain)

tip_data <- merge(
  data.frame(
    Strain    = tree_mid$tip.label,
    highlight = grepl(sample_name, tree_mid$tip.label)
  ),
  mdata, by = "Strain", all.x = TRUE
)

# If no full_label (tip without metadata match), use Strain as-is
tip_data$full_label <- ifelse(is.na(tip_data$full_label), tip_data$Strain, tip_data$full_label)

# ------------------------------
# 6. Infer subspp for tips without metadata (NA)
# ------------------------------
# Logic: For each tip without known subspecies, traverse up the tree
# level by level until finding the first clade with neighbors of
# known subspecies, and assign the most frequent one among them.
# "ambiguous = TRUE" indicates the neighboring clade mixes more than
# one subspecies, making the inference less reliable.
infer_subspp <- function(tip, tree, known_subspp) {
  tip_node     <- which(tree$tip.label == tip)
  current_node <- tip_node
  
  repeat {
    parent <- tree$edge[tree$edge[, 2] == current_node, 1]
    if (length(parent) == 0) {
      return(data.frame(inferred = NA, n_neighbors = 0, ambiguous = NA))
    }
    
    clade_tips <- extract.clade(tree, parent)$tip.label
    clade_tips <- setdiff(clade_tips, tip)
    
    neighbor_subspp <- known_subspp[intersect(clade_tips, names(known_subspp))]
    
    if (length(neighbor_subspp) > 0) {
      tbl       <- sort(table(neighbor_subspp), decreasing = TRUE)
      top       <- names(tbl)[1]
      ambiguous <- length(tbl) > 1
      return(data.frame(inferred = top, n_neighbors = length(neighbor_subspp), ambiguous = ambiguous))
    }
    
    current_node <- parent
  }
}

known_subspp <- setNames(tip_data$Subspp, tip_data$Strain)
known_subspp <- known_subspp[!is.na(known_subspp)]

na_tips <- tip_data$Strain[is.na(tip_data$Subspp)]

inferred_df <- lapply(na_tips, function(t) {
  data.frame(Strain = t, infer_subspp(t, tree_mid, known_subspp))
}) %>% bind_rows()

print(inferred_df)  # check n_neighbors / ambiguous before trusting inference

tip_data$Subspp[match(inferred_df$Strain, tip_data$Strain)] <- inferred_df$inferred

# ------------------------------
# 7. Detect maximal monophyletic subclades by subspecies
# ------------------------------
# Important: Several subspecies (e.g., fastidiosa) are NOT monophyletic
# in this tree -> a single MRCA() per group would produce a giant node
# encompassing other subspecies too. This function avoids that by returning
# AS MANY nodes as separate monophyletic blocks each group has.
get_max_clean_clades <- function(tips, tree) {
  tips <- intersect(tips, tree$tip.label)
  if (length(tips) == 0) return(data.frame(node = integer(0)))
  
  n_tips    <- length(tree$tip.label)
  remaining <- tips
  result    <- list()
  
  while (length(remaining) > 0) {
    tip           <- remaining[1]
    current_node  <- which(tree$tip.label == tip)
    best_node     <- current_node
    
    repeat {
      parent <- tree$edge[tree$edge[, 2] == current_node, 1]
      if (length(parent) == 0) break
      clade_tips <- extract.clade(tree, parent)$tip.label
      if (all(clade_tips %in% tips)) {
        best_node    <- parent
        current_node <- parent
      } else {
        break
      }
    }
    
    # If best_node is still a tip (size-1 block),
    # extract.clade() doesn't apply -> the tip itself is covered
    covered_tips <- if (best_node <= n_tips) tip else extract.clade(tree, best_node)$tip.label
    
    result[[length(result) + 1]] <- best_node
    remaining <- setdiff(remaining, covered_tips)
  }
  
  data.frame(node = unlist(result))
}

groupInfo <- split(tip_data$Strain, tip_data$Subspp)

clade_nodes_df <- lapply(names(groupInfo), function(g) {
  nodes_df <- get_max_clean_clades(groupInfo[[g]], tree_mid)
  if (nrow(nodes_df) == 0) return(NULL)
  data.frame(node = nodes_df$node, type = g)
}) %>% bind_rows() %>% filter(!is.na(node))

print(clade_nodes_df)  # check how many blocks per subspecies were found

clades <- setNames(as.integer(clade_nodes_df$node), clade_nodes_df$type)

# ------------------------------
# 8. Color branches by subspecies using groupClade()
# ------------------------------
p <- ggtree(tree_mid) %<+% tip_data
p <- groupClade(p, clades, group_name = "subtree") + aes(color = subtree)

# Note: Text labels are ALWAYS black (fixed color, not aes),
# to avoid inheriting \"subtree\" color -> prevents scale conflict
# that previously existed between aes(color = highlight) and scale_color_manual().
# sample-A (or other highlighted sample via `highlight`) is marked in bold.
p <- p +
  geom_tiplab(
    aes(label = full_label, fontface = ifelse(highlight, "bold", "plain")),
    color = "black", size = 2.5,
    align = TRUE, linetype = "dotted") +
  geom_tippoint(aes(subset = highlight), color = "red") +
  scale_color_manual(
    values = c(
      "fastidiosa"    = "#0072B2",
      "sandyi"        = "#6A3D9A",
      "morus"         = "#CC79A7",
      "pauca"         = "#E69F00",
      "multiplex"     = "#009E73"
    ),
    name = "Subspecies"
  ) +
  theme_tree2() +
  # theme(legend.position = "inside",
  #       legend.position.inside = c(0.05, 1),
  #       legend.justification = c(0, 1),
  #       legend.title = element_text(size = 15),
  #       legend.text = element_text(face = "italic", size = 14)) +
  theme(legend.position = "right",
        legend.justification = c(0, 1),
        legend.title = element_text(size = 15),
        legend.text = element_text(face = "italic", size = 14)) +
  hexpand(.10, direction = 1) +
  vexpand(.015, direction = -1) +

  labs(
    title = "Maximum-likelihood phylogenetic tree inferred using the best-fit partitioned model (MFP+MERGE)",
    subtitle = paste0(n_genes, " genes used")
  ) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

# highlight every tip matching sample_name (may be several, joined by "|")
sample_nodes <- grep(sample_name, tree_mid$tip.label)

for (node in sample_nodes) {
  p <- p + geom_hilight(node = node, fill = "red", alpha = .15, extend = 0.01)
}

# ------------------------------
# 9. Save as PDF
# ------------------------------
out_file <- paste0(tree_name, "_tree.pdf")

pdf(out_file, width = 12, height = 14.5)
print(p)
dev.off()

message("Tree plot saved as: ", out_file)