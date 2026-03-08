# CHIP Mutation Analysis Script
# Generates oncoplot and co-mutation analysis from smMIPTools output

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to called_mutations.txt input file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = ".",
              help = "Output directory [default: current directory]", metavar = "DIR"),
  make_option(c("-c", "--code"), type = "character", default = NULL,
              help = "Path to CHIP R code directory (for sourcing dependencies)", metavar = "DIR"),
  make_option(c("-v", "--vaf_threshold"), type = "numeric", default = 0.01,
              help = "VAF threshold for filtering [default: 0.01]", metavar = "NUM"),
  make_option(c("-p", "--pval_threshold"), type = "numeric", default = 0.05,
              help = "P-value threshold for filtering [default: 0.05]", metavar = "NUM"),
  make_option(c("-d", "--min_depth"), type = "integer", default = 0,
              help = "Minimum read depth to retain a variant. Variants below this threshold are flagged as LOW_DEPTH [default: 0 (disabled)]", metavar = "INT")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "CHIP Population Analysis - generates oncoplots, VAF heatmaps, and co-mutation analysis")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) {
  stop("Error: --input argument is required. Use -h for help.")
}

if (!file.exists(opt$input)) {
  stop(paste("Error: Input file not found:", opt$input))
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

cat("CHIP Analysis Parameters:\n")
cat("  Input file:", opt$input, "\n")
cat("  Output directory:", opt$output, "\n")
cat("  VAF threshold:", opt$vaf_threshold, "\n")
cat("  P-value threshold:", opt$pval_threshold, "\n")
cat("  Min depth filter:", ifelse(opt$min_depth > 0, opt$min_depth, "disabled"), "\n\n")

# ============================================================================
# SETUP AND DEPENDENCIES
# ============================================================================

# Install packages if needed (uncomment if required)
# install.packages(c("tidyverse", "ComplexHeatmap", "circlize", "RColorBrewer"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# ============================================================================
# DATA LOADING AND FILTERING
# ============================================================================

# Load data from input file
data <- read.delim(opt$input, stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "\"")

cat("Original data dimensions:", nrow(data), "rows x", ncol(data), "columns\n")

# Apply filtering criteria:
# 1. Flag mutations with quality flags (retained for review, not removed)
# 2. Allele frequency > VAF threshold (from command line)
# 3. P-value < P-value threshold (from command line)
#
# NOTE: Flagged mutations are kept but marked in 'flag_status' column.
# Removing all flagged mutations risks discarding true positives
# (e.g., variants near homopolymers that are real but flagged by
# conservative quality filters). Downstream analysis should filter
# on flag_status as appropriate for the specific use case.

data$flag_status <- ifelse(
  is.na(data$flags) | trimws(data$flags) == "",
  "PASS", "FLAGGED"
)

filtered_data <- data %>%
  filter(
    # Allele frequency > threshold
    allele.frequency > opt$vaf_threshold,
    # P-value < threshold
    `P-value` < opt$pval_threshold
  )

cat("After filtering:", nrow(filtered_data), "rows\n")

# Flag variants below minimum depth threshold (if enabled)
if (opt$min_depth > 0 && "coverage" %in% names(filtered_data)) {
  low_depth_idx <- which(filtered_data$coverage < opt$min_depth)
  if (length(low_depth_idx) > 0) {
    filtered_data$flag_status[low_depth_idx] <- paste0(filtered_data$flag_status[low_depth_idx], ";LOW_DEPTH")
    cat("Flagged", length(low_depth_idx), "variants with depth <", opt$min_depth, "as LOW_DEPTH\n")
  }
} else if (opt$min_depth > 0 && "total.coverage" %in% names(filtered_data)) {
  low_depth_idx <- which(filtered_data$total.coverage < opt$min_depth)
  if (length(low_depth_idx) > 0) {
    filtered_data$flag_status[low_depth_idx] <- paste0(filtered_data$flag_status[low_depth_idx], ";LOW_DEPTH")
    cat("Flagged", length(low_depth_idx), "variants with depth <", opt$min_depth, "as LOW_DEPTH\n")
  }
} else if (opt$min_depth > 0) {
  cat("WARNING: --min_depth specified but no coverage/total.coverage column found; skipping depth filter\n")
}

# Check if any mutations passed filtering
if (nrow(filtered_data) == 0) {
  cat("\nWARNING: No mutations passed filtering criteria.\n")
  cat("Consider adjusting --vaf_threshold or --pval_threshold\n")
  cat("No output files generated.\n")
  quit(status = 0)
}

# ============================================================================
# SAMPLE ID PROCESSING (MERGE REPLICATES)
# ============================================================================

# Extract base sample ID by removing replicate suffixes (A/B)
# e.g., "D79A,D79B" -> "D79"
filtered_data <- filtered_data %>%
  mutate(
    # First, take the first sample ID if comma-separated
    base_sample = str_extract(sample_ID, "^[^,]+"),
    # Then remove the trailing letter (A or B replicate indicator)
    base_sample = str_remove(base_sample, "[AB]$")
  )

cat("Unique samples after merging replicates:",
    n_distinct(filtered_data$base_sample), "\n")
cat("Unique genes:", n_distinct(filtered_data$gene), "\n")

# ============================================================================
# CREATE MUTATION MATRIX FOR ONCOPLOT
# ============================================================================

# Function to extract clean protein change from annotation
# e.g., "ENST00000539801:V617F; ENST00000381652:V617F" -> "V617F"
extract_protein_change <- function(protein_str) {
  if (is.na(protein_str) || protein_str == "NA" || protein_str == "") {
    return(NA)
  }
  # Take the first annotation and extract the change after the colon
  first_annot <- str_split(protein_str, ";")[[1]][1]
  if (str_detect(first_annot, ":")) {
    change <- str_extract(first_annot, ":[A-Z*]?[0-9]+[A-Z*]?$")
    if (!is.na(change)) {
      return(str_remove(change, "^:"))
    }
    # Try alternate pattern for simple annotations like "V430M"
    change <- str_extract(first_annot, "[A-Z][0-9]+[A-Z*]$")
    return(change)
  } else {
    # Simple format like "V430M" or "R834*"
    return(trimws(first_annot))
  }
}

# Add protein change annotations to filtered data
filtered_data <- filtered_data %>%
  rowwise() %>%
  mutate(
    protein_change = extract_protein_change(protein),
    # Create gene_variant label (e.g., "JAK2 V617F")
    gene_variant = if_else(
      !is.na(protein_change) & protein_change != "",
      paste0(gene, " ", protein_change),
      gene
    )
  ) %>%
  ungroup()

# Create a unique mutation identifier per sample-gene combination
# If same sample has multiple mutations in same gene, we'll classify by variant type
mutation_matrix_long <- filtered_data %>%
  select(base_sample, gene, variant_type, protein_change, alt) %>%
  distinct() %>%
  # Simplify variant types for visualization
  # Use alt column to distinguish insertions (+) from deletions (-)
  mutate(
    variant_class = case_when(
      str_detect(variant_type, "missense|nonsynonymous") ~ "Missense",
      str_detect(variant_type, "nonsense|stop_gained|stopgain") ~ "Nonsense",
      str_detect(variant_type, "frameshift") & alt == "+" ~ "Insertion",
      str_detect(variant_type, "frameshift") & alt == "-" ~ "Deletion",
      str_detect(variant_type, "frameshift") ~ "Indel",
      str_detect(variant_type, "splice") ~ "Splice_Site",
      str_detect(variant_type, "inframe") ~ "In_Frame",
      str_detect(variant_type, "^synonymous") ~ "Synonymous",  # Anchor to start to avoid matching nonsynonymous
      str_detect(variant_type, "intron") ~ "Intronic",
      str_detect(variant_type, "non_coding") ~ "Non_Coding",
      is.na(variant_type) & alt == "+" ~ "Insertion",
      is.na(variant_type) & alt == "-" ~ "Deletion",
      is.na(variant_type) ~ "Indel",
      TRUE ~ "Other"
    )
  )

# Create wide matrix for oncoplot (genes as rows, samples as columns)
# Value is the variant classification
mutation_matrix <- mutation_matrix_long %>%
  select(base_sample, gene, variant_class) %>%
  distinct() %>%
  # If multiple variant types per gene-sample, concatenate them
  group_by(base_sample, gene) %>%
  summarise(variant_class = paste(unique(variant_class), collapse = ";"),
            .groups = "drop") %>%
  pivot_wider(names_from = base_sample, values_from = variant_class,
              values_fill = "")

# Create protein annotation matrix (for displaying in cells)
protein_annotation_matrix <- mutation_matrix_long %>%
  select(base_sample, gene, protein_change) %>%
  filter(!is.na(protein_change)) %>%
  distinct() %>%
  group_by(base_sample, gene) %>%
  summarise(protein_label = paste(unique(protein_change), collapse = "\n"),
            .groups = "drop") %>%
  pivot_wider(names_from = base_sample, values_from = protein_label,
              values_fill = "")

# Create variant count matrix (number of distinct mutations per gene-sample)
variant_count_matrix <- filtered_data %>%
  select(base_sample, gene, chr, pos, ref, alt) %>%
  distinct() %>%
  group_by(base_sample, gene) %>%
  summarise(n_variants = n(), .groups = "drop") %>%
  pivot_wider(names_from = base_sample, values_from = n_variants,
              values_fill = 0)

# Convert to matrix format for ComplexHeatmap
genes <- mutation_matrix$gene
mat <- as.matrix(mutation_matrix[, -1])
rownames(mat) <- genes

# Convert protein annotation to matrix
protein_genes <- protein_annotation_matrix$gene
protein_mat <- as.matrix(protein_annotation_matrix[, -1])
rownames(protein_mat) <- protein_genes

# Ensure protein_mat has same dimensions as mat
# Fill missing combinations with empty strings
protein_mat_full <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                           dimnames = list(rownames(mat), colnames(mat)))
for (g in rownames(protein_mat)) {
  for (s in colnames(protein_mat)) {
    if (g %in% rownames(protein_mat_full) && s %in% colnames(protein_mat_full)) {
      protein_mat_full[g, s] <- protein_mat[g, s]
    }
  }
}

# Convert variant count to matrix
count_genes <- variant_count_matrix$gene
count_mat <- as.matrix(variant_count_matrix[, -1])
rownames(count_mat) <- count_genes

# Ensure count_mat has same dimensions as mat
count_mat_full <- matrix(0, nrow = nrow(mat), ncol = ncol(mat),
                         dimnames = list(rownames(mat), colnames(mat)))
for (g in rownames(count_mat)) {
  for (s in colnames(count_mat)) {
    if (g %in% rownames(count_mat_full) && s %in% colnames(count_mat_full)) {
      count_mat_full[g, s] <- count_mat[g, s]
    }
  }
}

cat("\nMutation matrix dimensions:", nrow(mat), "genes x", ncol(mat), "samples\n")

# ============================================================================
# CALCULATE MUTATION FREQUENCIES
# ============================================================================

# Calculate mutation frequency per gene (for ordering)
# Handle single-gene case where matrix becomes vector
if (nrow(mat) == 1) {
  gene_freq <- sum(mat != "") / ncol(mat) * 100
  names(gene_freq) <- rownames(mat)
  gene_order <- rownames(mat)
} else {
  gene_freq <- rowSums(mat != "") / ncol(mat) * 100
  gene_order <- names(sort(gene_freq, decreasing = TRUE))
}

# Filter to genes mutated in at least 1 sample and reorder
mat <- mat[gene_order, , drop = FALSE]
if (nrow(mat) == 1) {
  keep_rows <- sum(mat != "") > 0
} else {
  keep_rows <- rowSums(mat != "") > 0
}
mat <- mat[keep_rows, , drop = FALSE]

# Reorder annotation matrices to match
protein_mat_full <- protein_mat_full[rownames(mat), colnames(mat), drop = FALSE]
count_mat_full <- count_mat_full[rownames(mat), colnames(mat), drop = FALSE]

cat("Genes with at least 1 mutation:", nrow(mat), "\n")

# ============================================================================
# ONCOPLOT VISUALIZATION
# ============================================================================

# Define colors for variant types
variant_colors <- c(
  "Missense" = "#3498db",
  "Nonsense" = "#e74c3c",
  "Insertion" = "#9b59b6",
  "Deletion" = "#8e44ad",
  "Indel" = "#a569bd",
  "Splice_Site" = "#f39c12",
  "In_Frame" = "#1abc9c",
  "Synonymous" = "#95a5a6",
  "Intronic" = "#7f8c8d",
  "Non_Coding" = "#bdc3c7",
  "Other" = "#2c3e50"
)

# Define the alter_fun for oncoPrint (standard version without annotations)
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "#f0f0f0", col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Missense"], col = NA))
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Nonsense"], col = NA))
  },
  Insertion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Insertion"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Deletion"], col = NA))
  },
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Indel"], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Splice_Site"], col = NA))
  },
  In_Frame = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["In_Frame"], col = NA))
  },
  Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Synonymous"], col = NA))
  },
  Intronic = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Intronic"], col = NA))
  },
  Non_Coding = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Non_Coding"], col = NA))
  },
  Other = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.9,
              gp = gpar(fill = variant_colors["Other"], col = NA))
  }
)

# Get the variant types actually present in the data
present_variants <- unique(unlist(str_split(mat[mat != ""], ";")))
present_variants <- present_variants[present_variants %in% names(variant_colors)]

# Create legend
legend_colors <- variant_colors[present_variants]

# Generate oncoplot with protein annotations
pdf(file.path(opt$output, "CHIP_oncoplot.pdf"), width = 16, height = 12)

oncoplot <- oncoPrint(
  mat,
  name = "CHIP",
  alter_fun = alter_fun,
  col = variant_colors,
  column_title = "CHIP Mutations Oncoplot",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Alterations",
    at = present_variants,
    labels = present_variants
  ),
  row_names_gp = gpar(fontsize = 10),
  pct_gp = gpar(fontsize = 9),
  row_order = 1:nrow(mat),  # Keep frequency-based ordering
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  top_annotation = HeatmapAnnotation(
    `# Mutations` = anno_barplot(
      colSums(mat != ""),
      height = unit(2, "cm"),
      gp = gpar(fill = "#2c3e50")
    ),
    annotation_name_side = "left"
  ),
  right_annotation = rowAnnotation(
    `# Samples` = anno_barplot(
      rowSums(mat != ""),
      width = unit(2, "cm"),
      gp = gpar(fill = "#27ae60")
    )
  )
)

# Draw the oncoplot (no cell annotations)
draw(oncoplot)
dev.off()

# Also save as PNG for easy viewing
png(file.path(opt$output, "CHIP_oncoplot.png"), width = 1600, height = 1200, res = 100)
draw(oncoplot)
dev.off()

cat("\nOncoplot saved to:", file.path(opt$output, "CHIP_oncoplot.pdf"), "and CHIP_oncoplot.png\n")

# Save variant annotation summary
variant_summary <- filtered_data %>%
  select(base_sample, gene, protein_change, gene_variant, variant_type,
         SSCS.allele.frequency, allele.frequency) %>%
  distinct() %>%
  arrange(gene, base_sample)

write.csv(variant_summary, file.path(opt$output, "variant_annotation_summary.csv"), row.names = FALSE)
cat("Variant annotation summary saved to:", file.path(opt$output, "variant_annotation_summary.csv"), "\n")

# ============================================================================
# VAF (VARIANT ALLELE FREQUENCY) HEATMAP
# ============================================================================

# Create VAF matrix - for samples with multiple mutations in same gene, take max VAF
# Using SSCS.allele.frequency (Single-Strand Consensus Sequence) for more accurate VAF
vaf_matrix_long <- filtered_data %>%
  select(base_sample, gene, SSCS.allele.frequency) %>%
  filter(!is.na(SSCS.allele.frequency)) %>%  # Remove NA values
  group_by(base_sample, gene) %>%
  summarise(max_vaf = max(SSCS.allele.frequency, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(max_vaf))  # Remove any -Inf values

# Pivot to wide format (genes as rows, samples as columns)
vaf_matrix <- vaf_matrix_long %>%
  pivot_wider(names_from = base_sample, values_from = max_vaf, values_fill = 0)

# Convert to matrix
vaf_genes <- vaf_matrix$gene
vaf_mat <- as.matrix(vaf_matrix[, -1])
rownames(vaf_mat) <- vaf_genes

# Replace any non-finite values with 0
vaf_mat[!is.finite(vaf_mat)] <- 0

# Order genes by frequency (same as oncoplot)
vaf_mat <- vaf_mat[intersect(gene_order, rownames(vaf_mat)), , drop = FALSE]

# Order samples to match oncoplot if possible
if (all(colnames(mat) %in% colnames(vaf_mat))) {
  vaf_mat <- vaf_mat[, colnames(mat), drop = FALSE]
}

cat("\nVAF matrix dimensions:", nrow(vaf_mat), "genes x", ncol(vaf_mat), "samples\n")
cat("VAF range:", round(min(vaf_mat[vaf_mat > 0]), 4), "-", round(max(vaf_mat), 4), "\n")

# Define color scale for VAF (white = 0, gradient to dark red for high VAF)
vaf_col_fun <- colorRamp2(
  c(0, 0.05, 0.1, 0.2, 0.5),
  c("white", "#fee5d9", "#fcae91", "#fb6a4a", "#a50f15")
)

# Create VAF heatmap
pdf(file.path(opt$output, "VAF_heatmap.pdf"), width = 12, height = 8)

# Create annotation for row (gene) mutation counts
gene_counts <- rowSums(vaf_mat > 0)

# Create annotation for column (sample) total mutations
sample_mutation_counts <- colSums(vaf_mat > 0)

vaf_heatmap <- Heatmap(
  vaf_mat,
  name = "VAF",
  col = vaf_col_fun,
  column_title = "Variant Allele Frequency (VAF) by Gene and Sample",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = "Genes",
  row_title_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,  # Keep gene order from oncoplot
  cluster_columns = TRUE,  # Cluster samples by VAF pattern
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Add VAF value text for non-zero cells
    if (vaf_mat[i, j] > 0) {
      grid.text(sprintf("%.2f", vaf_mat[i, j]), x, y,
                gp = gpar(fontsize = 7, col = ifelse(vaf_mat[i, j] > 0.15, "white", "black")))
    }
  },
  top_annotation = HeatmapAnnotation(
    `# Mutations` = anno_barplot(
      sample_mutation_counts,
      height = unit(1.5, "cm"),
      gp = gpar(fill = "#2c3e50")
    ),
    annotation_name_side = "left"
  ),
  right_annotation = rowAnnotation(
    `# Samples` = anno_barplot(
      gene_counts,
      width = unit(1.5, "cm"),
      gp = gpar(fill = "#27ae60")
    )
  ),
  heatmap_legend_param = list(
    title = "VAF",
    at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5+")
  )
)

draw(vaf_heatmap)
dev.off()

cat("VAF heatmap saved to:", file.path(opt$output, "VAF_heatmap.pdf"), "\n")

# Also save as PNG for easy viewing
png(file.path(opt$output, "VAF_heatmap.png"), width = 1200, height = 800, res = 100)
draw(vaf_heatmap)
dev.off()
cat("VAF heatmap PNG saved to:", file.path(opt$output, "VAF_heatmap.png"), "\n")

# Also create a lollipop/dot plot version showing VAF distribution per gene
pdf(file.path(opt$output, "VAF_dotplot.pdf"), width = 10, height = 8)

vaf_plot_data <- filtered_data %>%
  mutate(gene = factor(gene, levels = rev(gene_order)))

p <- ggplot(vaf_plot_data, aes(x = SSCS.allele.frequency, y = gene)) +
  geom_jitter(aes(color = base_sample), size = 3, alpha = 0.8, height = 0.2) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0.02, linetype = "dotted", color = "gray70", linewidth = 0.5) +
  scale_x_continuous(
    limits = c(0, max(vaf_plot_data$SSCS.allele.frequency, na.rm = TRUE) * 1.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_brewer(palette = "Set3", name = "Sample") +
  labs(
    title = "SSCS Variant Allele Frequency Distribution by Gene",
    subtitle = "Each point represents a mutation; dashed line = 10% VAF, dotted line = 2% VAF",
    x = "SSCS Variant Allele Frequency (VAF)",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

print(p)
dev.off()

# Also save as PNG
png(file.path(opt$output, "VAF_dotplot.png"), width = 1000, height = 800, res = 100)
print(p)
dev.off()

cat("VAF dot plot saved to:", file.path(opt$output, "VAF_dotplot.pdf"), "and VAF_dotplot.png\n")

# ============================================================================
# CO-MUTATION ANALYSIS (PAIRWISE FISHER'S EXACT TEST)
# ============================================================================

# Create binary mutation matrix (1 = mutated, 0 = not mutated)
binary_mat <- (mat != "") * 1

# Skip co-mutation analysis if fewer than 2 genes
if (nrow(mat) < 2) {
  cat("\nSkipping co-mutation analysis (requires at least 2 genes, found", nrow(mat), ")\n")
  fisher_results <- data.frame()
  significant_results <- data.frame()
} else {

# Function to perform pairwise Fisher's exact test
pairwise_fisher <- function(binary_matrix) {
  genes <- rownames(binary_matrix)
  n_genes <- length(genes)
  n_samples <- ncol(binary_matrix)

  # Initialize results data frame
  results <- data.frame(
    Gene1 = character(),
    Gene2 = character(),
    Gene1_mutated = integer(),
    Gene2_mutated = integer(),
    Both_mutated = integer(),
    Neither_mutated = integer(),
    Odds_Ratio = numeric(),
    P_value = numeric(),
    Direction = character(),
    stringsAsFactors = FALSE
  )

  # Perform pairwise tests
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      gene1 <- genes[i]
      gene2 <- genes[j]

      # Create 2x2 contingency table
      both <- sum(binary_matrix[gene1, ] == 1 & binary_matrix[gene2, ] == 1)
      gene1_only <- sum(binary_matrix[gene1, ] == 1 & binary_matrix[gene2, ] == 0)
      gene2_only <- sum(binary_matrix[gene1, ] == 0 & binary_matrix[gene2, ] == 1)
      neither <- sum(binary_matrix[gene1, ] == 0 & binary_matrix[gene2, ] == 0)

      # Contingency table
      cont_table <- matrix(c(both, gene1_only, gene2_only, neither),
                           nrow = 2, byrow = TRUE)

      # Fisher's exact test
      fisher_result <- fisher.test(cont_table)

      # Determine direction
      direction <- ifelse(fisher_result$estimate > 1, "Co-occurrence",
                          ifelse(fisher_result$estimate < 1, "Mutual_Exclusivity",
                                 "No_Association"))

      results <- rbind(results, data.frame(
        Gene1 = gene1,
        Gene2 = gene2,
        Gene1_mutated = gene1_only + both,
        Gene2_mutated = gene2_only + both,
        Both_mutated = both,
        Neither_mutated = neither,
        Odds_Ratio = fisher_result$estimate,
        P_value = fisher_result$p.value,
        Direction = direction,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Adjust p-values for multiple testing (Benjamini-Hochberg)
  results$P_adjusted <- p.adjust(results$P_value, method = "BH")

  # Sort by p-value
  results <- results[order(results$P_value), ]

  return(results)
}

# Run pairwise Fisher's test
cat("\nPerforming pairwise Fisher's exact tests for co-mutation analysis...\n")
fisher_results <- pairwise_fisher(binary_mat)

# Save full results
write.csv(fisher_results, file.path(opt$output, "co_mutation_fisher_results.csv"), row.names = FALSE)
cat("Full Fisher's test results saved to:", file.path(opt$output, "co_mutation_fisher_results.csv"), "\n")

# ============================================================================
# SIGNIFICANT CO-MUTATION RESULTS
# ============================================================================

# Filter for significant results (adjusted p < 0.05)
significant_results <- fisher_results %>%
  filter(P_adjusted < 0.05)

cat("\nSignificant co-mutation pairs (adjusted p < 0.05):",
    nrow(significant_results), "\n")

if (nrow(significant_results) > 0) {
  cat("\nTop significant pairs:\n")
  print(head(significant_results, 20))

  write.csv(significant_results, file.path(opt$output, "significant_co_mutations.csv"), row.names = FALSE)
  cat("\nSignificant results saved to:", file.path(opt$output, "significant_co_mutations.csv"), "\n")
}

# ============================================================================
# CO-MUTATION HEATMAP (OPTIONAL VISUALIZATION)
# ============================================================================

# Create a co-occurrence matrix for visualization
create_cooccurrence_matrix <- function(fisher_results, genes) {
  n <- length(genes)
  co_mat <- matrix(0, nrow = n, ncol = n)
  rownames(co_mat) <- genes
  colnames(co_mat) <- genes

  # Fill matrix with -log10(p-value) * sign(log(OR))
  for (i in 1:nrow(fisher_results)) {
    g1 <- fisher_results$Gene1[i]
    g2 <- fisher_results$Gene2[i]

    if (g1 %in% genes && g2 %in% genes) {
      # Use signed -log10(p) to show direction
      value <- -log10(fisher_results$P_value[i] + 1e-300)
      if (fisher_results$Odds_Ratio[i] < 1) value <- -value

      co_mat[g1, g2] <- value
      co_mat[g2, g1] <- value
    }
  }

  return(co_mat)
}

# Create and plot co-occurrence heatmap
co_mat <- create_cooccurrence_matrix(fisher_results, rownames(binary_mat))

# Define color scale (blue = mutual exclusivity, red = co-occurrence)
max_val <- max(abs(co_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_val, 0, max_val), c("#3498db", "white", "#e74c3c"))

pdf(file.path(opt$output, "co_mutation_heatmap.pdf"), width = 12, height = 10)

Heatmap(
  co_mat,
  name = "Signed\n-log10(p)",
  col = col_fun,
  column_title = "Co-mutation Analysis\n(Red = Co-occurrence, Blue = Mutual Exclusivity)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE
)

dev.off()

cat("Co-mutation heatmap saved to:", file.path(opt$output, "co_mutation_heatmap.pdf"), "\n")

} # End of else block for co-mutation analysis (requires >= 2 genes)

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n============================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("============================================================\n")
cat("Total mutations after filtering:", nrow(filtered_data), "\n")
cat("Unique samples:", n_distinct(filtered_data$base_sample), "\n")
cat("Unique genes:", nrow(binary_mat), "\n")
cat("Total pairwise comparisons:", nrow(fisher_results), "\n")
cat("Significant pairs (adj. p < 0.05):", nrow(significant_results), "\n")

# Top mutated genes
cat("\nTop 10 most frequently mutated genes:\n")
gene_summary <- filtered_data %>%
  group_by(gene) %>%
  summarise(
    n_mutations = n(),
    n_samples = n_distinct(base_sample),
    pct_samples = n_distinct(base_sample) / n_distinct(filtered_data$base_sample) * 100
  ) %>%
  arrange(desc(n_samples))

print(head(gene_summary, 10))

# Save gene summary
write.csv(gene_summary, file.path(opt$output, "gene_mutation_summary.csv"), row.names = FALSE)
cat("\nGene mutation summary saved to:", file.path(opt$output, "gene_mutation_summary.csv"), "\n")

cat("\n============================================================\n")
cat("Analysis complete! Output files in:", opt$output, "\n")
cat("  - CHIP_oncoplot.pdf/png\n")
cat("  - VAF_heatmap.pdf/png\n")
cat("  - VAF_dotplot.pdf/png\n")
cat("  - co_mutation_heatmap.pdf\n")
cat("  - co_mutation_fisher_results.csv\n")
cat("  - significant_co_mutations.csv\n")
cat("  - gene_mutation_summary.csv\n")
cat("  - variant_annotation_summary.csv\n")
cat("============================================================\n")
