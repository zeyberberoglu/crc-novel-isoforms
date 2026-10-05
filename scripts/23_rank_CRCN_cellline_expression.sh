#!/bin/bash
#SBATCH --job-name=CRCN_cellline_rank
#SBATCH --output=logs/CRCN_cellline_rank_%j.out
#SBATCH --error=logs/CRCN_cellline_rank_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=compute

set -euo pipefail

module load miniforge
conda activate longread_de

BASE="/path/to/cell_lines"

INPUT="${BASE}/gene_level_absolute_TPM_heatmap/cellline_lncRNA_gene_TPM_matrix.csv"
OUTDIR="${BASE}/CRCN_cellline_expression_ranking"

mkdir -p "${OUTDIR}" "${BASE}/logs"

export INPUT OUTDIR

cat > "${OUTDIR}/rank_CRCN_cellline_expression.R" <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_file <- Sys.getenv("INPUT")
outdir <- Sys.getenv("OUTDIR")

if (!file.exists(input_file)) {
  stop("Input matrix does not exist: ", input_file)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Seven selected novel CRCN lncRNAs
targets <- c(
  "CRCN_LINC00280",
  "CRCN_LINC00335",
  "CRCN_FANCC-AS1",
  "CRCN_LINC00128",
  "CRCN_DENND3-AS1",
  "CRCN_KANTR-OT1",
  "CRCN_LINC00205"
)

dat <- fread(input_file, check.names = FALSE)

if (ncol(dat) < 2) {
  stop("The expression matrix contains fewer than two columns.")
}

cat("Input dimensions:", nrow(dat), "rows x", ncol(dat), "columns\n")
cat("Input columns:\n")
print(names(dat))

# Identify gene-name column
possible_gene_columns <- c(
  "gene_name",
  "Gene",
  "gene",
  "Gene_name",
  "gene_symbol",
  "GeneSymbol",
  "symbol"
)

gene_col <- possible_gene_columns[
  possible_gene_columns %in% names(dat)
][1]

if (is.na(gene_col)) {
  # Fall back to first column
  gene_col <- names(dat)[1]
  cat("No standard gene column name detected; using first column:", gene_col, "\n")
} else {
  cat("Detected gene column:", gene_col, "\n")
}

dat[, (gene_col) := as.character(get(gene_col))]

# Sample columns are all remaining numeric columns
candidate_sample_cols <- setdiff(names(dat), gene_col)

sample_cols <- candidate_sample_cols[
  vapply(dat[, ..candidate_sample_cols], is.numeric, logical(1))
]

# In case fread interpreted expression columns as text, attempt conversion
if (length(sample_cols) == 0) {
  for (column_name in candidate_sample_cols) {
    converted <- suppressWarnings(as.numeric(dat[[column_name]]))

    if (sum(!is.na(converted)) > 0) {
      dat[, (column_name) := converted]
    }
  }

  sample_cols <- candidate_sample_cols[
    vapply(dat[, ..candidate_sample_cols], is.numeric, logical(1))
  ]
}

if (length(sample_cols) == 0) {
  stop("No numeric cell-line expression columns were detected.")
}

cat("Number of detected cell-line columns:", length(sample_cols), "\n")

# Standardize alternative dash characters and whitespace
dat[, standardized_gene := trimws(get(gene_col))]
dat[, standardized_gene := gsub("\u2013|\u2014|\u2212", "-", standardized_gene)]

targets_standardized <- gsub("\u2013|\u2014|\u2212", "-", targets)

selected <- dat[standardized_gene %in% targets_standardized]

# Report candidates not found
found_genes <- unique(selected$standardized_gene)
missing_genes <- setdiff(targets_standardized, found_genes)

detection_report <- data.table(
  Gene = targets_standardized,
  Found_in_matrix = targets_standardized %in% found_genes
)

fwrite(
  detection_report,
  file.path(outdir, "CRCN_candidate_detection_report.csv")
)

cat("\nCandidate detection report:\n")
print(detection_report)

if (length(missing_genes) > 0) {
  warning(
    "The following genes were not found in the TPM matrix: ",
    paste(missing_genes, collapse = ", ")
  )
}

if (nrow(selected) == 0) {
  stop(
    "None of the seven CRCN genes were found. ",
    "Check CRCN_candidate_detection_report.csv and the gene identifiers."
  )
}

# If duplicate rows exist for a gene, sum TPM values to obtain one gene-level row
selected_gene <- selected[
  ,
  lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
  by = standardized_gene,
  .SDcols = sample_cols
]

setnames(selected_gene, "standardized_gene", "Gene")

# Create numeric matrix
expr_matrix <- as.matrix(selected_gene[, ..sample_cols])
storage.mode(expr_matrix) <- "numeric"

# Summary statistics
ranking <- data.table(
  Gene = selected_gene$Gene,
  Mean_TPM = rowMeans(expr_matrix, na.rm = TRUE),
  Median_TPM = apply(expr_matrix, 1, median, na.rm = TRUE),
  Max_TPM = apply(expr_matrix, 1, max, na.rm = TRUE),
  Min_TPM = apply(expr_matrix, 1, min, na.rm = TRUE),
  SD_TPM = apply(expr_matrix, 1, sd, na.rm = TRUE),
  Expressed_TPM_gt_0 = rowSums(expr_matrix > 0, na.rm = TRUE),
  Expressed_TPM_gt_1 = rowSums(expr_matrix > 1, na.rm = TRUE),
  Expressed_TPM_gt_5 = rowSums(expr_matrix > 5, na.rm = TRUE),
  Total_Cell_Lines = ncol(expr_matrix)
)

ranking[, Highest_Expressing_Cell_Line :=
          sample_cols[max.col(expr_matrix, ties.method = "first")]]

ranking[, Percent_Cell_Lines_TPM_gt_1 :=
          100 * Expressed_TPM_gt_1 / Total_Cell_Lines]

ranking[, Percent_Cell_Lines_TPM_gt_5 :=
          100 * Expressed_TPM_gt_5 / Total_Cell_Lines]

setorder(ranking, -Mean_TPM, -Median_TPM)

ranking[, Mean_TPM := round(Mean_TPM, 4)]
ranking[, Median_TPM := round(Median_TPM, 4)]
ranking[, Max_TPM := round(Max_TPM, 4)]
ranking[, Min_TPM := round(Min_TPM, 4)]
ranking[, SD_TPM := round(SD_TPM, 4)]
ranking[, Percent_Cell_Lines_TPM_gt_1 :=
          round(Percent_Cell_Lines_TPM_gt_1, 1)]
ranking[, Percent_Cell_Lines_TPM_gt_5 :=
          round(Percent_Cell_Lines_TPM_gt_5, 1)]

fwrite(
  ranking,
  file.path(outdir, "CRCN_lncRNAs_ranked_by_mean_cellline_TPM.csv")
)

# Save the TPM values for the seven genes
selected_output <- copy(selected_gene)
setcolorder(selected_output, c("Gene", sample_cols))

fwrite(
  selected_output,
  file.path(outdir, "CRCN_lncRNA_cellline_TPM_matrix.csv")
)

# Long-format table: useful for examining the best model for each gene
long_expression <- melt(
  selected_output,
  id.vars = "Gene",
  variable.name = "Cell_Line",
  value.name = "TPM"
)

setorder(long_expression, Gene, -TPM)

fwrite(
  long_expression,
  file.path(outdir, "CRCN_lncRNA_expression_long_format.csv")
)

# Top five cell lines for each candidate
top_cell_lines <- long_expression[
  order(Gene, -TPM),
  head(.SD, 5),
  by = Gene
]

fwrite(
  top_cell_lines,
  file.path(outdir, "top5_cell_lines_for_each_CRCN_lncRNA.csv")
)

# Plot 1: mean TPM with standard deviation
plot_data <- copy(ranking)
plot_data[, Gene := factor(Gene, levels = rev(Gene))]

p1 <- ggplot(
  plot_data,
  aes(x = Gene, y = Mean_TPM)
) +
  geom_col(width = 0.72, fill = "#4472C4") +
  geom_errorbar(
    aes(
      ymin = pmax(Mean_TPM - SD_TPM, 0),
      ymax = Mean_TPM + SD_TPM
    ),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = round(Mean_TPM, 2)),
    hjust = -0.15,
    size = 3.6
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Expression of novel CRCN lncRNAs across CRC cell lines",
    subtitle = "Bars show mean TPM; error bars show standard deviation",
    x = NULL,
    y = "Mean TPM"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic"),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave(
  file.path(outdir, "CRCN_lncRNAs_ranked_by_mean_TPM.pdf"),
  p1,
  width = 9,
  height = 6
)

ggsave(
  file.path(outdir, "CRCN_lncRNAs_ranked_by_mean_TPM.png"),
  p1,
  width = 9,
  height = 6,
  dpi = 300
)

# Plot 2: log2-transformed TPM distributions across cell lines
long_expression[, log2_TPM_plus_1 := log2(TPM + 1)]

gene_order <- ranking$Gene
long_expression[, Gene := factor(Gene, levels = rev(gene_order))]

p2 <- ggplot(
  long_expression,
  aes(x = Gene, y = log2_TPM_plus_1)
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    fill = "#A9C4EB"
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 1.5,
    alpha = 0.65
  ) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    title = "Distribution of novel CRCN lncRNA expression",
    subtitle = "Each point represents one colorectal cancer cell line",
    x = NULL,
    y = expression(log[2](TPM + 1))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic")
  )

ggsave(
  file.path(outdir, "CRCN_lncRNA_cellline_expression_distributions.pdf"),
  p2,
  width = 9,
  height = 6
)

ggsave(
  file.path(outdir, "CRCN_lncRNA_cellline_expression_distributions.png"),
  p2,
  width = 9,
  height = 6,
  dpi = 300
)

cat("\n============================================\n")
cat("Analysis completed successfully.\n")
cat("Output directory:", outdir, "\n")
cat("============================================\n\n")

cat("Ranked candidates:\n")
print(ranking)
EOF

Rscript "${OUTDIR}/rank_CRCN_cellline_expression.R"
