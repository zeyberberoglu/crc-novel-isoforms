#!/bin/bash
#SBATCH --job-name=gene_zscore_heatmap
#SBATCH --output=logs/gene_zscore_heatmap_%j.out
#SBATCH --error=logs/gene_zscore_heatmap_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH --partition=compute

set -euo pipefail

# ============================================================
# Working directory and environment
# ============================================================

BASE_DIR="/path/to/cell_lines"

cd "${BASE_DIR}"

module load miniforge
conda activate longread_de

# ============================================================
# Input and output paths
# ============================================================

SALMON_DIR="${BASE_DIR}/salmon_quant"

ANNOTATION_FILE="${BASE_DIR}/cell_line_annotations (1).csv"

GFF_FILE="/path/to/final_extended_annotation.gff"

OUT_DIR="${BASE_DIR}/gene_level_zscore_heatmap"

# Number of most variable genes to show.
TOP_N=100

mkdir -p "${OUT_DIR}" "${BASE_DIR}/logs"

export SALMON_DIR
export ANNOTATION_FILE
export GFF_FILE
export OUT_DIR
export TOP_N

# ============================================================
# Write the R script
# ============================================================

cat > "${OUT_DIR}/make_gene_level_zscore_heatmap.R" <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(pheatmap)
  library(RColorBrewer)
})

salmon_dir      <- Sys.getenv("SALMON_DIR")
annotation_file <- Sys.getenv("ANNOTATION_FILE")
gff_file        <- Sys.getenv("GFF_FILE")
out_dir         <- Sys.getenv("OUT_DIR")
top_n           <- as.integer(Sys.getenv("TOP_N"))

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Salmon directory: ", salmon_dir)
message("Annotation file: ", annotation_file)
message("GFF file: ", gff_file)
message("Output directory: ", out_dir)
message("Top variable genes requested: ", top_n)

# ============================================================
# 1. Check inputs
# ============================================================

if (!dir.exists(salmon_dir)) {
  stop("Salmon directory does not exist: ", salmon_dir)
}

if (!file.exists(annotation_file)) {
  stop("Annotation file does not exist: ", annotation_file)
}

if (!file.exists(gff_file)) {
  stop("GFF file does not exist: ", gff_file)
}

# ============================================================
# 2. Find Salmon quant.sf files
# ============================================================

quant_files <- list.files(
  salmon_dir,
  pattern = "^quant\\.sf$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(quant_files) == 0) {
  stop("No quant.sf files were found under: ", salmon_dir)
}

sample_names <- basename(
  dirname(quant_files)
)

names(quant_files) <- sample_names

if (anyDuplicated(sample_names)) {

  duplicated_samples <- unique(
    sample_names[duplicated(sample_names)]
  )

  stop(
    "Duplicated Salmon sample folder names: ",
    paste(duplicated_samples, collapse = ", ")
  )
}

message("Completed Salmon samples found: ", length(quant_files))

# ============================================================
# 3. Import transcript-to-gene mapping from final GFF
# ============================================================

message("Importing final_extended_annotation.gff...")

gff <- import(gff_file)

feature_type <- as.character(
  mcols(gff)$type
)

transcript_features <- gff[
  feature_type %in% c(
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA"
  )
]

if (length(transcript_features) == 0) {
  stop(
    "No transcript-level features were found in the GFF. ",
    "Check the feature types in final_extended_annotation.gff."
  )
}

get_first_attribute <- function(gr, possible_names) {

  available_names <- intersect(
    possible_names,
    colnames(mcols(gr))
  )

  if (length(available_names) == 0) {
    return(
      rep(NA_character_, length(gr))
    )
  }

  values <- as.character(
    mcols(gr)[[available_names[1]]]
  )

  values[values == ""] <- NA_character_

  values
}

transcript_id <- get_first_attribute(
  transcript_features,
  c(
    "transcript_id",
    "transcriptId",
    "ID",
    "Name"
  )
)

gene_id <- get_first_attribute(
  transcript_features,
  c(
    "gene_id",
    "geneID",
    "gene",
    "Parent"
  )
)

gene_name <- get_first_attribute(
  transcript_features,
  c(
    "gene_name",
    "geneName",
    "gene_symbol"
  )
)

# Remove common GFF prefixes.
transcript_id <- sub(
  "^transcript:",
  "",
  transcript_id
)

gene_id <- sub(
  "^gene:",
  "",
  gene_id
)

gene_id <- sub(
  "^transcript:",
  "",
  gene_id
)

tx2gene <- data.table(
  transcript_id = transcript_id,
  gene_id = gene_id,
  gene_name = gene_name
)

tx2gene <- tx2gene[
  !is.na(transcript_id) &
  transcript_id != ""
]

# Prefer gene name, then gene ID.
tx2gene[
  is.na(gene_name) | gene_name == "",
  gene_name := gene_id
]

# Final fallback.
tx2gene[
  is.na(gene_name) | gene_name == "",
  gene_name := transcript_id
]

tx2gene[
  is.na(gene_id) | gene_id == "",
  gene_id := gene_name
]

tx2gene <- unique(
  tx2gene,
  by = "transcript_id"
)

message(
  "Transcript-to-gene mappings obtained: ",
  nrow(tx2gene)
)

fwrite(
  tx2gene,
  file.path(
    out_dir,
    "transcript_to_gene_mapping_used.csv"
  )
)

# ============================================================
# 4. Read Salmon TPM and sum transcripts to genes
# ============================================================

gene_expression_list <- vector(
  mode = "list",
  length = length(quant_files)
)

mapping_summary_list <- vector(
  mode = "list",
  length = length(quant_files)
)

for (i in seq_along(quant_files)) {

  sample_name <- names(quant_files)[i]
  quant_file  <- quant_files[i]

  message(
    "[",
    i,
    "/",
    length(quant_files),
    "] Reading ",
    sample_name
  )

  quant <- fread(
    quant_file,
    select = c(
      "Name",
      "TPM"
    )
  )

  setnames(
    quant,
    old = c(
      "Name",
      "TPM"
    ),
    new = c(
      "transcript_id",
      "TPM"
    )
  )

  quant[
    ,
    transcript_id := sub(
      "^transcript:",
      "",
      transcript_id
    )
  ]

  quant <- merge(
    quant,
    tx2gene,
    by = "transcript_id",
    all.x = TRUE,
    sort = FALSE
  )

  mapping_summary_list[[i]] <- data.table(
    sample = sample_name,
    total_salmon_transcripts = nrow(quant),
    mapped_transcripts = sum(!is.na(quant$gene_name)),
    unmapped_transcripts = sum(is.na(quant$gene_name))
  )

  quant <- quant[
    !is.na(gene_name) &
    gene_name != ""
  ]

  # Sum all transcript TPMs belonging to the same gene.
  gene_quant <- quant[
    ,
    .(
      TPM = sum(TPM, na.rm = TRUE)
    ),
    by = gene_name
  ]

  gene_quant[
    ,
    sample := sample_name
  ]

  gene_expression_list[[i]] <- gene_quant
}

mapping_summary <- rbindlist(
  mapping_summary_list
)

fwrite(
  mapping_summary,
  file.path(
    out_dir,
    "salmon_mapping_summary.csv"
  )
)

gene_expression_long <- rbindlist(
  gene_expression_list,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================
# 5. Build gene × cell-line TPM matrix
# ============================================================

gene_tpm_table <- dcast(
  gene_expression_long,
  gene_name ~ sample,
  value.var = "TPM",
  fun.aggregate = sum,
  fill = 0
)

gene_names <- gene_tpm_table$gene_name

gene_tpm_table[
  ,
  gene_name := NULL
]

gene_tpm_matrix <- as.matrix(
  gene_tpm_table
)

rownames(gene_tpm_matrix) <- gene_names

storage.mode(gene_tpm_matrix) <- "numeric"

# Remove genes with no expression in any cell line.
gene_tpm_matrix <- gene_tpm_matrix[
  rowSums(
    gene_tpm_matrix,
    na.rm = TRUE
  ) > 0,
  ,
  drop = FALSE
]

message(
  "Expressed genes in gene-level TPM matrix: ",
  nrow(gene_tpm_matrix)
)

fwrite(
  data.table(
    gene_name = rownames(gene_tpm_matrix),
    gene_tpm_matrix
  ),
  file.path(
    out_dir,
    "cellline_gene_TPM_matrix.csv"
  )
)

# ============================================================
# 6. Log2 transformation
# ============================================================

log2_tpm_matrix <- log2(
  gene_tpm_matrix + 1
)

# ============================================================
# 7. Select the most variable genes
# ============================================================

gene_variance <- apply(
  log2_tpm_matrix,
  1,
  var,
  na.rm = TRUE
)

gene_variance[
  !is.finite(gene_variance)
] <- 0

ordered_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)

selected_genes <- head(
  ordered_genes,
  min(
    top_n,
    length(ordered_genes)
  )
)

heatmap_log2_matrix <- log2_tpm_matrix[
  selected_genes,
  ,
  drop = FALSE
]

writeLines(
  selected_genes,
  file.path(
    out_dir,
    "genes_used_in_heatmap.txt"
  )
)

message(
  "Genes selected for the heatmap: ",
  nrow(heatmap_log2_matrix)
)

# ============================================================
# 8. Calculate row z-scores
# ============================================================

zscore_matrix <- t(
  scale(
    t(heatmap_log2_matrix),
    center = TRUE,
    scale = TRUE
  )
)

# Remove genes with undefined z-scores.
valid_rows <- apply(
  zscore_matrix,
  1,
  function(x) {
    all(is.finite(x))
  }
)

if (any(!valid_rows)) {

  removed_genes <- rownames(zscore_matrix)[
    !valid_rows
  ]

  writeLines(
    removed_genes,
    file.path(
      out_dir,
      "genes_removed_due_to_invalid_zscore.txt"
    )
  )
}

zscore_matrix <- zscore_matrix[
  valid_rows,
  ,
  drop = FALSE
]

# Match the scale of the transcript-level heatmap.
zscore_matrix[
  zscore_matrix > 2
] <- 2

zscore_matrix[
  zscore_matrix < -2
] <- -2

fwrite(
  data.table(
    gene_name = rownames(zscore_matrix),
    zscore_matrix
  ),
  file.path(
    out_dir,
    "cellline_gene_zscore_matrix.csv"
  )
)

# ============================================================
# 9. Read exact annotation columns
# ============================================================

annotation <- fread(
  annotation_file,
  check.names = FALSE
)

required_annotation_columns <- c(
  "cell_line",
  "MSI_status",
  "origin",
  "CMS",
  "driver",
  "KRAS_variant"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annotation)
)

if (length(missing_annotation_columns) > 0) {

  stop(
    "These required annotation columns are missing: ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    ),
    "\nColumns found in the file: ",
    paste(
      names(annotation),
      collapse = ", "
    )
  )
}

message(
  "Annotation columns confirmed: ",
  paste(
    required_annotation_columns,
    collapse = ", "
  )
)

# ============================================================
# 10. Standardise cell-line names for matching
# ============================================================

standardise_cell_line <- function(x) {

  x <- trimws(
    as.character(x)
  )

  # Example:
  # CRC_HCT116_PRJNA523380 -> HCT116
  x <- sub(
    "^CRC_",
    "",
    x,
    ignore.case = TRUE
  )

  x <- sub(
    "_PRJNA[0-9]+$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- toupper(x)

  x <- gsub(
    "[^A-Z0-9]",
    "",
    x
  )

  x
}

annotation[
  ,
  sample_key := standardise_cell_line(
    cell_line
  )
]

salmon_sample_table <- data.table(
  salmon_sample = colnames(zscore_matrix)
)

salmon_sample_table[
  ,
  sample_key := standardise_cell_line(
    salmon_sample
  )
]

# Check for duplicated annotation names after standardisation.
duplicated_keys <- annotation[
  duplicated(sample_key) |
  duplicated(sample_key, fromLast = TRUE),
  unique(sample_key)
]

if (length(duplicated_keys) > 0) {

  message(
    "Warning: duplicated standardised annotation names: ",
    paste(
      duplicated_keys,
      collapse = ", "
    )
  )
}

annotation <- unique(
  annotation,
  by = "sample_key"
)

# ============================================================
# 11. Match annotations to Salmon samples
# ============================================================

matched_annotation <- merge(
  salmon_sample_table,
  annotation,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)

# Restore exact expression matrix order.
matched_annotation <- matched_annotation[
  match(
    colnames(zscore_matrix),
    salmon_sample
  )
]

metadata_columns <- c(
  "KRAS_variant",
  "driver",
  "CMS",
  "origin",
  "MSI_status"
)

# Replace blank or missing annotation values with Unknown.
for (column_name in metadata_columns) {

  matched_annotation[
    is.na(get(column_name)) |
    trimws(as.character(get(column_name))) == "",
    (column_name) := "Unknown"
  ]
}

# Report samples that did not match any annotation row.
unmatched_samples <- matched_annotation[
  is.na(cell_line),
  salmon_sample
]

if (length(unmatched_samples) > 0) {

  message(
    "Samples without an annotation match: ",
    paste(
      unmatched_samples,
      collapse = ", "
    )
  )

  writeLines(
    unmatched_samples,
    file.path(
      out_dir,
      "samples_without_annotation_match.txt"
    )
  )

} else {

  message(
    "All Salmon samples matched the annotation file."
  )
}

# Build pheatmap annotation data frame.
annotation_col <- as.data.frame(
  matched_annotation[
    ,
    ..metadata_columns
  ]
)

rownames(annotation_col) <- matched_annotation$salmon_sample

for (column_name in colnames(annotation_col)) {

  annotation_col[[column_name]] <- factor(
    annotation_col[[column_name]]
  )
}

# Ensure row order exactly matches heatmap columns.
annotation_col <- annotation_col[
  colnames(zscore_matrix),
  ,
  drop = FALSE
]

fwrite(
  data.table(
    salmon_sample = rownames(annotation_col),
    annotation_col
  ),
  file.path(
    out_dir,
    "cellline_annotations_used.csv"
  )
)

message(
  "Annotation tracks used in heatmap: ",
  paste(
    colnames(annotation_col),
    collapse = ", "
  )
)

for (column_name in colnames(annotation_col)) {

  message("")
  message("Counts for ", column_name, ":")

  print(
    table(
      annotation_col[[column_name]],
      useNA = "ifany"
    )
  )
}

# ============================================================
# 12. Annotation colours
# ============================================================

make_annotation_palette <- function(
  values,
  palette_name = "Set3"
) {

  categories <- levels(
    factor(values)
  )

  number_categories <- length(categories)

  if (number_categories == 1) {

    colours <- "#BDBDBD"

  } else {

    base_number <- max(
      3,
      min(
        12,
        number_categories
      )
    )

    base_colours <- brewer.pal(
      base_number,
      palette_name
    )

    colours <- colorRampPalette(
      base_colours
    )(
      number_categories
    )
  }

  names(colours) <- categories

  colours
}

annotation_colors <- list(

  KRAS_variant = make_annotation_palette(
    annotation_col$KRAS_variant,
    "Set3"
  ),

  driver = c(
    KRAS = "#E41A1C",
    BRAF = "#377EB8",
    WT = "#4DAF4A",
    Unknown = "#BDBDBD"
  ),

  CMS = c(
    CMS1 = "#66C2A5",
    CMS2 = "#FC8D62",
    CMS3 = "#8DA0CB",
    CMS4 = "#E78AC3",
    Unknown = "#BDBDBD"
  ),

  origin = c(
    primary = "#984EA3",
    metastatic = "#F781BF",
    Unknown = "#BDBDBD"
  ),

  MSI_status = c(
    MSS = "#1B9E77",
    MSI = "#D95F02",
    Unknown = "#BDBDBD"
  )
)

# Add colours automatically for any unexpected categories.
for (annotation_name in names(annotation_colors)) {

  observed_categories <- levels(
    annotation_col[[annotation_name]]
  )

  categories_without_colour <- setdiff(
    observed_categories,
    names(
      annotation_colors[[annotation_name]]
    )
  )

  if (length(categories_without_colour) > 0) {

    extra_palette <- make_annotation_palette(
      annotation_col[[annotation_name]],
      "Set3"
    )

    annotation_colors[[annotation_name]][
      categories_without_colour
    ] <- extra_palette[
      categories_without_colour
    ]
  }
}

# ============================================================
# 13. Heatmap colours and dimensions
# ============================================================

heatmap_colours <- colorRampPalette(
  c(
    "#313695",
    "#FFFFFF",
    "#D73027"
  )
)(101)

heatmap_breaks <- seq(
  -2,
  2,
  length.out = 102
)

pdf_width <- max(
  16,
  8 + ncol(zscore_matrix) * 0.22
)

pdf_height <- max(
  13,
  5 + nrow(zscore_matrix) * 0.11
)

# ============================================================
# 14. Plot gene-level row z-score heatmap
# ============================================================

heatmap_arguments <- list(

  mat = zscore_matrix,

  scale = "none",

  color = heatmap_colours,
  breaks = heatmap_breaks,

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",

  annotation_col = annotation_col,
  annotation_colors = annotation_colors,

  annotation_names_col = TRUE,
  annotation_legend = TRUE,

  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize = 8,
  fontsize_row = 4.5,
  fontsize_col = 6,

  angle_col = 90,

  border_color = "#BDBDBD",

  treeheight_row = 50,
  treeheight_col = 50,

  main = "Cell line gene expression heatmap, row Z-score"
)

message("")
message("Creating gene-level z-score PDF heatmap...")

pdf(
  file.path(
    out_dir,
    "cellline_gene_heatmap_zscore.pdf"
  ),
  width = pdf_width,
  height = pdf_height,
  onefile = FALSE
)

do.call(
  pheatmap,
  heatmap_arguments
)

dev.off()

message("Creating gene-level z-score PNG heatmap...")

png(
  file.path(
    out_dir,
    "cellline_gene_heatmap_zscore.png"
  ),
  width = pdf_width,
  height = pdf_height,
  units = "in",
  res = 300
)

do.call(
  pheatmap,
  heatmap_arguments
)

dev.off()

# ============================================================
# 15. Final summary
# ============================================================

summary_table <- data.table(
  item = c(
    "Completed Salmon samples",
    "Expressed genes",
    "Genes shown in heatmap",
    "Annotation tracks"
  ),
  value = c(
    length(quant_files),
    nrow(gene_tpm_matrix),
    nrow(zscore_matrix),
    ncol(annotation_col)
  )
)

fwrite(
  summary_table,
  file.path(
    out_dir,
    "gene_zscore_heatmap_summary.tsv"
  ),
  sep = "\t"
)

message("")
message("Gene-level z-score heatmap completed.")
message("Cell lines plotted: ", ncol(zscore_matrix))
message("Genes plotted: ", nrow(zscore_matrix))
message(
  "Annotation tracks: ",
  paste(
    colnames(annotation_col),
    collapse = ", "
  )
)
message("")
message(
  "PDF: ",
  file.path(
    out_dir,
    "cellline_gene_heatmap_zscore.pdf"
  )
)
message(
  "PNG: ",
  file.path(
    out_dir,
    "cellline_gene_heatmap_zscore.png"
  )
)
EOF

# ============================================================
# Execute R script
# ============================================================

Rscript "${OUT_DIR}/make_gene_level_zscore_heatmap.R"

echo
echo "Gene-level z-score heatmap completed."
echo
echo "Main output:"
echo "${OUT_DIR}/cellline_gene_heatmap_zscore.pdf"
