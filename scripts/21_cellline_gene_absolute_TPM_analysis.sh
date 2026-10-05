#!/bin/bash
#SBATCH --job-name=gene_TPM_heatmap
#SBATCH --output=logs/gene_TPM_heatmap_%j.out
#SBATCH --error=logs/gene_TPM_heatmap_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH --partition=compute

set -euo pipefail

BASE_DIR="/path/to/cell_lines"

cd "${BASE_DIR}"

module load miniforge
conda activate longread_de

SALMON_DIR="${BASE_DIR}/salmon_quant"

ANNOTATION_FILE="${BASE_DIR}/cell_line_annotations (1).csv"

GFF_FILE="/path/to/final_extended_annotation.gff"

OUT_DIR="${BASE_DIR}/gene_level_absolute_TPM_heatmap"

TOP_N=40

mkdir -p "${OUT_DIR}" "${BASE_DIR}/logs"

export SALMON_DIR
export ANNOTATION_FILE
export GFF_FILE
export OUT_DIR
export TOP_N

cat > "${OUT_DIR}/make_gene_absolute_TPM_heatmap.R" <<'EOF'
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
  stop("No quant.sf files found in: ", salmon_dir)
}

sample_names <- basename(dirname(quant_files))
names(quant_files) <- sample_names

if (anyDuplicated(sample_names)) {
  stop("Duplicated Salmon sample directory names were detected.")
}

message("Completed Salmon samples: ", length(quant_files))

# ============================================================
# 3. Import transcript and gene information from the GFF
# ============================================================

message("Importing final extended annotation...")

gff <- import(gff_file)

feature_type <- as.character(mcols(gff)$type)

transcript_features <- gff[
  feature_type %in% c(
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA"
  )
]

if (length(transcript_features) == 0) {
  stop("No transcript-level features were found in the GFF.")
}

get_first_attribute <- function(gr, possible_names) {

  matching_names <- intersect(
    possible_names,
    colnames(mcols(gr))
  )

  if (length(matching_names) == 0) {
    return(rep(NA_character_, length(gr)))
  }

  result <- as.character(
    mcols(gr)[[matching_names[1]]]
  )

  result[result == ""] <- NA_character_

  result
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

# Search possible biotype columns.
gene_biotype <- get_first_attribute(
  transcript_features,
  c(
    "gene_type",
    "gene_biotype",
    "geneType",
    "biotype",
    "transcript_type",
    "transcript_biotype"
  )
)

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
  gene_name = gene_name,
  gene_biotype = gene_biotype
)

tx2gene <- tx2gene[
  !is.na(transcript_id) &
  transcript_id != ""
]

# Prefer gene name, otherwise use gene ID.
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

tx2gene[
  is.na(gene_biotype) | gene_biotype == "",
  gene_biotype := "Unknown"
]

tx2gene <- unique(
  tx2gene,
  by = "transcript_id"
)

message("Transcript-to-gene mappings: ", nrow(tx2gene))

message(
  "Biotypes detected: ",
  paste(
    sort(unique(tx2gene$gene_biotype)),
    collapse = ", "
  )
)

fwrite(
  tx2gene,
  file.path(
    out_dir,
    "transcript_to_gene_mapping_used.csv"
  )
)

# ============================================================
# 4. Define lncRNA biotypes
# ============================================================

lncrna_biotypes <- c(
  "lncRNA",
  "lincRNA",
  "long_noncoding",
  "long_non_coding",
  "antisense",
  "sense_intronic",
  "sense_overlapping",
  "processed_transcript",
  "3prime_overlapping_ncRNA",
  "bidirectional_promoter_lncRNA",
  "macro_lncRNA",
  "non_coding",
  "ncRNA"
)

tx2gene[
  ,
  is_lncRNA := (
    gene_biotype %in% lncrna_biotypes |
    grepl(
      "lnc|linc|antisense|non.?coding",
      gene_biotype,
      ignore.case = TRUE
    ) |
    grepl(
      "^CRCN_",
      gene_name,
      ignore.case = TRUE
    )
  )
]

message(
  "Transcripts classified as lncRNA: ",
  sum(tx2gene$is_lncRNA)
)

# ============================================================
# 5. Read Salmon TPM and aggregate transcripts to genes
# ============================================================

gene_expression_list <- vector(
  "list",
  length(quant_files)
)

mapping_summary_list <- vector(
  "list",
  length(quant_files)
)

for (i in seq_along(quant_files)) {

  sample_name <- names(quant_files)[i]

  message(
    "[",
    i,
    "/",
    length(quant_files),
    "] ",
    sample_name
  )

  quant <- fread(
    quant_files[i],
    select = c("Name", "TPM")
  )

  setnames(
    quant,
    c("Name", "TPM"),
    c("transcript_id", "TPM")
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
    total_transcripts = nrow(quant),
    mapped_transcripts = sum(!is.na(quant$gene_name)),
    unmapped_transcripts = sum(is.na(quant$gene_name)),
    mapped_lncRNA_transcripts = sum(
      quant$is_lncRNA %in% TRUE,
      na.rm = TRUE
    )
  )

  # Keep only transcripts classified as lncRNA.
  quant <- quant[
    !is.na(gene_name) &
    gene_name != "" &
    is_lncRNA %in% TRUE
  ]

  # Sum transcript TPM values belonging to the same gene.
  gene_quant <- quant[
    ,
    .(
      TPM = sum(TPM, na.rm = TRUE),
      gene_id = paste(
        unique(gene_id),
        collapse = ";"
      ),
      gene_biotype = paste(
        unique(gene_biotype),
        collapse = ";"
      )
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
    "salmon_lncRNA_mapping_summary.csv"
  )
)

gene_expression_long <- rbindlist(
  gene_expression_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(gene_expression_long) == 0) {
  stop(
    "No lncRNA expression records remained after filtering. ",
    "Inspect transcript_to_gene_mapping_used.csv and the GFF biotype columns."
  )
}

# ============================================================
# 6. Build lncRNA gene × cell-line TPM matrix
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

# Remove genes with zero TPM in every cell line.
gene_tpm_matrix <- gene_tpm_matrix[
  rowSums(
    gene_tpm_matrix,
    na.rm = TRUE
  ) > 0,
  ,
  drop = FALSE
]

message(
  "Expressed lncRNA genes: ",
  nrow(gene_tpm_matrix)
)

fwrite(
  data.table(
    gene_name = rownames(gene_tpm_matrix),
    gene_tpm_matrix
  ),
  file.path(
    out_dir,
    "cellline_lncRNA_gene_TPM_matrix.csv"
  )
)

# ============================================================
# 7. Select top 40 lncRNAs by mean TPM
# ============================================================

mean_tpm <- rowMeans(
  gene_tpm_matrix,
  na.rm = TRUE
)

median_tpm <- apply(
  gene_tpm_matrix,
  1,
  median,
  na.rm = TRUE
)

max_tpm <- apply(
  gene_tpm_matrix,
  1,
  max,
  na.rm = TRUE
)

expressed_cell_lines <- rowSums(
  gene_tpm_matrix > 0
)

ranking_table <- data.table(
  gene_name = rownames(gene_tpm_matrix),
  mean_TPM = mean_tpm,
  median_TPM = median_tpm,
  maximum_TPM = max_tpm,
  expressed_cell_lines = expressed_cell_lines
)

setorder(
  ranking_table,
  -mean_TPM
)

fwrite(
  ranking_table,
  file.path(
    out_dir,
    "lncRNA_genes_ranked_by_mean_TPM.csv"
  )
)

number_to_select <- min(
  top_n,
  nrow(ranking_table)
)

top_genes <- ranking_table[
  seq_len(number_to_select),
  gene_name
]

top_40_table <- ranking_table[
  gene_name %in% top_genes
]

fwrite(
  top_40_table,
  file.path(
    out_dir,
    "top40_lncRNAs_by_mean_TPM.csv"
  )
)

writeLines(
  top_genes,
  file.path(
    out_dir,
    "top40_lncRNA_gene_names.txt"
  )
)

heatmap_matrix <- gene_tpm_matrix[
  top_genes,
  ,
  drop = FALSE
]

message(
  "lncRNA genes plotted: ",
  nrow(heatmap_matrix)
)

message(
  "Mean TPM range among selected genes: ",
  round(min(rowMeans(heatmap_matrix)), 4),
  " to ",
  round(max(rowMeans(heatmap_matrix)), 4)
)

fwrite(
  data.table(
    gene_name = rownames(heatmap_matrix),
    heatmap_matrix
  ),
  file.path(
    out_dir,
    "top40_lncRNA_absolute_TPM_matrix.csv"
  )
)

# ============================================================
# 8. Read exact cell-line annotations
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
    "Missing annotation columns: ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
}

# ============================================================
# 9. Match annotation names to Salmon sample names
# ============================================================

standardise_cell_line <- function(x) {

  x <- trimws(as.character(x))

  # CRC_HCT116_PRJNA523380 becomes HCT116.
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

sample_table <- data.table(
  salmon_sample = colnames(heatmap_matrix)
)

sample_table[
  ,
  sample_key := standardise_cell_line(
    salmon_sample
  )
]

annotation <- unique(
  annotation,
  by = "sample_key"
)

matched_annotation <- merge(
  sample_table,
  annotation,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)

matched_annotation <- matched_annotation[
  match(
    colnames(heatmap_matrix),
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

for (column_name in metadata_columns) {

  matched_annotation[
    is.na(get(column_name)) |
    trimws(as.character(get(column_name))) == "",
    (column_name) := "Unknown"
  ]
}

unmatched_samples <- matched_annotation[
  is.na(cell_line),
  salmon_sample
]

if (length(unmatched_samples) > 0) {

  writeLines(
    unmatched_samples,
    file.path(
      out_dir,
      "samples_without_annotation_match.txt"
    )
  )

  message(
    "Samples without annotation matches: ",
    paste(unmatched_samples, collapse = ", ")
  )

} else {

  message("All cell lines matched the annotation file.")
}

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

annotation_col <- annotation_col[
  colnames(heatmap_matrix),
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

# ============================================================
# 10. Annotation colours
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
    )(number_categories)
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

# Add colours for unexpected categories.
for (annotation_name in names(annotation_colors)) {

  observed_categories <- levels(
    annotation_col[[annotation_name]]
  )

  missing_colours <- setdiff(
    observed_categories,
    names(annotation_colors[[annotation_name]])
  )

  if (length(missing_colours) > 0) {

    extra_colours <- make_annotation_palette(
      annotation_col[[annotation_name]],
      "Set3"
    )

    annotation_colors[[annotation_name]][missing_colours] <-
      extra_colours[missing_colours]
  }
}

# ============================================================
# 11. Raw TPM heatmap colour scale
# ============================================================

# These are raw TPM values, with no z-score and no log transform.
#
# Quantile-based breaks prevent one extreme value from making
# the rest of the heatmap appear as a single colour. The matrix
# itself remains raw TPM.

maximum_break <- as.numeric(
  quantile(
    heatmap_matrix,
    probs = 0.99,
    na.rm = TRUE
  )
)

if (!is.finite(maximum_break) || maximum_break <= 0) {
  maximum_break <- max(
    heatmap_matrix,
    na.rm = TRUE
  )
}

if (!is.finite(maximum_break) || maximum_break <= 0) {
  stop("The selected TPM matrix does not contain positive values.")
}

heatmap_breaks <- seq(
  0,
  maximum_break,
  length.out = 102
)

plot_matrix <- heatmap_matrix

# Values above the 99th percentile receive the maximum colour.
# The exported matrix remains completely unchanged.
plot_matrix[
  plot_matrix > maximum_break
] <- maximum_break

heatmap_colours <- colorRampPalette(
  c(
    "#FFFFFF",
    "#FFF7BC",
    "#FEC44F",
    "#D95F0E",
    "#7F0000"
  )
)(101)

pdf_width <- max(
  16,
  8 + ncol(plot_matrix) * 0.22
)

pdf_height <- max(
  11,
  5 + nrow(plot_matrix) * 0.24
)

# ============================================================
# 12. Plot absolute TPM heatmap
# ============================================================

heatmap_arguments <- list(

  mat = plot_matrix,

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
  fontsize_row = 7,
  fontsize_col = 6,

  angle_col = 90,

  border_color = "#D9D9D9",

  treeheight_row = 50,
  treeheight_col = 50,

  main = paste0(
    "Absolute expression of top ",
    nrow(plot_matrix),
    " lncRNA genes across colorectal cancer cell lines (TPM)"
  )
)

message("Creating PDF heatmap...")

pdf(
  file.path(
    out_dir,
    "cellline_top40_lncRNA_absolute_TPM_heatmap.pdf"
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

message("Creating PNG heatmap...")

png(
  file.path(
    out_dir,
    "cellline_top40_lncRNA_absolute_TPM_heatmap.png"
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
# 13. Summary
# ============================================================

summary_table <- data.table(
  item = c(
    "Completed Salmon samples",
    "Expressed lncRNA genes",
    "lncRNAs shown in heatmap",
    "Selection criterion",
    "Expression scale",
    "Colour scale upper limit"
  ),
  value = c(
    length(quant_files),
    nrow(gene_tpm_matrix),
    nrow(heatmap_matrix),
    "Highest mean TPM across cell lines",
    "Raw TPM, no log transformation and no z-score",
    maximum_break
  )
)

fwrite(
  summary_table,
  file.path(
    out_dir,
    "absolute_TPM_heatmap_summary.tsv"
  ),
  sep = "\t"
)

message("")
message("Absolute TPM heatmap completed.")
message("Cell lines plotted: ", ncol(heatmap_matrix))
message("lncRNA genes plotted: ", nrow(heatmap_matrix))
message("Selection: top genes ranked by mean raw TPM")
message("Transformation: none")
message("")
message(
  "PDF: ",
  file.path(
    out_dir,
    "cellline_top40_lncRNA_absolute_TPM_heatmap.pdf"
  )
)
EOF

Rscript "${OUT_DIR}/make_gene_absolute_TPM_heatmap.R"

echo
echo "Absolute gene-level TPM heatmap completed."
echo
echo "Main output:"
echo "${OUT_DIR}/cellline_top40_lncRNA_absolute_TPM_heatmap.pdf"
