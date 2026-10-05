Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(EnhancedVolcano)
  library(rtracklayer)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(tibble)
})

setDTthreads(24)
register(MulticoreParam(workers = 24))

# ==============================
# User-defined paths
# ==============================

quant_dir <- "/path/to/PRJNA387172/salmon_quant/"
out_dir <- "/path/to/PRJNA387172/differential_expression"
gff_path <- "/path/to/final_extended_annotation.gff"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

setwd(out_dir)

# ==============================
# Locate Salmon quantifications
# ==============================

sample_names <- list.files(
  quant_dir,
  pattern = "^(CRC|NAT)"
)

files <- file.path(
  quant_dir,
  sample_names,
  "quant.sf"
)

names(files) <- sample_names

ok <- file.exists(files)

if (any(!ok)) {
  message(
    "Missing quant.sf, dropping: ",
    paste(sample_names[!ok], collapse = ", ")
  )
}

files <- files[ok]
sample_names <- sample_names[ok]

# ==============================
# Import annotation
# ==============================

gff <- import(gff_path)

biotype_col <- if ("gene_type" %in% names(mcols(gff))) {
  "gene_type"
} else {
  "gene_biotype"
}

name_col <- if ("gene_name" %in% names(mcols(gff))) {
  "gene_name"
} else {
  "Name"
}

tx2gene <- as.data.frame(
  mcols(gff[gff$type == "transcript"])
) %>%
  select(
    transcript_id,
    gene_id
  ) %>%
  filter(
    !is.na(transcript_id),
    !is.na(gene_id)
  ) %>%
  distinct()

gene_meta <- as.data.frame(
  mcols(
    gff[gff$type %in% c("gene", "transcript")]
  )
) %>%
  select(
    gene_id,
    all_of(biotype_col),
    all_of(name_col)
  ) %>%
  filter(!is.na(gene_id)) %>%
  distinct(
    gene_id,
    .keep_all = TRUE
  )

saveRDS(
  tx2gene,
  "tx2gene.rds"
)

saveRDS(
  gene_meta,
  "gene_meta.rds"
)

# ==============================
# Build paired sample metadata
# ==============================

samples <- data.frame(
  row.names = sample_names,
  condition = factor(
    ifelse(
      grepl("^CRC", sample_names),
      "CRC",
      "NAT"
    ),
    levels = c("NAT", "CRC")
  ),
  patient = sub(
    "^(CRC|NAT)_",
    "",
    sample_names
  )
)

paired_patients <- intersect(
  samples$patient[
    samples$condition == "CRC"
  ],
  samples$patient[
    samples$condition == "NAT"
  ]
)

if (length(paired_patients) < 2) {
  stop(
    "Not enough complete CRC/NAT pairs for paired DESeq2 analysis."
  )
}

unpaired <- setdiff(
  samples$patient,
  paired_patients
)

if (length(unpaired)) {
  message(
    "Dropping unpaired patients: ",
    paste(
      sort(unpaired),
      collapse = ", "
    )
  )
}

samples <- samples[
  samples$patient %in% paired_patients,
]

samples$patient <- factor(
  samples$patient
)

samples <- samples[
  order(
    samples$patient,
    samples$condition
  ),
]

message(
  sprintf(
    "Complete pairs: %d patients, %d samples",
    length(paired_patients),
    nrow(samples)
  )
)

saveRDS(
  samples,
  "samples_metadata.rds"
)

# ==============================
# tximport + DESeq2
# ==============================

files <- files[
  rownames(samples)
]

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

stopifnot(
  identical(
    rownames(samples),
    colnames(txi$counts)
  )
)

saveRDS(
  txi,
  "txi_salmon.rds"
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ patient + condition
)

keep <- rowSums(
  counts(dds) >= 10
) >= length(paired_patients)

dds <- dds[
  keep,
]

dds <- DESeq(
  dds,
  parallel = TRUE
)

vsd <- vst(
  dds,
  blind = FALSE
)

saveRDS(
  vsd,
  "vsd.rds"
)

# ==============================
# VST PCA
# ==============================

vsd_mat <- assay(vsd)

rv_vsd <- apply(
  vsd_mat,
  1,
  var
)

select_vsd <- order(
  rv_vsd,
  decreasing = TRUE
)[
  seq_len(
    min(
      500,
      length(rv_vsd)
    )
  )
]

vsd_top <- vsd_mat[
  select_vsd,
]

pca_vsd <- prcomp(
  t(vsd_top),
  scale. = TRUE
)

pca_vsd_data <- data.frame(
  PC1 = pca_vsd$x[, 1],
  PC2 = pca_vsd$x[, 2],
  condition = samples$condition
)

pdf(
  "QC_PCA_VST_scaled.pdf",
  width = 6,
  height = 5
)

print(
  ggplot(
    pca_vsd_data,
    aes(
      x = PC1,
      y = PC2,
      color = condition
    )
  ) +
    geom_point(size = 3) +
    theme_bw() +
    ggtitle(
      "PCA on Scaled VST Data"
    )
)

dev.off()

# ==============================
# TPM PCA
# ==============================

tpm_mat <- txi$abundance

rv_tpm <- apply(
  tpm_mat,
  1,
  var
)

select_tpm <- order(
  rv_tpm,
  decreasing = TRUE
)[
  seq_len(
    min(
      500,
      length(rv_tpm)
    )
  )
]

tpm_top <- tpm_mat[
  select_tpm,
]

pca_tpm <- prcomp(
  t(tpm_top),
  scale. = TRUE
)

pca_tpm_data <- data.frame(
  PC1 = pca_tpm$x[, 1],
  PC2 = pca_tpm$x[, 2],
  condition = samples$condition
)

pdf(
  "QC_PCA_TPM_scaled.pdf",
  width = 6,
  height = 5
)

print(
  ggplot(
    pca_tpm_data,
    aes(
      x = PC1,
      y = PC2,
      color = condition
    )
  ) +
    geom_point(size = 3) +
    theme_bw() +
    ggtitle(
      "PCA on Scaled TPM Data"
    )
)

dev.off()

# ==============================
# Sample distance
# ==============================

sample_dist <- dist(
  t(assay(vsd))
)

pdf(
  "QC_sample_distance.pdf",
  width = 14,
  height = 13
)

pheatmap(
  as.matrix(sample_dist),
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  col = colorRampPalette(
    rev(
      brewer.pal(
        9,
        "Blues"
      )
    )
  )(255),
  annotation_col = samples["condition"]
)

dev.off()

# ==============================
# Dispersion
# ==============================

pdf(
  "QC_dispersion.pdf",
  width = 6,
  height = 5
)

plotDispEsts(dds)

dev.off()

# ==============================
# LFC shrinkage
# ==============================

res_shrink <- lfcShrink(
  dds,
  coef = "condition_CRC_vs_NAT",
  type = "apeglm",
  parallel = TRUE
)

pdf(
  "QC_MA_shrunk.pdf",
  width = 6,
  height = 5
)

plotMA(
  res_shrink,
  ylim = c(-5, 5)
)

dev.off()

# ==============================
# Gene annotation
# ==============================

id2name <- setNames(
  gene_meta[[name_col]],
  gene_meta$gene_id
)

id2name[
  is.na(id2name) |
    id2name == ""
] <- names(id2name)[
  is.na(id2name) |
    id2name == ""
]

get_labels <- function(res) {

  unname(
    ifelse(
      rownames(res) %in% names(id2name),
      id2name[rownames(res)],
      rownames(res)
    )
  )

}

lncrna_types <- c(
  "lncRNA",
  "lincRNA",
  "antisense",
  "antisense_RNA",
  "macro_lncRNA",
  "bidirectional_promoter_lncRNA"
)

lncrna_ids <- gene_meta$gene_id[
  gene_meta[[biotype_col]] %in% lncrna_types
]

annotate_res <- function(res) {

  as.data.frame(res) %>%
    rownames_to_column(
      "gene_id"
    ) %>%
    mutate(
      gene_name = unname(
        id2name[gene_id]
      )
    ) %>%
    select(
      gene_id,
      gene_name,
      everything()
    )

}

res_lncrna <- res_shrink[
  rownames(res_shrink) %in% lncrna_ids,
]

write.csv(
  annotate_res(res_shrink),
  "DESeq2_paired_all_genes_shrunken.csv",
  row.names = FALSE
)

write.csv(
  annotate_res(res_lncrna),
  "DESeq2_paired_lncRNA_shrunken.csv",
  row.names = FALSE
)

# ==============================
# Volcano plots
# ==============================

volcano <- function(
  res,
  title,
  file
) {

  pdf(
    file,
    width = 7,
    height = 8
  )

  print(
    EnhancedVolcano(
      res,
      lab = get_labels(res),
      x = "log2FoldChange",
      y = "padj",
      pCutoff = 0.05,
      FCcutoff = 1,
      title = title,
      subtitle = "paired, apeglm-shrunken LFC",
      pointSize = 1.5,
      labSize = 3,
      drawConnectors = TRUE,
      widthConnectors = 0.3,
      colAlpha = 0.6,
      max.overlaps = 30
    )
  )

  dev.off()

}

volcano(
  res_shrink,
  sprintf(
    "CRC vs NAT - paired all genes (n=%d)",
    length(paired_patients)
  ),
  "Volcano_paired_all_genes.pdf"
)

volcano(
  res_lncrna,
  sprintf(
    "CRC vs NAT - paired lncRNA (n=%d)",
    length(paired_patients)
  ),
  "Volcano_paired_lncRNA.pdf"
)

saveRDS(
  dds,
  "dds_paired.rds"
)

saveRDS(
  res_shrink,
  "res_shrink_paired.rds"
)
