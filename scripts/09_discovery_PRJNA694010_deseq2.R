Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(tximport)
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

register(MulticoreParam(workers = 8))

quant_dir <- "/path/to/project/diff_exp_salmon/quant/"

sample_names <- list.files(quant_dir)

files <- file.path(
  quant_dir,
  sample_names,
  "quant.sf"
)

names(files) <- sample_names

stopifnot(all(file.exists(files)))

gff <- import(
  "/path/to/project/final_extended_annotation.gff"
)

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
  select(transcript_id, gene_id) %>%
  filter(
    !is.na(transcript_id),
    !is.na(gene_id)
  ) %>%
  distinct()

gene_meta <- as.data.frame(
  mcols(gff[gff$type %in% c("gene", "transcript")])
) %>%
  select(
    gene_id,
    all_of(biotype_col),
    all_of(name_col)
  ) %>%
  filter(!is.na(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

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

stopifnot(
  identical(
    rownames(samples),
    colnames(txi$counts)
  )
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ condition
)

keep <- rowSums(
  counts(dds) >= 10
) >= min(table(samples$condition))

dds <- dds[keep, ]

dds <- DESeq(
  dds,
  parallel = TRUE
)

vsd <- vst(
  dds,
  blind = FALSE
)

pdf(
  "QC_PCA.pdf",
  width = 6,
  height = 5
)

print(
  plotPCA(
    vsd,
    intgroup = "condition"
  ) +
    theme_bw()
)

dev.off()

sample_dist <- dist(
  t(assay(vsd))
)

pdf(
  "QC_sample_distance.pdf",
  width = 10,
  height = 9
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

pdf(
  "QC_dispersion.pdf",
  width = 6,
  height = 5
)

plotDispEsts(dds)

dev.off()

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
    rownames_to_column("gene_id") %>%
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
  "DESeq2_all_genes_shrunken.csv",
  row.names = FALSE
)

write.csv(
  annotate_res(res_lncrna),
  "DESeq2_lncRNA_shrunken.csv",
  row.names = FALSE
)

volcano <- function(res, title, file) {

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
      subtitle = "apeglm-shrunken LFC",
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
  "CRC vs NAT - all genes",
  "Volcano_all_genes.pdf"
)

volcano(
  res_lncrna,
  "CRC vs NAT - lncRNA",
  "Volcano_lncRNA.pdf"
)

saveRDS(
  dds,
  "dds_full.rds"
)

saveRDS(
  res_shrink,
  "res_shrink_full.rds"
)
