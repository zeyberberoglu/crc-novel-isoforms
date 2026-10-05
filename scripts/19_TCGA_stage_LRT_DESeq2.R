suppressMessages({
  library(tximport)
  library(data.table)
  library(DESeq2)
})

# ==============================
# User-defined paths
# ==============================

salmon_dir <- "/path/to/TCGA/salmon_quant"
gff_file   <- "/path/to/final_extended_annotation.gff"
lnc_file   <- "/path/to/robust_upregulated_lncRNAs_for_survival.csv"
clin_file  <- "/path/to/TCGA/survival/clinical_parsed.csv"
out_dir    <- "/path/to/TCGA/survival/stage_LRT_DESeq2"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ==============================
# Helper functions
# ==============================

strip <- function(x) {

  ifelse(
    grepl("^ENS", x),
    sub("\\.[0-9]+$", "", x),
    x
  )
}

clean_stage <- function(x) {

  x <- toupper(
    as.character(x)
  )

  ifelse(
    grepl(
      "STAGE IV",
      x
    ),
    "Stage IV",

    ifelse(
      grepl(
        "STAGE III",
        x
      ),
      "Stage III",

      ifelse(
        grepl(
          "STAGE II",
          x
        ),
        "Stage II",

        ifelse(
          grepl(
            "STAGE I",
            x
          ),
          "Stage I",
          NA
        )
      )
    )
  )
}

# ==============================
# Locate Salmon quantifications
# ==============================

samp <- list.dirs(
  salmon_dir,
  recursive = FALSE,
  full.names = TRUE
)

samp <- samp[
  file.exists(
    file.path(
      samp,
      "quant.sf"
    )
  )
]

files <- file.path(
  samp,
  "quant.sf"
)

names(files) <- basename(
  samp
)

# ==============================
# Parse transcript annotation
# ==============================

txr <- fread(
  cmd = paste(
    "awk '$3 == \"transcript\"'",
    shQuote(gff_file)
  ),
  sep = "\t",
  header = FALSE,
  quote = "",
  select = 9L
)

setnames(
  txr,
  "V9"
)

txr[, tid := sub(
  ".*transcript_id=([^;]+).*",
  "\\1",
  V9
)]

txr[, gid := sub(
  ".*gene_id=([^;]+).*",
  "\\1",
  V9
)]

# ==============================
# Restrict to target lncRNAs
# ==============================

lnc <- fread(
  lnc_file
)

targets <- strip(
  unique(
    lnc$ens_id
  )
)

ex <- fread(
  files[1],
  header = TRUE
)

g2 <- unique(
  txr[, .(
    k = strip(tid),
    gid
  )]
)

m <- g2[
  match(
    strip(ex$Name),
    k
  )
]

tx2gene <- data.frame(
  tx = ex$Name,
  gene = m$gid
)

tx2gene <- tx2gene[
  !is.na(tx2gene$gene) &
    strip(tx2gene$gene) %in% targets,
]

# ==============================
# tximport
# ==============================

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)

# ==============================
# Clinical stage information
# ==============================

clin <- fread(
  clin_file
)

clin[, patient_id := submitter_id]

clin[, stage := clean_stage(
  ajcc_stage
)]

clin <- clin[
  !is.na(stage)
]

clin[, stage := factor(
  stage,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)]

# ==============================
# Match TCGA samples to patients
# ==============================

sample_table <- data.table(
  sample = names(files)
)

sample_table[
  ,
  patient_id := substr(
    sample,
    1,
    12
  )
]

sample_table <- merge(
  sample_table,
  clin[
    ,
    .(
      patient_id,
      stage
    )
  ],
  by = "patient_id"
)

sample_table <- sample_table[
  !duplicated(
    patient_id
  )
]

keep_samples <- sample_table$sample

txi$counts <- txi$counts[
  ,
  keep_samples,
  drop = FALSE
]

txi$abundance <- txi$abundance[
  ,
  keep_samples,
  drop = FALSE
]

txi$length <- txi$length[
  ,
  keep_samples,
  drop = FALSE
]

rownames(
  sample_table
) <- sample_table$sample

# ==============================
# DESeq2 LRT
# ==============================

dds <- DESeqDataSetFromTximport(
  txi,
  colData = as.data.frame(
    sample_table
  ),
  design = ~ stage
)

dds <- dds[
  rowSums(
    counts(dds)
  ) >= 10,
]

dds <- DESeq(
  dds,
  test = "LRT",
  reduced = ~ 1
)

# ==============================
# Results
# ==============================

res <- as.data.frame(
  results(dds)
)

res$ens_id <- rownames(
  res
)

res$ens_id_stripped <- strip(
  res$ens_id
)

gene_map <- unique(
  lnc[
    ,
    .(
      ens_id,
      gene_name
    )
  ]
)

gene_map[
  ,
  ens_id_stripped := strip(
    ens_id
  )
]

res <- merge(
  res,
  gene_map,
  by = "ens_id_stripped",
  all.x = TRUE
)

res <- res[
  order(
    res$padj
  ),
]

# ==============================
# Save outputs
# ==============================

fwrite(
  res,
  file.path(
    out_dir,
    "TCGA_stage_LRT_DESeq2_results.csv"
  )
)

saveRDS(
  dds,
  file.path(
    out_dir,
    "TCGA_stage_LRT_dds.rds"
  )
)

cat(
  "Samples used:",
  ncol(dds),
  "\n"
)

cat(
  "Genes tested:",
  nrow(res),
  "\n"
)

cat(
  "padj < 0.05:",
  sum(
    res$padj < 0.05,
    na.rm = TRUE
  ),
  "\n"
)
