suppressMessages({
  library(tximport)
  library(data.table)
  library(jsonlite)
})

# ==============================
# User-defined paths
# ==============================

salmon_dir <- "/path/to/TCGA/salmon_quant"
lnc_file <- "/path/to/robust_upregulated_lncRNAs_for_survival.csv"
gff_file <- "/path/to/final_extended_annotation.gff"
clin_json <- "/path/to/clinical_data.json"
out_dir <- "/path/to/TCGA/survival"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

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
    file.path(samp, "quant.sf")
  )
]

files <- file.path(
  samp,
  "quant.sf"
)

names(files) <- basename(samp)

# ==============================
# Build transcript-to-gene map
# ==============================

ex <- fread(
  files[1],
  header = TRUE
)

tmp <- tempfile()

system2(
  "grep",
  c(
    "-P",
    "\\ttranscript\\t",
    gff_file
  ),
  stdout = tmp
)

txr <- fread(
  tmp,
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

strip <- function(x) {
  ifelse(
    grepl("^ENS", x),
    sub("\\.[0-9]+$", "", x),
    x
  )
}

# ==============================
# Restrict to target lncRNAs
# ==============================

lnc <- fread(lnc_file)

targets <- strip(
  unique(lnc$ens_id)
)

keep_tx <- txr[
  strip(gid) %in% targets,
  tid
]

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
# Import TPM values
# ==============================

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)

tpm <- as.data.frame(
  txi$abundance
)

resolve <- function(g) {

  if (g %in% rownames(tpm)) {
    return(g)
  }

  h <- rownames(tpm)[
    strip(rownames(tpm)) == strip(g)
  ]

  if (length(h)) {
    h[1]
  } else {
    NA_character_
  }
}

map <- data.table(
  ens_id = unique(lnc$ens_id)
)

map[, matched := vapply(
  ens_id,
  resolve,
  ""
)]

fwrite(
  map,
  file.path(
    out_dir,
    "gene_id_mapping.csv"
  )
)

keep <- map[
  !is.na(matched)
]

expr <- tpm[
  keep$matched,
  ,
  drop = FALSE
]

rownames(expr) <- keep$ens_id

expr <- data.frame(
  ens_id = rownames(expr),
  gene_name = lnc$gene_name[
    match(
      rownames(expr),
      lnc$ens_id
    )
  ],
  expr,
  check.names = FALSE
)

fwrite(
  expr,
  file.path(
    out_dir,
    "lncRNA_TPM_matrix.csv"
  )
)

# ==============================
# Parse TCGA clinical data
# ==============================

cj <- fromJSON(
  clin_json,
  simplifyVector = FALSE
)

clin <- rbindlist(
  lapply(
    cj,
    function(x) {

      d <- x$diagnoses[[1]]
      dm <- x$demographic

      data.table(
        submitter_id = x$submitter_id,
        project_id = x$project$project_id,
        vital_status = dm$vital_status,
        days_to_death =
          if (!is.null(dm$days_to_death))
            dm$days_to_death
          else
            NA_real_,
        days_to_last_follow_up =
          if (!is.null(d$days_to_last_follow_up))
            d$days_to_last_follow_up
          else
            NA_real_,
        age_at_diagnosis =
          if (!is.null(d$age_at_diagnosis))
            d$age_at_diagnosis
          else
            NA_real_,
        gender = dm$gender,
        ajcc_stage =
          if (!is.null(d$ajcc_pathologic_stage))
            d$ajcc_pathologic_stage
          else
            NA_character_
      )
    }
  ),
  fill = TRUE
)

fwrite(
  clin,
  file.path(
    out_dir,
    "clinical_parsed.csv"
  )
)

cat(
  "samples:",
  length(files),
  "| lncRNAs matched:",
  nrow(keep),
  "| NA:",
  sum(is.na(map$matched)),
  "\n"
)
