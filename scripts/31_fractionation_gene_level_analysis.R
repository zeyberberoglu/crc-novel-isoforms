#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(DESeq2)
})

options(stringsAsFactors = FALSE)

base_dir   <- Sys.getenv("BASE_DIR")
salmon_dir <- Sys.getenv("SALMON_DIR")
gff_file   <- Sys.getenv("GFF")
out_dir    <- Sys.getenv("OUT_DIR")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots_pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots_png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Prioritised CRCN lncRNAs
###############################################################################

target_genes <- c(
    "CRCN_LINC00280",
    "CRCN_LINC00335",
    "CRCN_FANCC-AS1",
    "CRCN_LINC00128",
    "CRCN_DENND3-AS1",
    "CRCN_KANTR-OT1",
    "CRCN_LINC00205"
)

###############################################################################
# Helper functions
###############################################################################

extract_attribute <- function(attribute_string, key) {
    patterns <- c(
        paste0("(?:^|;)", key, "=([^;]+)"),
        paste0("(?:^|;)\\s*", key, "\\s+\"([^\"]+)\""),
        paste0("(?:^|;)", key, "\\s+([^;]+)")
    )

    result <- rep(NA_character_, length(attribute_string))

    for (pattern in patterns) {
        missing_rows <- which(is.na(result))

        if (length(missing_rows) == 0) {
            break
        }

        matches <- regexec(
            pattern,
            attribute_string[missing_rows],
            perl = TRUE
        )

        extracted <- regmatches(
            attribute_string[missing_rows],
            matches
        )

        values <- vapply(
            extracted,
            function(x) {
                if (length(x) >= 2) {
                    trimws(x[2])
                } else {
                    NA_character_
                }
            },
            character(1)
        )

        result[missing_rows[!is.na(values)]] <- values[!is.na(values)]
    }

    gsub('"', "", result, fixed = TRUE)
}

get_significance_label <- function(padj) {
    if (is.na(padj)) {
        return("")
    }

    if (padj < 0.0001) {
        return("****")
    }

    if (padj < 0.001) {
        return("***")
    }

    if (padj < 0.01) {
        return("**")
    }

    if (padj < 0.05) {
        return("*")
    }

    ""
}

safe_filename <- function(x) {
    gsub("[^A-Za-z0-9_.-]", "_", x)
}

###############################################################################
# Read transcript-to-gene mapping from the final extended annotation
###############################################################################

message("Reading annotation: ", gff_file)

if (!file.exists(gff_file)) {
    stop("Annotation file was not found: ", gff_file)
}

gff <- fread(
    gff_file,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    comment.char = "#",
    data.table = TRUE,
    showProgress = FALSE
)

if (ncol(gff) < 9) {
    stop("The annotation does not appear to contain nine GFF/GTF columns.")
}

setnames(
    gff,
    1:9,
    c(
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes"
    )
)

transcript_rows <- gff[
    feature %chin% c("transcript", "mRNA")
]

transcript_rows[, transcript_id := extract_attribute(attributes, "transcript_id")]

missing_tx <- is.na(transcript_rows$transcript_id)

if (any(missing_tx)) {
    transcript_rows$transcript_id[missing_tx] <-
        extract_attribute(
            transcript_rows$attributes[missing_tx],
            "ID"
        )
}

transcript_rows[, gene_id := extract_attribute(attributes, "gene_id")]
transcript_rows[, gene_name := extract_attribute(attributes, "gene_name")]

missing_gene_id <- is.na(transcript_rows$gene_id)

if (any(missing_gene_id)) {
    transcript_rows$gene_id[missing_gene_id] <-
        extract_attribute(
            transcript_rows$attributes[missing_gene_id],
            "Parent"
        )
}

transcript_rows[
    is.na(gene_name) | gene_name == "",
    gene_name := gene_id
]

tx2gene <- unique(
    transcript_rows[
        !is.na(transcript_id) &
        transcript_id != "" &
        !is.na(gene_name) &
        gene_name != "",
        .(
            transcript_id,
            gene_id,
            gene_name
        )
    ],
    by = "transcript_id"
)

message("Transcript-to-gene mappings: ", nrow(tx2gene))

found_targets <- intersect(target_genes, unique(tx2gene$gene_name))
missing_targets <- setdiff(target_genes, found_targets)

message(
    "Target genes found in annotation: ",
    paste(found_targets, collapse = ", ")
)

if (length(missing_targets) > 0) {
    warning(
        "These target gene names were not found in the annotation: ",
        paste(missing_targets, collapse = ", ")
    )
}

###############################################################################
# Locate and parse Salmon samples
###############################################################################

quant_files <- list.files(
    salmon_dir,
    pattern = "quant\\.sf$",
    recursive = TRUE,
    full.names = TRUE
)

if (length(quant_files) == 0) {
    stop("No quant.sf files were found under: ", salmon_dir)
}

sample_table <- data.table(
    quant_file = quant_files,
    sample = basename(dirname(quant_files))
)

sample_table[
    ,
    assay := fifelse(
        grepl(
            "QuantSeq|3prime|3_prime|3-prime",
            sample,
            ignore.case = TRUE
        ),
        "3prime",
        fifelse(
            grepl(
                "FullRNASeq|TotalRNA|Total_RNA",
                sample,
                ignore.case = TRUE
            ),
            "total",
            NA_character_
        )
    )
]

sample_table[
    ,
    fraction := fifelse(
        grepl("(^|_)C($|_)", sample),
        "Cytoplasmic",
        fifelse(
            grepl("(^|_)N($|_)", sample),
            "Nuclear",
            NA_character_
        )
    )
]

sample_table[
    ,
    cell_line := sub(
        "_[CN]($|_).*",
        "",
        sample
    )
]

sample_table[
    ,
    replicate := as.integer(
        sub(
            ".*_rep([0-9]+).*",
            "\\1",
            sample
        )
    )
]

sample_table[
    !grepl("_rep[0-9]+", sample),
    replicate := 1L
]

sample_table <- sample_table[
    !is.na(assay) &
    !is.na(fraction) &
    !is.na(cell_line)
]

sample_table[
    ,
    fraction := factor(
        fraction,
        levels = c("Cytoplasmic", "Nuclear")
    )
]

fwrite(
    sample_table,
    file.path(out_dir, "tables", "samples_used.csv")
)

message("Samples identified:")
print(
    sample_table[
        ,
        .(
            sample,
            cell_line,
            fraction,
            replicate,
            assay
        )
    ]
)

###############################################################################
# Import Salmon values and aggregate transcripts to genes
###############################################################################

all_expression <- vector("list", nrow(sample_table))

for (i in seq_len(nrow(sample_table))) {
    sample_name <- sample_table$sample[i]
    quant_file  <- sample_table$quant_file[i]

    message(
        "Importing ",
        i,
        "/",
        nrow(sample_table),
        ": ",
        sample_name
    )

    q <- fread(quant_file)

    required_columns <- c(
        "Name",
        "Length",
        "EffectiveLength",
        "TPM",
        "NumReads"
    )

    missing_columns <- setdiff(required_columns, names(q))

    if (length(missing_columns) > 0) {
        stop(
            "Missing columns in ",
            quant_file,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }

    setnames(q, "Name", "transcript_id")

    q <- merge(
        q,
        tx2gene,
        by = "transcript_id",
        all.x = FALSE,
        all.y = FALSE
    )

    gene_q <- q[
        ,
        .(
            TPM = sum(TPM, na.rm = TRUE),
            NumReads = sum(NumReads, na.rm = TRUE)
        ),
        by = .(
            gene_id,
            gene_name
        )
    ]

    gene_q[, sample := sample_name]
    gene_q[, cell_line := sample_table$cell_line[i]]
    gene_q[, fraction := as.character(sample_table$fraction[i])]
    gene_q[, replicate := sample_table$replicate[i]]
    gene_q[, assay := sample_table$assay[i]]

    all_expression[[i]] <- gene_q
}

expression_long <- rbindlist(
    all_expression,
    use.names = TRUE,
    fill = TRUE
)

expression_long[
    ,
    log2_TPM := log2(TPM + 1)
]

fwrite(
    expression_long,
    file.path(
        out_dir,
        "tables",
        "all_gene_fractionation_expression.csv"
    )
)

target_expression <- expression_long[
    gene_name %chin% target_genes
]

fwrite(
    target_expression,
    file.path(
        out_dir,
        "tables",
        "CRCN_target_fractionation_expression.csv"
    )
)

###############################################################################
# Paired DESeq2 analysis for 3-prime RNA-seq
#
# Analysis is performed separately for each cell line using:
#
# design = ~ replicate + fraction
#
# Nuclear is compared against Cytoplasmic.
###############################################################################

three_prime_metadata <- unique(
    sample_table[
        assay == "3prime",
        .(
            sample,
            cell_line,
            fraction,
            replicate
        )
    ]
)

stats_list <- list()
stats_index <- 1L

# Store DESeq2 size-factor-normalised counts for plotting
normalised_count_list <- list()
normalised_count_index <- 1L

for (current_cell_line in unique(three_prime_metadata$cell_line)) {
    message("Running paired DESeq2 for: ", current_cell_line)

    md <- copy(
        three_prime_metadata[
            cell_line == current_cell_line
        ]
    )

    paired_replicates <- md[
        ,
        .(
            n_fractions = uniqueN(fraction)
        ),
        by = replicate
    ][
        n_fractions == 2,
        replicate
    ]

    md <- md[
        replicate %in% paired_replicates
    ]

    if (length(unique(md$replicate)) < 2) {
        warning(
            "Skipping DESeq2 for ",
            current_cell_line,
            ": fewer than two complete paired replicates."
        )

        for (gene in target_genes) {
            stats_list[[stats_index]] <- data.table(
                gene_name = gene,
                cell_line = current_cell_line,
                n_pairs = length(unique(md$replicate)),
                baseMean = NA_real_,
                log2FoldChange = NA_real_,
                pvalue = NA_real_,
                padj = NA_real_,
                mean_cytoplasmic_TPM = NA_real_,
                mean_nuclear_TPM = NA_real_,
                all_TPM_zero = NA,
                significance = ""
            )

            stats_index <- stats_index + 1L
        }

        next
    }

    md[
        ,
        fraction := factor(
            fraction,
            levels = c("Cytoplasmic", "Nuclear")
        )
    ]

    md[
        ,
        replicate := factor(replicate)
    ]

    counts_long <- expression_long[
        assay == "3prime" &
        cell_line == current_cell_line &
        sample %chin% md$sample,
        .(
            count = sum(NumReads, na.rm = TRUE)
        ),
        by = .(
            gene_name,
            sample
        )
    ]

    count_matrix <- dcast(
        counts_long,
        gene_name ~ sample,
        value.var = "count",
        fill = 0
    )

    gene_names <- count_matrix$gene_name
    count_matrix[, gene_name := NULL]

    count_matrix <- as.matrix(count_matrix)
    rownames(count_matrix) <- gene_names

    count_matrix <- round(count_matrix)

    md <- md[
        match(colnames(count_matrix), sample)
    ]

    if (!identical(md$sample, colnames(count_matrix))) {
        stop(
            "Metadata and count matrix sample order could not be matched for ",
            current_cell_line
        )
    }

    coldata <- data.frame(
        row.names = md$sample,
        replicate = md$replicate,
        fraction = md$fraction
    )

    keep_genes <- rowSums(count_matrix) > 0
    count_matrix_filtered <- count_matrix[keep_genes, , drop = FALSE]

    if (nrow(count_matrix_filtered) < 2) {
        warning(
            "Skipping DESeq2 for ",
            current_cell_line,
            ": insufficient expressed genes."
        )
        next
    }

    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix_filtered,
        colData = coldata,
        design = ~ replicate + fraction
    )

    dds <- DESeq(
        dds,
        quiet = TRUE
    )

    ###########################################################################
    # Extract DESeq2 size-factor-normalised counts for visualisation
    ###########################################################################

    norm_matrix <- counts(
        dds,
        normalized = TRUE
    )

    norm_dt <- as.data.table(
        as.data.frame(norm_matrix),
        keep.rownames = "gene_name"
    )

    norm_long <- melt(
        norm_dt,
        id.vars = "gene_name",
        variable.name = "sample",
        value.name = "DESeq2_normalised_count"
    )

    norm_long[
        ,
        sample := as.character(sample)
    ]

    norm_long <- merge(
        norm_long,
        md[
            ,
            .(
                sample,
                cell_line,
                fraction,
                replicate
            )
        ],
        by = "sample",
        all.x = TRUE
    )

    normalised_count_list[[normalised_count_index]] <- norm_long
    normalised_count_index <- normalised_count_index + 1L

    result <- results(
        dds,
        contrast = c(
            "fraction",
            "Nuclear",
            "Cytoplasmic"
        ),
        alpha = 0.05
    )

    result_table <- as.data.table(
        as.data.frame(result),
        keep.rownames = "gene_name"
    )

    for (gene in target_genes) {
        target_values <- expression_long[
            assay == "3prime" &
            cell_line == current_cell_line &
            sample %chin% md$sample &
            gene_name == gene
        ]

        mean_cyto <- target_values[
            fraction == "Cytoplasmic",
            mean(TPM, na.rm = TRUE)
        ]

        mean_nuclear <- target_values[
            fraction == "Nuclear",
            mean(TPM, na.rm = TRUE)
        ]

        if (is.nan(mean_cyto)) {
            mean_cyto <- 0
        }

        if (is.nan(mean_nuclear)) {
            mean_nuclear <- 0
        }

        all_zero <- nrow(target_values) == 0 ||
            all(target_values$TPM <= 0, na.rm = TRUE)

        gene_result <- result_table[
            gene_name == gene
        ]

        if (nrow(gene_result) == 0) {
            gene_result <- data.table(
                baseMean = NA_real_,
                log2FoldChange = NA_real_,
                pvalue = NA_real_,
                padj = NA_real_
            )
        }

        final_padj <- gene_result$padj[1]

        # Critical safeguard:
        # never display significance if every TPM value is zero.
        if (all_zero) {
            final_padj <- NA_real_
        }

        significance <- get_significance_label(final_padj)

        stats_list[[stats_index]] <- data.table(
            gene_name = gene,
            cell_line = current_cell_line,
            n_pairs = length(unique(md$replicate)),
            baseMean = gene_result$baseMean[1],
            log2FoldChange = gene_result$log2FoldChange[1],
            pvalue = gene_result$pvalue[1],
            padj = final_padj,
            mean_cytoplasmic_TPM = mean_cyto,
            mean_nuclear_TPM = mean_nuclear,
            all_TPM_zero = all_zero,
            significance = significance
        )

        stats_index <- stats_index + 1L
    }
}

###############################################################################
# Combine DESeq2-normalised counts
###############################################################################

normalised_counts <- rbindlist(
    normalised_count_list,
    use.names = TRUE,
    fill = TRUE
)

fwrite(
    normalised_counts,
    file.path(
        out_dir,
        "tables",
        "CRCN_3prime_DESeq2_normalised_counts.csv"
    )
)

stats_table <- rbindlist(
    stats_list,
    use.names = TRUE,
    fill = TRUE
)

stats_table[
    ,
    localisation := fifelse(
        !is.na(log2FoldChange) & log2FoldChange > 0,
        "Nuclear enriched",
        fifelse(
            !is.na(log2FoldChange) & log2FoldChange < 0,
            "Cytoplasmic enriched",
            "Undetermined"
        )
    )
]

fwrite(
    stats_table,
    file.path(
        out_dir,
        "tables",
        "CRCN_3prime_paired_DESeq2_statistics.csv"
    )
)

###############################################################################
# Plotting
###############################################################################

plot_data <- copy(target_expression)

###############################################################################
# Add DESeq2-normalised counts to the 3-prime samples
###############################################################################

target_norm_counts <- normalised_counts[
    gene_name %chin% target_genes,
    .(
        gene_name,
        sample,
        DESeq2_normalised_count
    )
]

plot_data <- merge(
    plot_data,
    target_norm_counts,
    by = c(
        "gene_name",
        "sample"
    ),
    all.x = TRUE
)

###############################################################################
# Expression value used for plotting:
# 3-prime RNA-seq = DESeq2 size-factor-normalised counts
# Total RNA-seq   = TPM
###############################################################################

plot_data[
    ,
    plot_value := fifelse(
        assay == "3prime",
        log2(DESeq2_normalised_count + 1),
        log2(TPM + 1)
    )
]

plot_data[
    ,
    Assay := fifelse(
        assay == "3prime",
        "3\u2032 RNA-seq\n(DESeq2-normalised counts)",
        "Total RNA-seq\n(TPM)"
    )
]

plot_data[
    ,
    Assay := factor(
        Assay,
        levels = c(
            "3\u2032 RNA-seq\n(DESeq2-normalised counts)",
            "Total RNA-seq\n(TPM)"
        )
    )
]

plot_data[
    ,
    fraction := factor(
        fraction,
        levels = c(
            "Cytoplasmic",
            "Nuclear"
        )
    )
]

cell_line_order <- c(
    "1CT",
    "SW480",
    "SW620"
)

remaining_cell_lines <- setdiff(
    unique(plot_data$cell_line),
    cell_line_order
)

plot_data[
    ,
    cell_line := factor(
        cell_line,
        levels = c(
            intersect(cell_line_order, unique(cell_line)),
            sort(remaining_cell_lines)
        )
    )
]

fraction_colours <- c(
    "Cytoplasmic" = "#2C7BB6",
    "Nuclear" = "#D7191C"
)

fraction_shapes <- c(
    "Cytoplasmic" = 16,
    "Nuclear" = 17
)

summary_rows <- list()

for (gene in target_genes) {
    gene_data <- plot_data[
        gene_name == gene
    ]

    if (nrow(gene_data) == 0) {
        warning(
            "No quantified expression values were found for ",
            gene
        )
        next
    }

    gene_stats <- stats_table[
        gene_name == gene
    ]

    overall_direction <- "Localisation undetermined"

    significant_stats <- gene_stats[
        significance != ""
    ]

    if (nrow(significant_stats) > 0) {
        median_lfc <- median(
            significant_stats$log2FoldChange,
            na.rm = TRUE
        )

        if (is.finite(median_lfc)) {
            if (median_lfc > 0) {
                overall_direction <- "Nuclear enriched"
            } else if (median_lfc < 0) {
                overall_direction <- "Cytoplasmic enriched"
            }
        }
    } else {
        valid_lfc <- gene_stats[
            !is.na(log2FoldChange),
            log2FoldChange
        ]

        if (length(valid_lfc) > 0) {
            median_lfc <- median(
                valid_lfc,
                na.rm = TRUE
            )

            if (is.finite(median_lfc)) {
                if (median_lfc > 0) {
                    overall_direction <- "Nuclear enriched"
                } else if (median_lfc < 0) {
                    overall_direction <- "Cytoplasmic enriched"
                }
            }
        }
    }

    maximum_y <- max(
        gene_data[
            assay == "3prime",
            plot_value
        ],
        na.rm = TRUE
    )

    if (!is.finite(maximum_y)) {
        maximum_y <- 0
    }

    y_padding <- max(
        0.35,
        maximum_y * 0.12
    )

    bracket_y <- maximum_y + y_padding
    star_y <- bracket_y + y_padding * 0.22
    plot_upper_limit <- star_y + y_padding * 0.65

    annotation_data <- gene_stats[
        significance != "" &
        !is.na(padj) &
        !all_TPM_zero,
        .(
            cell_line,
            significance,
            padj
        )
    ]

    annotation_data[
        ,
        cell_line := factor(
            cell_line,
            levels = levels(gene_data$cell_line)
        )
    ]

    annotation_data[
        ,
        x_numeric := as.numeric(cell_line)
    ]

    annotation_data[
        ,
        `:=`(
            xmin = x_numeric - 0.17,
            xmax = x_numeric + 0.17,
            y = bracket_y,
            y_star = star_y
        )
    ]

    p <- ggplot(
        gene_data,
        aes(
            x = cell_line,
            y = plot_value,
            colour = fraction,
            shape = fraction
        )
    ) +
        geom_point(
            position = position_jitter(
                width = 0.055,
                height = 0
            ),
            size = 3.6,
            alpha = 0.95
        ) +
        stat_summary(
            aes(group = fraction),
            fun = mean,
            geom = "crossbar",
            width = 0.24,
            linewidth = 0.75,
            colour = "black",
            fatten = 0,
            position = position_dodge(width = 0.32)
        ) +
        facet_wrap(
            ~ Assay,
            scales = "free",
            nrow = 1
        ) +
        scale_colour_manual(
            values = fraction_colours,
            drop = FALSE
        ) +
        scale_shape_manual(
            values = fraction_shapes,
            drop = FALSE
        ) +
        scale_y_continuous(
            expand = expansion(
                mult = c(0.02, 0.18)
            )
        ) +
        labs(
            title = gene,
            subtitle = overall_direction,
            x = NULL,
            y = expression(
                log[2] * "(expression + 1)"
            ),
            colour = "Fraction",
            shape = "Fraction",
            caption = paste0(
                "3\u2032 RNA-seq expression is shown as log2(DESeq2 size-factor-normalised count + 1); ",
                "stars indicate paired DESeq2 analysis ",
                "(design: ~ replicate + fraction; BH-adjusted P values). ",
                "Total RNA-seq expression is shown as log2(TPM + 1) and is descriptive ",
                "because n = 1 per fraction and cell line."
            )
        ) +
        theme_classic(base_size = 15) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 22,
                hjust = 0.5
            ),
            plot.subtitle = element_text(
                size = 15,
                hjust = 0.5,
                margin = margin(
                    b = 12
                )
            ),
            strip.background = element_blank(),
            strip.text = element_text(
                face = "bold",
                size = 16,
                margin = margin(
                    b = 8
                )
            ),
            axis.title.y = element_text(
                size = 15
            ),
            axis.text.x = element_text(
                face = "bold",
                size = 13
            ),
            axis.text.y = element_text(
                size = 12
            ),
            legend.position = "top",
            legend.title = element_text(
                face = "bold"
            ),
            legend.text = element_text(
                size = 12
            ),
            plot.caption = element_text(
                hjust = 0,
                size = 9.5,
                margin = margin(
                    t = 12
                )
            ),
            panel.spacing.x = unit(
                1.2,
                "cm"
            ),
            plot.margin = margin(
                15,
                20,
                15,
                15
            )
        )

    # Add brackets and stars only to the 3-prime RNA-seq facet.
    if (nrow(annotation_data) > 0) {
        annotation_data[, Assay := factor(
            "3\u2032 RNA-seq",
            levels = levels(gene_data$Assay)
        )]

        p <- p +
            geom_segment(
                data = annotation_data,
                aes(
                    x = xmin,
                    xend = xmax,
                    y = y,
                    yend = y
                ),
                inherit.aes = FALSE,
                colour = "black",
                linewidth = 0.7
            ) +
            geom_segment(
                data = annotation_data,
                aes(
                    x = xmin,
                    xend = xmin,
                    y = y,
                    yend = y - y_padding * 0.12
                ),
                inherit.aes = FALSE,
                colour = "black",
                linewidth = 0.7
            ) +
            geom_segment(
                data = annotation_data,
                aes(
                    x = xmax,
                    xend = xmax,
                    y = y,
                    yend = y - y_padding * 0.12
                ),
                inherit.aes = FALSE,
                colour = "black",
                linewidth = 0.7
            ) +
            geom_text(
                data = annotation_data,
                aes(
                    x = x_numeric,
                    y = y_star,
                    label = significance
                ),
                inherit.aes = FALSE,
                colour = "black",
                fontface = "bold",
                size = 5.5
            )
    }

    file_stub <- safe_filename(gene)

    ggsave(
        filename = file.path(
            out_dir,
            "plots_pdf",
            paste0(
                file_stub,
                "_fractionation_v2.pdf"
            )
        ),
        plot = p,
        width = 13.5,
        height = 7.5,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        filename = file.path(
            out_dir,
            "plots_png",
            paste0(
                file_stub,
                "_fractionation_v2.png"
            )
        ),
        plot = p,
        width = 13.5,
        height = 7.5,
        units = "in",
        dpi = 350,
        bg = "white"
    )

    summary_rows[[gene]] <- data.table(
        gene_name = gene,
        plotted = TRUE,
        localisation_summary = overall_direction,
        maximum_log2_TPM = maximum_y,
        significant_cell_lines = paste(
            annotation_data$cell_line,
            collapse = ";"
        )
    )

    message("Created plots for: ", gene)
}

plot_summary <- rbindlist(
    summary_rows,
    use.names = TRUE,
    fill = TRUE
)

fwrite(
    plot_summary,
    file.path(
        out_dir,
        "tables",
        "CRCN_fractionation_plot_summary.csv"
    )
)

###############################################################################
# Print CRCN_LINC00280 specifically
###############################################################################

message("")
message("CRCN_LINC00280 expression values:")

print(
    target_expression[
        gene_name == "CRCN_LINC00280",
        .(
            sample,
            cell_line,
            fraction,
            replicate,
            assay,
            TPM,
            NumReads,
            log2_TPM
        )
    ][
        order(
            assay,
            cell_line,
            replicate,
            fraction
        )
    ]
)

message("")
message("CRCN_LINC00280 paired DESeq2 statistics:")

print(
    stats_table[
        gene_name == "CRCN_LINC00280"
    ]
)

message("")
message("Analysis complete.")
message("Results directory: ", out_dir)
