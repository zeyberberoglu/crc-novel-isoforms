import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import (
    mannwhitneyu,
    kruskal
)

from statsmodels.stats.multitest import multipletests

# ==============================
# User-defined paths
# ==============================

OUT = "/path/to/TCGA/survival"

PLOT_DIR = os.path.join(
    OUT,
    "stage_TPM_plots"
)

os.makedirs(
    PLOT_DIR,
    exist_ok=True
)

# ==============================
# Load data
# ==============================

expr = pd.read_csv(
    os.path.join(
        OUT,
        "lncRNA_TPM_matrix.csv"
    )
)

clin = pd.read_csv(
    os.path.join(
        OUT,
        "clinical_parsed.csv"
    )
)

# ==============================
# Clean tumour stage
# ==============================

def clean_stage(x):

    if pd.isna(x):
        return np.nan

    x = str(x).upper()

    if "STAGE IV" in x:
        return "Stage IV"

    elif "STAGE III" in x:
        return "Stage III"

    elif "STAGE II" in x:
        return "Stage II"

    elif "STAGE I" in x:
        return "Stage I"

    return np.nan


def p_to_star(p):

    if p < 0.0001:
        return "****"

    elif p < 0.001:
        return "***"

    elif p < 0.01:
        return "**"

    elif p < 0.05:
        return "*"

    else:
        return "ns"


clin["stage_clean"] = (
    clin["ajcc_stage"]
    .apply(clean_stage)
)

clin = clin.dropna(
    subset=[
        "stage_clean"
    ]
)

clin = clin.set_index(
    "submitter_id"
)

# ==============================
# Prepare expression matrix
# ==============================

meta = expr[
    [
        "ens_id",
        "gene_name"
    ]
]

mat = (
    expr
    .set_index("ens_id")
    .drop(
        columns=["gene_name"]
    )
)

mat.columns = [
    c[:12]
    for c in mat.columns
]

mat = (
    mat.T
    .groupby(level=0)
    .mean()
    .T
)

common = [
    s
    for s in mat.columns
    if s in clin.index
]

mat = mat[
    common
]

clin = clin.loc[
    common
]

stage_order = [
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
]

summary_rows = []

# ==============================
# Analyze each lncRNA
# ==============================

for ens_id, row in mat.iterrows():

    gene_name = meta.loc[
        meta["ens_id"] == ens_id,
        "gene_name"
    ].iloc[0]

    df = pd.DataFrame(
        {
            "patient": mat.columns,
            "TPM": np.log2(
                row.values.astype(float)
                + 1
            ),
            "stage": clin[
                "stage_clean"
            ].values
        }
    ).dropna()

    df = df[
        df["stage"].isin(
            stage_order
        )
    ]

    groups = []
    labels = []

    for st in stage_order:

        vals = df.loc[
            df["stage"] == st,
            "TPM"
        ].values.astype(float)

        if len(vals) > 0:
            groups.append(vals)
            labels.append(st)

    if len(groups) < 2:
        continue

    try:
        overall_p = kruskal(
            *groups
        ).pvalue

    except Exception:
        overall_p = np.nan

    # ==============================
    # Pairwise comparisons
    # ==============================

    comparisons = []
    pvals = []

    for i in range(
        len(groups)
    ):

        for j in range(
            i + 1,
            len(groups)
        ):

            if (
                len(groups[i]) >= 3
                and len(groups[j]) >= 3
            ):

                p = mannwhitneyu(
                    groups[i],
                    groups[j],
                    alternative="two-sided"
                ).pvalue

                comparisons.append(
                    (
                        i,
                        j
                    )
                )

                pvals.append(p)

    if len(pvals) > 0:

        padj = multipletests(
            pvals,
            method="fdr_bh"
        )[1]

    else:

        padj = []

    # ==============================
    # Plot
    # ==============================

    fig, ax = plt.subplots(
        figsize=(5.5, 4.5),
        dpi=300
    )

    positions = np.arange(
        1,
        len(groups) + 1
    )

    ax.boxplot(
        groups,
        positions=positions,
        widths=0.55,
        showfliers=False,
        patch_artist=True
    )

    for pos, vals in zip(
        positions,
        groups
    ):

        jitter = np.random.normal(
            pos,
            0.05,
            size=len(vals)
        )

        ax.scatter(
            jitter,
            vals,
            s=14,
            alpha=0.65
        )

    ax.set_xticks(
        positions
    )

    ax.set_xticklabels(
        labels,
        rotation=30,
        ha="right"
    )

    ax.set_ylabel(
        "log2(TPM + 1)"
    )

    ax.set_xlabel(
        "Tumour stage"
    )

    ax.set_title(
        str(gene_name),
        fontstyle="italic"
    )

    ax.spines[
        [
            "top",
            "right"
        ]
    ].set_visible(
        False
    )

    ymax = max(
        [
            np.max(g)
            for g in groups
        ]
    )

    if ymax <= 0:
        ymax = 1

    y = ymax * 1.10
    step = ymax * 0.12

    for comp, q in zip(
        comparisons,
        padj
    ):

        if q >= 0.05:
            continue

        i, j = comp

        x1 = positions[i]
        x2 = positions[j]

        ax.plot(
            [
                x1,
                x1,
                x2,
                x2
            ],
            [
                y,
                y + step * 0.2,
                y + step * 0.2,
                y
            ],
            lw=1.2,
            c="black"
        )

        ax.text(
            (x1 + x2) / 2,
            y + step * 0.25,
            p_to_star(q),
            ha="center",
            va="bottom",
            fontsize=10
        )

        y += step

    ax.set_ylim(
        0,
        y + step
    )

    if not np.isnan(
        overall_p
    ):

        ax.text(
            0.98,
            0.98,
            f"Kruskal P = {overall_p:.2e}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9
        )

    plt.tight_layout()

    safe_gene = "".join(
        ch
        if ch.isalnum()
        else "_"
        for ch in str(gene_name)
    )

    fig.savefig(
        os.path.join(
            PLOT_DIR,
            f"{safe_gene}_TPM_by_stage.pdf"
        ),
        bbox_inches="tight"
    )

    fig.savefig(
        os.path.join(
            PLOT_DIR,
            f"{safe_gene}_TPM_by_stage.png"
        ),
        bbox_inches="tight"
    )

    plt.close(fig)

    summary_rows.append(
        {
            "ens_id": ens_id,
            "gene_name": gene_name,
            "n_patients": len(df),
            "kruskal_p": overall_p,
            "min_pairwise_FDR":
                np.min(padj)
                if len(padj) > 0
                else np.nan
        }
    )

# ==============================
# Save summary
# ==============================

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    os.path.join(
        OUT,
        "stage_TPM_plot_statistics.csv"
    ),
    index=False
)

print("Done")

print(
    "Patients matched:",
    len(common)
)

print(
    "Plots saved in:",
    PLOT_DIR
)

print(
    "Summary saved:",
    os.path.join(
        OUT,
        "stage_TPM_plot_statistics.csv"
    )
)
