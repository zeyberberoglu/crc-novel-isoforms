#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import kruskal
except ImportError:
    kruskal = None


def attrs(text):
    out = {}
    for key, value in re.findall(r'([\w.-]+)\s+"([^"]+)"', text):
        out[key] = value
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if "=" in field:
            key, value = field.split("=", 1)
            out.setdefault(key.strip(), value.strip())
    return out


def strip_version(x):
    return re.sub(r"\.\d+$", "", str(x).strip())


def read_mapping(annotation):
    rows = []
    with open(annotation) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 9 or f[2].lower() not in {"transcript", "mrna"}:
                continue
            a = attrs(f[8])
            tx = a.get("transcript_id") or a.get("ID")
            gene = a.get("gene_id") or a.get("Parent") or tx
            name = a.get("gene_name") or a.get("Name") or gene
            if tx:
                tx = re.sub(r"^transcript:", "", tx)
                gene = re.sub(r"^gene:", "", gene)
                rows.append((strip_version(tx), gene, strip_version(gene), name))
    mapping = pd.DataFrame(
        rows, columns=["tx_key", "gene_id", "gene_key", "gene_name"]
    ).drop_duplicates("tx_key")
    if mapping.empty:
        raise SystemExit("No transcript-to-gene mappings were found.")
    return mapping


def read_candidates(path):
    df = pd.read_csv(path, sep=None, engine="python")
    priorities = ["gene_name", "symbol", "candidate", "gene_id", "gene"]
    column = None
    for priority in priorities:
        for candidate_column in df.columns:
            if priority in str(candidate_column).lower():
                column = candidate_column
                break
        if column is not None:
            break
    if column is None:
        column = df.columns[0]
    return list(dict.fromkeys(
        df[column].dropna().astype(str).str.strip().tolist()
    ))


def sample_gene_tpm(quant_file, mapping):
    q = pd.read_csv(quant_file, sep="\t", usecols=["Name", "TPM"])
    q["tx_key"] = q["Name"].map(strip_version)
    q = q.merge(mapping, on="tx_key", how="left")
    q["gene_id"] = q["gene_id"].fillna(q["Name"])
    q["gene_key"] = q["gene_key"].fillna(q["Name"].map(strip_version))
    q["gene_name"] = q["gene_name"].fillna(q["gene_id"])
    return q.groupby(
        ["gene_id", "gene_key", "gene_name"], as_index=False
    )["TPM"].sum()


def collect(tissue, directory, mapping):
    rows = []
    quant_files = sorted(Path(directory).rglob("quant.sf"))
    if not quant_files:
        raise SystemExit(f"No quant.sf files found for {tissue}: {directory}")
    print(f"{tissue}: {len(quant_files)} quant.sf files")
    for quant_file in quant_files:
        sample = f"{tissue}_{quant_file.parent.name}"
        g = sample_gene_tpm(quant_file, mapping)
        g["Sample"] = sample
        g["Tissue"] = tissue
        rows.append(g)
    return pd.concat(rows, ignore_index=True)


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if not len(vals):
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    idx = np.where(ok)[0]
    out[idx[order]] = np.minimum(adj, 1)
    return out


def heatmap(matrix, output, zscore=False, title=""):
    plot = np.log2(matrix + 1)
    if zscore:
        sd = plot.std(axis=1).replace(0, np.nan)
        plot = plot.sub(plot.mean(axis=1), axis=0).div(sd, axis=0).fillna(0)

    fig, ax = plt.subplots(
        figsize=(max(8, matrix.shape[1] * 0.18),
                 max(7, matrix.shape[0] * 0.35))
    )
    image = ax.imshow(plot.to_numpy(), aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(plot.index)))
    ax.set_yticklabels(plot.index, fontsize=7)
    ax.set_xticks(range(len(plot.columns)))
    ax.set_xticklabels(plot.columns, rotation=90, fontsize=5)
    ax.set_title(title)
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--colon", required=True)
    parser.add_argument("--stomach", required=True)
    parser.add_argument("--kidney", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    mapping = read_mapping(args.annotation)
    candidates = read_candidates(args.candidates)

    all_expr = pd.concat([
        collect("Colon_mucosa", args.colon, mapping),
        collect("Stomach", args.stomach, mapping),
        collect("Kidney", args.kidney, mapping)
    ], ignore_index=True)

    all_expr["name_key"] = all_expr["gene_name"].astype(str).str.upper()
    all_expr["id_key"] = all_expr["gene_id"].astype(str).str.upper()
    all_expr["gene_key_upper"] = all_expr["gene_key"].astype(str).str.upper()

    matched = []
    report = []

    for candidate in candidates:
        key = candidate.upper()
        key_no_version = strip_version(candidate).upper()

        subset = all_expr[
            (all_expr["name_key"] == key)
            | (all_expr["id_key"] == key)
            | (all_expr["gene_key_upper"] == key_no_version)
        ].copy()

        report.append({
            "Candidate": candidate,
            "Matched": not subset.empty,
            "Matched_gene_names": (
                ";".join(sorted(subset["gene_name"].astype(str).unique()))
                if not subset.empty else ""
            ),
            "Matched_gene_ids": (
                ";".join(sorted(subset["gene_id"].astype(str).unique()))
                if not subset.empty else ""
            )
        })

        if not subset.empty:
            subset["Candidate"] = candidate
            matched.append(subset)

    report_df = pd.DataFrame(report)
    report_df.to_csv(
        output / "candidate_matching_report.tsv", sep="\t", index=False
    )

    if not matched:
        raise SystemExit("No candidates matched the annotation.")

    expression = pd.concat(matched, ignore_index=True)
    expression = expression.groupby(
        ["Candidate", "Sample", "Tissue"], as_index=False
    )["TPM"].sum()

    expression.to_csv(
        output / "candidate_expression_long.tsv", sep="\t", index=False
    )

    sample_matrix = expression.pivot_table(
        index="Candidate",
        columns="Sample",
        values="TPM",
        fill_value=0
    )
    sample_matrix.to_csv(
        output / "candidate_sample_TPM_matrix.tsv", sep="\t"
    )

    # Per-tissue metrics.
    summary = expression.groupby(["Candidate", "Tissue"]).agg(
        n_samples=("TPM", "size"),
        mean_TPM=("TPM", "mean"),
        median_TPM=("TPM", "median"),
        sd_TPM=("TPM", "std"),
        maximum_TPM=("TPM", "max"),
        samples_TPM_gt_0=("TPM", lambda s: int((s > 0).sum())),
        samples_TPM_ge_0_1=("TPM", lambda s: int((s >= 0.1).sum())),
        samples_TPM_ge_1=("TPM", lambda s: int((s >= 1).sum()))
    ).reset_index()

    summary["percent_samples_TPM_gt_0"] = (
        100 * summary["samples_TPM_gt_0"] / summary["n_samples"]
    )
    summary["percent_samples_TPM_ge_0_1"] = (
        100 * summary["samples_TPM_ge_0_1"] / summary["n_samples"]
    )
    summary["percent_samples_TPM_ge_1"] = (
        100 * summary["samples_TPM_ge_1"] / summary["n_samples"]
    )

    # CV is based on raw TPM: standard deviation divided by mean.
    summary["CV_raw_TPM"] = np.where(
        summary["mean_TPM"] > 0,
        summary["sd_TPM"] / summary["mean_TPM"],
        np.nan
    )

    summary.to_csv(
        output / "candidate_tissue_expression_summary.tsv",
        sep="\t", index=False
    )

    median_matrix = summary.pivot(
        index="Candidate", columns="Tissue", values="median_TPM"
    ).fillna(0)
    maximum_matrix = summary.pivot(
        index="Candidate", columns="Tissue", values="maximum_TPM"
    ).fillna(0)
    detection_matrix = summary.pivot(
        index="Candidate",
        columns="Tissue",
        values="percent_samples_TPM_ge_1"
    ).fillna(0)
    cv_matrix = summary.pivot(
        index="Candidate", columns="Tissue", values="CV_raw_TPM"
    )

    median_matrix.to_csv(
        output / "candidate_tissue_median_TPM_matrix.tsv", sep="\t"
    )
    maximum_matrix.to_csv(
        output / "candidate_tissue_maximum_TPM_matrix.tsv", sep="\t"
    )
    detection_matrix.to_csv(
        output / "candidate_tissue_percent_samples_TPM_ge_1_matrix.tsv",
        sep="\t"
    )
    cv_matrix.to_csv(
        output / "candidate_tissue_CV_matrix.tsv", sep="\t"
    )

    # Candidate-level tissue specificity:
    # fraction of total tissue-median expression contributed by the highest tissue.
    # Range: about 1/3 for equal expression across three tissues, to 1 for one-tissue-specific.
    specificity_rows = []

    for candidate in median_matrix.index:
        tissue_medians = median_matrix.loc[candidate]
        total_median = tissue_medians.sum()
        highest_tissue = tissue_medians.idxmax()
        highest_median = tissue_medians.max()

        specificity_score = (
            highest_median / total_median if total_median > 0 else np.nan
        )

        overall_values = expression.loc[
            expression["Candidate"] == candidate, "TPM"
        ]

        overall_mean = overall_values.mean()
        overall_sd = overall_values.std()
        overall_cv = (
            overall_sd / overall_mean if overall_mean > 0 else np.nan
        )

        specificity_rows.append({
            "Candidate": candidate,
            "Highest_expression_tissue": highest_tissue,
            "Highest_tissue_median_TPM": highest_median,
            "Sum_of_tissue_median_TPM": total_median,
            "Tissue_specificity_score": specificity_score,
            "Overall_median_TPM": overall_values.median(),
            "Overall_maximum_TPM": overall_values.max(),
            "Overall_percent_samples_TPM_gt_0": 100 * (overall_values > 0).mean(),
            "Overall_percent_samples_TPM_ge_0_1": 100 * (overall_values >= 0.1).mean(),
            "Overall_percent_samples_TPM_ge_1": 100 * (overall_values >= 1).mean(),
            "Overall_CV_raw_TPM": overall_cv
        })

    specificity = pd.DataFrame(specificity_rows)

    # Kruskal-Wallis across the three healthy tissues.
    stats = []
    for candidate, subset in expression.groupby("Candidate"):
        groups = [
            group["TPM"].to_numpy()
            for _, group in subset.groupby("Tissue")
        ]

        if kruskal is not None and len(groups) >= 2:
            try:
                statistic, p_value = kruskal(*groups)
            except ValueError:
                statistic, p_value = np.nan, np.nan
        else:
            statistic, p_value = np.nan, np.nan

        stats.append({
            "Candidate": candidate,
            "Kruskal_Wallis_statistic_raw_TPM": statistic,
            "Kruskal_Wallis_pvalue_raw_TPM": p_value
        })

    stats = pd.DataFrame(stats)
    stats["Kruskal_Wallis_padj_BH"] = bh(
        stats["Kruskal_Wallis_pvalue_raw_TPM"]
    )

    final_summary = specificity.merge(stats, on="Candidate", how="left")

    # Add per-tissue metrics in a wide, easy-to-read form.
    for tissue in ["Colon_mucosa", "Stomach", "Kidney"]:
        tissue_summary = summary[summary["Tissue"] == tissue].set_index("Candidate")
        for source, label in [
            ("median_TPM", "median_TPM"),
            ("maximum_TPM", "maximum_TPM"),
            ("percent_samples_TPM_gt_0", "percent_samples_TPM_gt_0"),
            ("percent_samples_TPM_ge_0_1", "percent_samples_TPM_ge_0_1"),
            ("percent_samples_TPM_ge_1", "percent_samples_TPM_ge_1"),
            ("CV_raw_TPM", "CV_raw_TPM")
        ]:
            final_summary[f"{tissue}_{label}"] = final_summary["Candidate"].map(
                tissue_summary[source]
            )

    final_summary = final_summary.sort_values(
        ["Overall_percent_samples_TPM_ge_1", "Overall_median_TPM"],
        ascending=[True, True]
    )

    final_summary.to_csv(
        output / "FINAL_candidate_healthy_tissue_summary.tsv",
        sep="\t", index=False
    )

    stats.to_csv(
        output / "candidate_cross_tissue_statistics.tsv",
        sep="\t", index=False
    )

    heatmap(
        sample_matrix,
        output / "sample_heatmap_log2TPM.png",
        False,
        "Healthy tissues: log2(TPM + 1)"
    )
    heatmap(
        sample_matrix,
        output / "sample_heatmap_gene_zscore.png",
        True,
        "Healthy tissues: gene-wise z-score of log2(TPM + 1)"
    )
    heatmap(
        median_matrix,
        output / "tissue_median_heatmap_log2TPM.png",
        False,
        "Tissue median expression: log2(TPM + 1)"
    )
    heatmap(
        maximum_matrix,
        output / "tissue_maximum_heatmap_log2TPM.png",
        False,
        "Tissue maximum expression: log2(TPM + 1)"
    )
    heatmap(
        detection_matrix,
        output / "tissue_detection_frequency_TPM_ge_1.png",
        False,
        "Detection frequency by tissue"
    )

    with open(output / "analysis_summary.txt", "w") as handle:
        handle.write(f"Candidates supplied: {len(candidates)}\n")
        handle.write(f"Candidates matched: {int(report_df['Matched'].sum())}\n")
        handle.write("\nSamples per tissue:\n")
        handle.write(
            expression.groupby("Tissue")["Sample"].nunique().to_string()
        )
        handle.write("\n\nMain final table:\n")
        handle.write("FINAL_candidate_healthy_tissue_summary.tsv\n")
        handle.write(
            "\nTissue specificity score: highest tissue median TPM divided "
            "by the sum of median TPM across all three tissues.\n"
        )
        handle.write(
            "CV: standard deviation of raw TPM divided by mean raw TPM.\n"
        )
        handle.write(
            "Kruskal-Wallis test: performed on raw TPM values.\n"
        )

    print(f"Finished. Results: {output}")


if __name__ == "__main__":
    main()
