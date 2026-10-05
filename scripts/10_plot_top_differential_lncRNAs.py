import pandas as pd
import matplotlib.pyplot as plt

# Load DE results
df = pd.read_csv("DESeq2_lncRNA_shrunken.csv")

# Remove rows without adjusted p-values
df = df.dropna(subset=["padj"])

# Keep significant lncRNAs
sig = df[
    (df["padj"] < 0.05)
    & (abs(df["log2FoldChange"]) > 1)
].copy()

# Top 10 upregulated
up = sig.sort_values(
    "log2FoldChange",
    ascending=False
).head(10)

# Top 10 downregulated
down = sig.sort_values(
    "log2FoldChange"
).head(10)

plot_df = pd.concat(
    [down, up]
)

# Use gene name when available
plot_df["label"] = (
    plot_df["gene_name"]
    .fillna(plot_df["gene_id"])
)

plt.figure(
    figsize=(10, 8)
)

plt.barh(
    plot_df["label"],
    plot_df["log2FoldChange"]
)

plt.axvline(
    0,
    color="black"
)

plt.xlabel(
    "Log2 Fold Change (CRC vs NAT)"
)

plt.ylabel(
    "lncRNA"
)

plt.title(
    "Top differentially expressed lncRNAs"
)

plt.tight_layout()

plt.savefig(
    "Top_DE_lncRNAs_barplot.png",
    dpi=300
)

plt.savefig(
    "Top_DE_lncRNAs_barplot.pdf"
)

print(
    "Saved Top_DE_lncRNAs_barplot.png"
)
