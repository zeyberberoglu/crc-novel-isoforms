#!/usr/bin/env python3

import os
import re
import glob
import pandas as pd


BASE_DIR = "/path/to/PRJNA397552_fractionation"
METADATA_DIR = os.path.join(BASE_DIR, "metadata")
OUTPUT = os.path.join(METADATA_DIR, "fractionation_manifest.tsv")


def find_metadata_file():
    candidates = []

    for pattern in ("*.tsv", "*.csv", "*.txt"):
        candidates.extend(glob.glob(os.path.join(METADATA_DIR, pattern)))

    candidates = [
        path for path in candidates
        if os.path.abspath(path) != os.path.abspath(OUTPUT)
    ]

    for path in candidates:
        try:
            sep = "\t" if path.endswith((".tsv", ".txt")) else ","
            test = pd.read_csv(path, sep=sep, nrows=5)

            normalised = {
                col.lower().replace(" ", "_"): col
                for col in test.columns
            }

            has_run = any(
                name in normalised
                for name in ["run", "run_accession"]
            )

            has_title = any(
                name in normalised
                for name in [
                    "sample_title",
                    "experiment_title",
                    "sample_name",
                    "title"
                ]
            )

            if has_run and has_title:
                return path

        except Exception:
            continue

    raise FileNotFoundError(
        "Could not automatically locate the original ENA/SRA metadata "
        "file inside: {}".format(METADATA_DIR)
    )


def find_column(df, possible_names, required=True):
    lookup = {
        col.lower().strip().replace(" ", "_"): col
        for col in df.columns
    }

    for name in possible_names:
        key = name.lower().strip().replace(" ", "_")

        if key in lookup:
            return lookup[key]

    if required:
        raise ValueError(
            "Could not find any of these columns: {}".format(
                ", ".join(possible_names)
            )
        )

    return None


def clean(value):
    if pd.isna(value):
        return ""

    return str(value).strip()


def classify_sample(title):
    """
    Returns:
        Sample, LibraryType, CellLine, Fraction, Replicate

    Selected design:
      QuantSeq:
        1CT, SW480 and SW620
        C and N
        rep1 and rep2

      Full RNA-seq:
        1CT and SW620 only
        C and N
    """

    title = clean(title)

    # Exact untreated QuantSeq names:
    # 1CT_C_rep1
    # SW480_N_rep2
    quant_match = re.fullmatch(
        r"(1CT|SW480|SW620)_([CN])_rep([12])",
        title,
        flags=re.IGNORECASE
    )

    if quant_match:
        cell_line = quant_match.group(1).upper()

        if cell_line == "1CT":
            cell_line = "1CT"

        fraction = quant_match.group(2).upper()
        replicate = int(quant_match.group(3))

        sample = "{}_{}_rep{}_QuantSeq3prime".format(
            cell_line,
            fraction,
            replicate
        )

        return (
            sample,
            "QuantSeq3prime",
            cell_line,
            fraction,
            replicate
        )

    # Full RNA-seq is available only for 1CT and SW620.
    full_match = re.fullmatch(
        r"(1CT|SW620)_(Cytoplasmic|Nuclear)_FullRNASeq",
        title,
        flags=re.IGNORECASE
    )

    if full_match:
        cell_line = full_match.group(1).upper()

        if cell_line == "1CT":
            cell_line = "1CT"

        fraction_word = full_match.group(2).lower()
        fraction = "C" if fraction_word == "cytoplasmic" else "N"

        sample = "{}_{}_FullRNASeq".format(
            cell_line,
            fraction
        )

        return (
            sample,
            "FullRNASeq",
            cell_line,
            fraction,
            1
        )

    return None


os.makedirs(METADATA_DIR, exist_ok=True)

input_file = find_metadata_file()

print("Reading metadata from:")
print(input_file)

separator = "\t" if input_file.endswith((".tsv", ".txt")) else ","
metadata = pd.read_csv(input_file, sep=separator, dtype=str)

run_col = find_column(
    metadata,
    ["Run", "run_accession"]
)

sample_title_col = find_column(
    metadata,
    ["SampleTitle", "sample_title", "sample_name", "title"]
)

experiment_title_col = find_column(
    metadata,
    ["ExperimentTitle", "experiment_title"],
    required=False
)

layout_col = find_column(
    metadata,
    ["LibraryLayout", "library_layout"],
    required=False
)

ftp_col = find_column(
    metadata,
    ["FastqFTP", "fastq_ftp"],
    required=False
)

md5_col = find_column(
    metadata,
    ["FastqMD5", "fastq_md5"],
    required=False
)

selected_rows = []

for _, row in metadata.iterrows():
    run = clean(row[run_col])
    sample_title = clean(row[sample_title_col])

    classification = classify_sample(sample_title)

    if classification is None:
        continue

    (
        sample,
        library_type,
        cell_line,
        fraction,
        replicate
    ) = classification

    layout = (
        clean(row[layout_col]).upper()
        if layout_col is not None
        else "SINGLE"
    )

    if not layout:
        layout = "SINGLE"

    if layout != "SINGLE":
        print(
            "Skipping non-single-end sample: {} {}".format(
                run,
                sample_title
            )
        )
        continue

    fastq_ftp = (
        clean(row[ftp_col])
        if ftp_col is not None
        else ""
    )

    fastq_md5 = (
        clean(row[md5_col])
        if md5_col is not None
        else ""
    )

    experiment_title = (
        clean(row[experiment_title_col])
        if experiment_title_col is not None
        else ""
    )

    selected_rows.append({
        "Run": run,
        "Sample": sample,
        "LibraryType": library_type,
        "CellLine": cell_line,
        "Fraction": fraction,
        "Replicate": replicate,
        "Layout": layout,
        "FastqFTP": fastq_ftp,
        "FastqMD5": fastq_md5,
        "SampleTitle": sample_title,
        "ExperimentTitle": experiment_title
    })


manifest = pd.DataFrame(selected_rows)

if manifest.empty:
    raise RuntimeError("ERROR: No matching samples were found.")

manifest = manifest.drop_duplicates(
    subset=["Run", "Sample"]
)

library_order = {
    "FullRNASeq": 0,
    "QuantSeq3prime": 1
}

cell_order = {
    "1CT": 0,
    "SW480": 1,
    "SW620": 2
}

fraction_order = {
    "C": 0,
    "N": 1
}

manifest["_library_order"] = manifest["LibraryType"].map(library_order)
manifest["_cell_order"] = manifest["CellLine"].map(cell_order)
manifest["_fraction_order"] = manifest["Fraction"].map(fraction_order)

manifest = manifest.sort_values(
    by=[
        "_library_order",
        "_cell_order",
        "_fraction_order",
        "Replicate"
    ]
)

manifest = manifest.drop(
    columns=[
        "_library_order",
        "_cell_order",
        "_fraction_order"
    ]
)

manifest.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

print("\nSelected samples:\n")
print(manifest.to_string(index=False))

print("\nGroup counts:\n")
print(
    manifest.groupby(
        ["LibraryType", "CellLine", "Fraction"]
    ).size().to_string()
)

print("\nSummary:")
print(
    "Full RNA-seq samples:",
    sum(manifest["LibraryType"] == "FullRNASeq")
)
print(
    "QuantSeq 3-prime samples:",
    sum(manifest["LibraryType"] == "QuantSeq3prime")
)
print("Total samples:", len(manifest))
print("Manifest written to:", OUTPUT)
