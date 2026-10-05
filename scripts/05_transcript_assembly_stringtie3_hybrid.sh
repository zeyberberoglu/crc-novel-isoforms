#!/bin/bash
#SBATCH --job-name=stringtie3_mix
#SBATCH --output=slurm-%j-stringtie3_mix.out
#SBATCH --error=slurm-%j-stringtie3_mix.err
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
#SBATCH --time=120:00:00

set -euo pipefail

# ==============================
# Conda environment
# ==============================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate stringtie_env

# ==============================
# User-defined paths
# ==============================

BASE="/path/to/chitebuster_runs"
LONG_BAMDIR="/path/to/long_read_bams"
MAP="/path/to/sample_to_srr.tsv"
OUTDIR="/path/to/stringtie3_hybrid_output"

REFERENCE_DIR="/path/to/reference"
GTF="${REFERENCE_DIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3"

mkdir -p "$OUTDIR"

# ==============================
# Run hybrid StringTie3
# ==============================

while read -r sample srr
do
    short_bam="${BASE}/chitebuster_output_${sample}/${sample}/2_star/${sample}_Aligned.sortedByCoord.out.bam"

    long_bam="${LONG_BAMDIR}/${srr}.sorted.bam"

    if [[ ! -f "$long_bam" ]]; then
        echo "Long-read BAM missing for ${sample}, skipping..."
        echo "Expected: ${long_bam}"
        continue
    fi

    if [[ ! -f "$short_bam" ]]; then
        echo "Short-read BAM missing for ${sample}, skipping..."
        echo "Expected: ${short_bam}"
        continue
    fi

    echo "Running hybrid StringTie3 for ${sample}"
    echo "Short BAM: ${short_bam}"
    echo "Long BAM: ${long_bam}"

    stringtie --mix \
        "$short_bam" \
        "$long_bam" \
        -G "$GTF" \
        -o "${OUTDIR}/${sample}.hybrid.gtf" \
        -A "${OUTDIR}/${sample}_gene_abund.tab" \
        -p 8 \
        -f 0.05 \
        -j 2 \
        -c 2 \
        -s 5 \
        --rf \
        -N \
        -l "${sample}"

done < "$MAP"

echo "All hybrid StringTie3 runs completed."
