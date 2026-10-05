#!/bin/bash
#SBATCH --job-name=sqanti3_hybrid
#SBATCH --output=slurm-%j-sqanti3_hybrid.out
#SBATCH --error=slurm-%j-sqanti3_hybrid.err
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
#SBATCH --time=120:00:00

set -euo pipefail

# ==============================
# Conda environment
# ==============================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sqanti3

# ==============================
# User-defined paths
# ==============================

INDIR="/path/to/stringtie3_hybrid_output"
OUTDIR="/path/to/sqanti3_hybrid_output"

REFERENCE_DIR="/path/to/reference"
REF_GTF="${REFERENCE_DIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3"
REF_FASTA="${REFERENCE_DIR}/GRCh38.p14.genome.fa"

SHORT_BASE="/path/to/chitebuster_runs"

SQANTI3_DIR="/path/to/SQANTI3"

mkdir -p "$OUTDIR"

# ==============================
# Run SQANTI3 QC
# ==============================

for gtf in "${INDIR}"/*.hybrid.gtf
do
    sample=$(basename "$gtf" .hybrid.gtf)

    sample_out="${OUTDIR}/${sample}"
    mkdir -p "$sample_out"

    short_r1="${SHORT_BASE}/chitebuster_output_${sample}/${sample}/1_fastp/${sample}_reads1.fastq.gz"
    short_r2="${SHORT_BASE}/chitebuster_output_${sample}/${sample}/1_fastp/${sample}_reads2.fastq.gz"

    echo "Running SQANTI3 for ${sample}"

    python "${SQANTI3_DIR}/sqanti3_qc.py" \
        --gtf "$gtf" \
        --cpus 8 \
        "$REF_GTF" \
        "$REF_FASTA" \
        --short_reads "${short_r1},${short_r2}" \
        -o "$sample" \
        -d "$sample_out" \
        -n 8 \
        --isoAnnotLite \
        --report both

done

echo "All SQANTI3 jobs completed."
