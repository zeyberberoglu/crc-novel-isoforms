#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=salmon_quant_%j.out
#SBATCH --error=salmon_quant_%j.err
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G

set -euo pipefail

module load miniforge
source activate longread_de

WORK="/path/to/project/diff_exp_salmon"
GFF="/path/to/project/final_extended_annotation.gff"
GENOME="/path/to/reference/GRCh38.p14.genome.fa"
FASTQ_ROOT="/path/to/chitebuster_runs"

INDEX="${WORK}/salmon_index"
BUILD="${WORK}/salmon_build"
THREADS=24

mkdir -p "${WORK}/quant" "$BUILD"
cd "$WORK"

if [[ ! -s "${INDEX}/info.json" ]]; then

    gffread \
        -w "${BUILD}/transcripts.fa" \
        -g "$GENOME" \
        "$GFF"

    grep "^>" "$GENOME" \
        | cut -d ' ' -f1 \
        | sed 's/^>//' \
        > "${BUILD}/decoys.txt"

    cat "${BUILD}/transcripts.fa" "$GENOME" \
        | pigz -p "$THREADS" \
        > "${BUILD}/gentrome.fa.gz"

    salmon index \
        -t "${BUILD}/gentrome.fa.gz" \
        -d "${BUILD}/decoys.txt" \
        -i "$INDEX" \
        -p "$THREADS"

fi

for D in "${FASTQ_ROOT}"/chitebuster_output_*; do

    S=$(basename "$D" | sed 's/chitebuster_output_//')

    R1="${D}/${S}/1_fastp/${S}_reads1.fastq.gz"
    R2="${D}/${S}/1_fastp/${S}_reads2.fastq.gz"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "Skipping ${S} because FASTQ files are missing"
        continue
    fi

    if [[ -s "${WORK}/quant/${S}/quant.sf" ]]; then
        echo "Skipping ${S} because quant.sf already exists"
        continue
    fi

    salmon quant \
        -i "$INDEX" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS" \
        --seqBias \
        --gcBias \
        --posBias \
        -o "${WORK}/quant/${S}"

done

Rscript "${WORK}/salmon_deseq2.R"
