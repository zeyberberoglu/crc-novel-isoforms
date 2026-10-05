#!/bin/bash
#SBATCH -o slurm-%j-star_eastr.out
#SBATCH --job-name=star_eastr
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --mem=90G

# ==============================
# User-defined paths
# ==============================

BASE_DIR="/path/to/project"
REFERENCE_DIR="/path/to/reference"

FASTQ_DIR="${BASE_DIR}/fastp_trimmed"

STAR_INDEX="${REFERENCE_DIR}/star_hg38"
REFERENCE_FA="${REFERENCE_DIR}/GRCh38.p14.genome.fa"
ANNOTATION="${REFERENCE_DIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3"

STAR_DIR="${BASE_DIR}/star_alignment"
EASTR_DIR="${BASE_DIR}/eastr_output"
LOG_DIR="${BASE_DIR}/logs_star_eastr"

# ==============================
# Conda environment
# ==============================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eastr_env

# ==============================
# Sample settings
# ==============================

SAMPLE="SRR13518175"
CORES=8

R1="${FASTQ_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"

SAMPLE_STAR_DIR="${STAR_DIR}/${SAMPLE}"
BAM_FILE="${SAMPLE_STAR_DIR}/Aligned.sortedByCoord.out.bam"
LOG_FILE="${LOG_DIR}/${SAMPLE}.log"

mkdir -p "$SAMPLE_STAR_DIR" "$EASTR_DIR" "$LOG_DIR"

cd "$BASE_DIR" || exit 1

{
    echo "--- Starting STAR for ${SAMPLE} at $(date) ---"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "ERROR: Missing R1 or R2 for ${SAMPLE}"
        exit 1
    fi

    STAR \
        --runThreadN "$CORES" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --limitBAMsortRAM 103079215104 \
        --twopassMode Basic \
        --sjdbGTFfile "$ANNOTATION" \
        --outFilterType BySJout \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${SAMPLE_STAR_DIR}/" \
        --winAnchorMultimapNmax 200 \
        --outFilterMultimapNmax 100 \
        --sjdbOverhang 150

    if [ $? -ne 0 ]; then
        echo "ERROR: STAR failed for ${SAMPLE}"
        exit 1
    fi

    samtools index "$BAM_FILE"

    if [ $? -ne 0 ]; then
        echo "ERROR: samtools index failed for ${SAMPLE}"
        exit 1
    fi

    echo "--- Running EASTR for ${SAMPLE} at $(date) ---"

    eastr \
        --bam "$BAM_FILE" \
        --reference "$REFERENCE_FA" \
        --output "${EASTR_DIR}/${SAMPLE}"

    if [ $? -ne 0 ]; then
        echo "ERROR: EASTR failed for ${SAMPLE}"
        exit 1
    fi

    echo "--- Finished STAR + EASTR for ${SAMPLE} at $(date) ---"

} >> "$LOG_FILE" 2>&1
