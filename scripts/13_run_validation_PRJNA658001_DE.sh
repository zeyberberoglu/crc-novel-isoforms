#!/bin/bash
#SBATCH -o slurm-%j-CRC_DE.out
#SBATCH --time=240:00:00
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G

module load miniforge
conda activate longread_de

# ==============================
# User-defined paths
# ==============================

BASE_DIR="/path/to/PRJNA658001"
DE_SCRIPT="${BASE_DIR}/differential_expression/12_validation_PRJNA658001_deseq2.R"

cd "$BASE_DIR" || exit 1

Rscript "$DE_SCRIPT"
