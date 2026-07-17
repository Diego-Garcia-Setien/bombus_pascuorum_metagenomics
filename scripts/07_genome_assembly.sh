#!/bin/bash
#SBATCH --job-name=bp_megahit_tfm
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=12000
#SBATCH --array=1-93%93

set -euo pipefail

##################################
# Cargar software
##################################
module load Miniforge3/24.11.3-2
conda activate /scratch/lchueca/conda-env/megahit

CPU=$SLURM_CPUS_PER_TASK

######################################
# Script: script_megahit.sh
#
# Ensamblar las lecturas cortas en contigs
######################################

WORKDIR=$(pwd)

INPUT_DIR="$WORKDIR/data/microbiome_reads"
OUTPUT_DIR="$WORKDIR/data/megahit_results"
SAMPLES_LIST="$WORKDIR/data/samples.txt"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

###############################################################################
# Detectar la muestra automáticamente
###############################################################################
R1=$(find "$INPUT_DIR" -maxdepth 1 -name "*_microbiota_1.fastq.gz" \
      | sort \
      | sed -n "${SLURM_ARRAY_TASK_ID}p")

SAMPLE=$(basename "$R1")
SAMPLE=${SAMPLE%_microbiota_1.fastq.gz}

R1="$INPUT_DIR/${SAMPLE}_microbiota_1.fastq.gz"
R2="$INPUT_DIR/${SAMPLE}_microbiota_2.fastq.gz"
OUT_DIR="$OUTPUT_DIR/${SAMPLE}"

megahit -1 "$R1" -2 "$R2" -o "$OUT_DIR" -t "$CPU" --presets meta-large
