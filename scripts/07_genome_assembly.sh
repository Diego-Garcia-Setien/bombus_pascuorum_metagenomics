#!/bin/bash
#SBATCH --job-name=07_genome_assembly
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:00:00
#SBATCH --mem=32000
#SBATCH --array=1-93%25

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

INPUT_DIR="$WORKDIR/data/03.MicrobiomeReads"
OUTPUT_DIR="$WORKDIR/data/06.MegahitResults"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

###############################################################################
# Detectar la muestra automáticamente (una subcarpeta por muestra, mismo
# esquema que 01_quality_check.sh / 02_fastp.sh / 05_host_depletion.sh)
###############################################################################

cd "$INPUT_DIR"

SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=${SAMPLE#./}

if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: Sample not found."
    exit 1
fi

SAMPLE_DIR="$INPUT_DIR/$SAMPLE"

R1=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*_microbiome_R1.fastq.gz" | head -1 || true)
R2=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*_microbiome_R2.fastq.gz" | head -1 || true)

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: FASTQ files not found."
    echo "$SAMPLE_DIR"
    exit 1
fi

OUT_DIR="$OUTPUT_DIR/${SAMPLE}"

megahit -1 "$R1" -2 "$R2" -o "$OUT_DIR" -t "$CPU" --presets meta-large
