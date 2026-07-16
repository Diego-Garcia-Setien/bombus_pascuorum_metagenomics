#!/bin/bash

#SBATCH --job-name=02_fastp
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out

#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=12000
#SBATCH --array=1-93%93

###############################################################################
# Script: 02_fastp.sh
#
# Description:
#   Adapter trimming and quality filtering with fastp.
#
# Input:
#   data/01.RawReads/
#
# Output:
#   data/02.CleanReads/
#
# Reports:
#   data/QC/02.Fastp/
#
###############################################################################

set -euo pipefail

#######################################
# Load software
#######################################

module load Miniforge3/24.11.3-2

conda activate /scratch/lchueca/conda-env/fastp

#######################################
# Settings
#######################################

CPU=$SLURM_CPUS_PER_TASK

#######################################
# Directories
#######################################

WORKDIR=$(pwd)

INPUT_DIR="$WORKDIR/data/01.RawReads"
OUTPUT_DIR="$WORKDIR/data/02.CleanReads"
REPORT_DIR="$WORKDIR/data/QC/02.Fastp"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"

#######################################
# Select sample
#######################################

cd "$INPUT_DIR"

SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=${SAMPLE#./}

if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: Sample not found."
    exit 1
fi

SAMPLE_DIR="$INPUT_DIR/$SAMPLE"

mkdir -p "$OUTPUT_DIR/$SAMPLE"

echo
echo "========================================="
echo "Processing: $SAMPLE"
echo "========================================="
echo

#######################################
# Locate FASTQ files
#######################################

cd "$SAMPLE_DIR"

R1=$(find . -maxdepth 1 -name "*_1.fq.gz" | head -1 || true)
R2=$(find . -maxdepth 1 -name "*_2.fq.gz" | head -1 || true)

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: FASTQ files not found."
    exit 1
fi

#######################################
# Output files
#######################################

OUT_R1="$OUTPUT_DIR/$SAMPLE/$(basename "$R1")"
OUT_R2="$OUTPUT_DIR/$SAMPLE/$(basename "$R2")"

FAILED="$OUTPUT_DIR/$SAMPLE/${SAMPLE}.failed.fq.gz"

HTML="$REPORT_DIR/${SAMPLE}.html"
JSON="$REPORT_DIR/${SAMPLE}.json"

#######################################
# Run fastp
#######################################

fastp \
    --thread "$CPU" \
    --in1 "$R1" \
    --in2 "$R2" \
    --out1 "$OUT_R1" \
    --out2 "$OUT_R2" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --trim_poly_x \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --length_required 50 \
    --correction \
    --failed_out "$FAILED" \
    --html "$HTML" \
    --json "$JSON" \
    --report_title "$SAMPLE"

echo
echo "Finished $SAMPLE"
echo
