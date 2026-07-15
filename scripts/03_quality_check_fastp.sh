#!/bin/bash
#SBATCH --job-name=03_quality_check
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=8000
#SBATCH --array=1-93%93      # Adjust to the total number of samples

###############################################################################
# Script: 01_quality_check.sh
#
# Description:
#   Perform quality assessment of raw Illumina FASTQ files using FastQC.
#   Each array task processes one sample. Once every array task has
#   finished, a single MultiQC run is triggered automatically.
#
#   Completion is tracked via per-task marker files combined with a
#   lock file (flock), so MultiQC is guaranteed to run exactly once,
#   only after all FastQC tasks have completed, regardless of the
#   order in which array tasks finish.
#
# Usage:
#   sbatch 01_quality_check.sh
#
###############################################################################
set -euo pipefail

# Load modules
module load FastQC/0.12.1-Java-11
module load Miniforge3/24.11.3-2

# Activate MultiQC environment
conda activate /scratch/lchueca/conda-env/MultiQC

# Settings
CPU=$SLURM_CPUS_PER_TASK
TOTAL_TASKS=$SLURM_ARRAY_TASK_COUNT

# Directories
WORKDIR=$(pwd)
CLEANDATA_DIR="$WORKDIR/data/02.CleanReads"
FASTQC_DIR="$WORKDIR/data/QC/03.FastQC_MultiQC"
FASTP_REP="$WORKDIR/data/QC/Fastp/"
FASTP_DIR="$WORKDIR/data/QC/03.Fastp_MultiQC"
DONE_DIR="$FASTQC_DIR/.done"

# Create output directories
mkdir -p "$FASTQC_DIR" "$FASTP_DIR" "$DONE_DIR"

# Sanity check: warn if the array size does not match the number of samples
N_SAMPLES=$(find "$CLEANDATA_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
if [ "$TOTAL_TASKS" -ne "$N_SAMPLES" ]; then
    echo "WARNING: --array size ($TOTAL_TASKS) does not match the number of samples ($N_SAMPLES)."
    echo "Update the #SBATCH --array directive accordingly."
fi

# Get current sample
cd "$CLEANDATA_DIR"
SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | \
         sort | \
         sed -n ${SLURM_ARRAY_TASK_ID}p)
SAMPLE=${SAMPLE#./}

echo
echo "=========================================="
echo "Running FastQC"
echo "Sample: $SAMPLE"
echo "=========================================="
echo

cd "$CLEANDATA_DIR/$SAMPLE"
fastqc \
    --threads "$CPU" \
    *.fq.gz \
    --outdir "$FASTQC_DIR"

echo
echo "FastQC completed for $SAMPLE"
echo

# Mark this task as finished
touch "$DONE_DIR/${SLURM_ARRAY_TASK_ID}.done"

# Check whether all array tasks have finished. The lock ensures only one
# task can evaluate/trigger MultiQC at a time, even if several tasks
# finish almost simultaneously.
(
    flock -n 200 || exit 0

    N_DONE=$(find "$DONE_DIR" -name "*.done" | wc -l)
    if [ "$N_DONE" -eq "$TOTAL_TASKS" ] && [ ! -f "$FASTQC_DIR/.multiqc_done" ]; then
        echo
        echo "=========================================="
        echo "All $TOTAL_TASKS FastQC tasks finished."
        echo "Running MultiQC"
        echo "=========================================="
        echo

        multiqc \
            "$FASTQC_DIR" \
            --outdir "$FASTQC_DIR" \
            --force

        touch "$FASTQC_DIR/.multiqc_done"
        rm -rf "$DONE_DIR"
    fi
) 200>"$FASTQC_DIR/.multiqc.lock"

#######################################
# Running MultiQC for fastp stats
#######################################

echo
echo "Generating MultiQC report..."
echo

multiqc \
    "$FASTP_REP" \
    --outdir "$FASTP_DIR" \
    --force

echo
echo "MultiQC report created successfully."
echo

echo
echo "Done."
