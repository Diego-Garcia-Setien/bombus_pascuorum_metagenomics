#!/bin/bash

#SBATCH --job-name=07_kraken2
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out

#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=24000
#SBATCH --array=1-93%40

###############################################################################
#
# Kraken2 + Bracken
#
# Input:
#   data/03.MicrobiomeReads/<sample>/
#
# Output:
#
#   data/05.MicrobiotaTaxonomy/
#       *.kraken
#       *.report
#       summary/*.summary.tsv
#       summary/Kraken2_summary.tsv
#
#   data/05.BrackenTaxonomy/
#       *.bracken
#       *.bracken.report
#
# NOTE on memory:
#   The Kraken2 DB (/data/lchueca/databases/kraken_std) is ~93 GB, almost as
#   large as a full compute node's RAM (~94 GB). --memory-mapping makes
#   Kraken2 use mmap/page cache instead of loading the DB into the process's
#   own RAM, so the job can run within a modest --mem request.
#
###############################################################################

set -euo pipefail

module load Miniforge3/24.11.3-2

conda activate /scratch/lchueca/conda-env/kraken2

CPU=16

WORKDIR=$(pwd)

INPUT_DIR="$WORKDIR/data/03.MicrobiomeReads"

DB="/data/lchueca/databases/kraken_std"

KRAKEN_DIR="$WORKDIR/data/05.MicrobiotaTaxonomy"

BRACKEN_DIR="$WORKDIR/data/05.BrackenTaxonomy"

SUMMARY_DIR="$WORKDIR/data/05.MicrobiotaTaxonomy/summary"

mkdir -p "$KRAKEN_DIR"
mkdir -p "$BRACKEN_DIR"
mkdir -p "$SUMMARY_DIR"

###############################################################################
# Detect sample automatically (one subdirectory per sample, same layout
# as 01_quality_check.sh / 02_fastp.sh / 05_host_depletion.sh)
###############################################################################

TOTAL_TASKS=$SLURM_ARRAY_TASK_COUNT

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

echo
echo "=========================================="
echo "Sample: $SAMPLE"
echo "=========================================="
echo

###############################################################################
# Kraken2
###############################################################################

kraken2 \
    --db "$DB" \
    --paired \
    --gzip-compressed \
    --memory-mapping \
    --threads "$CPU" \
    --use-names \
    --report "$KRAKEN_DIR/${SAMPLE}.report" \
    --output "$KRAKEN_DIR/${SAMPLE}.kraken" \
    "$R1" \
    "$R2"

###############################################################################
# Bracken
###############################################################################

bracken \
    -d "$DB" \
    -i "$KRAKEN_DIR/${SAMPLE}.report" \
    -o "$BRACKEN_DIR/${SAMPLE}.bracken" \
    -w "$BRACKEN_DIR/${SAMPLE}.bracken.report" \
    -r 150 \
    -l S

###############################################################################
# Summary statistics
###############################################################################

TOTAL=$(awk 'END{print NR}' "$KRAKEN_DIR/${SAMPLE}.kraken")

CLASSIFIED=$(awk '$1=="C"{n++} END{print n+0}' "$KRAKEN_DIR/${SAMPLE}.kraken")

UNCLASSIFIED=$(awk '$1=="U"{n++} END{print n+0}' "$KRAKEN_DIR/${SAMPLE}.kraken")

CLASSIFIED_PERCENT=$(awk "BEGIN {printf \"%.2f\",100*$CLASSIFIED/$TOTAL}")

UNCLASSIFIED_PERCENT=$(awk "BEGIN {printf \"%.2f\",100*$UNCLASSIFIED/$TOTAL}")

echo -e "Sample\tTotal_reads\tClassified\tClassified_percent\tUnclassified\tUnclassified_percent" \
> "$SUMMARY_DIR/${SAMPLE}.summary.tsv"

echo -e "${SAMPLE}\t${TOTAL}\t${CLASSIFIED}\t${CLASSIFIED_PERCENT}\t${UNCLASSIFIED}\t${UNCLASSIFIED_PERCENT}" \
>> "$SUMMARY_DIR/${SAMPLE}.summary.tsv"

echo
echo "Finished $SAMPLE"
echo

###############################################################################
# Create global summary once every array task has finished.
#
# Completion is tracked via per-task marker files combined with a lock
# file (flock), so the global summary is guaranteed to be built exactly
# once, only after all tasks have completed, regardless of the order in
# which array tasks finish (same pattern as 01_quality_check.sh).
###############################################################################

DONE_DIR="$SUMMARY_DIR/.done"
mkdir -p "$DONE_DIR"

touch "$DONE_DIR/${SLURM_ARRAY_TASK_ID}.done"

(
    flock -n 200 || exit 0

    N_DONE=$(find "$DONE_DIR" -name "*.done" | wc -l)
    if [ "$N_DONE" -eq "$TOTAL_TASKS" ] && [ ! -f "$SUMMARY_DIR/.global_summary_done" ]; then
        echo
        echo "=========================================="
        echo "All $TOTAL_TASKS Kraken2/Bracken tasks finished."
        echo "Building global summary"
        echo "=========================================="
        echo

        head -1 "$SUMMARY_DIR/"*.summary.tsv | head -1 \
            > "$SUMMARY_DIR/Kraken2_summary.tsv"

        tail -q -n +2 "$SUMMARY_DIR/"*.summary.tsv \
            >> "$SUMMARY_DIR/Kraken2_summary.tsv"

        touch "$SUMMARY_DIR/.global_summary_done"
        rm -rf "$DONE_DIR"
    fi
) 200>"$SUMMARY_DIR/.summary.lock"
