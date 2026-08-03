#!/bin/bash
#SBATCH --job-name=11_bowtie2_mapping
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=16000
#SBATCH --array=1-93%25

set -euo pipefail

######################################
# Script: 11_bowtie2_mapping.sh
#
# Maps every sample's microbiome reads against the dereplicated MAG
# catalog produced by 10.2_checkm2.sh (dRep's dereplicated_genomes/),
# so their abundance across all 93 samples can be quantified
# afterwards with 11.2_coverm.sh.
#
# This is a per-sample SLURM array (same pattern as 07/08/09), since
# unlike GTDB-Tk/CheckM2/dRep, mapping has no large fixed cost that
# would be wasted by repeating it per sample - it genuinely benefits
# from running in parallel across the cluster.
#
# The MAG catalog + its Bowtie2 index are shared by every task, so
# they are only built once: array task 1 builds them, every other
# task waits for a "index.done" marker file before mapping (SLURM
# array tasks all start around the same time, so a naive "build if
# missing" in every task would race).
#
# IMPORTANT: contig IDs inside the dereplicated genomes come straight
# from each sample's own megahit assembly (e.g. "k127_34299"), and
# different samples independently reuse the same IDs. Concatenating
# the genomes as-is produces duplicate reference names, which
# bowtie2/samtools reject outright. So each genome is first copied
# into mag_catalog/renamed_genomes/ with contigs renamed to
# "<genome_filename>__<original_contig_id>" (globally unique), and
# THAT copy is what gets concatenated into the Bowtie2 reference.
# 11.2_coverm.sh must point --genome-fasta-directory at this same
# renamed_genomes/ directory, not at dRep's original output, so its
# contig names match what's in the BAM headers.
#
# Input:
#   data/03.MicrobiomeReads/<sample>/<sample>_microbiome_R{1,2}.fastq.gz
#   data/08.CheckM2Results/drep_output/dereplicated_genomes/*.fa
#
# Output:
#   data/09.CoverMResults/mag_catalog/renamed_genomes/*.fa (contig IDs made unique)
#   data/09.CoverMResults/mag_catalog/combined_MAGs.fasta + Bowtie2 index
#   data/09.CoverMResults/bam/<sample>.sorted.bam(.bai)
######################################

WORKDIR=$(pwd)

READS_DIR="$WORKDIR/data/03.MicrobiomeReads"
DEREP_DIR="$WORKDIR/data/08.CheckM2Results/drep_output/dereplicated_genomes"
OUTPUT_DIR="$WORKDIR/data/09.CoverMResults"
CATALOG_DIR="$OUTPUT_DIR/mag_catalog"
RENAMED_DIR="$CATALOG_DIR/renamed_genomes"
BAM_DIR="$OUTPUT_DIR/bam"
COMBINED_REF="$CATALOG_DIR/combined_MAGs.fasta"
INDEX="$CATALOG_DIR/mag_catalog"
INDEX_DONE="$CATALOG_DIR/index.done"

mkdir -p logs "$CATALOG_DIR" "$BAM_DIR"

CPU=$SLURM_CPUS_PER_TASK

module load Bowtie2/2.5.5-GCC-14.2.0
module load SAMtools/1.18-GCC-12.3.0

###############################################################################
# 0) Build the combined MAG reference + Bowtie2 index ONCE (array task 1),
#    every other task waits for the "index.done" marker.
###############################################################################

if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
    if [[ ! -f "$INDEX_DONE" ]]; then
        N_GENOMES=$(find "$DEREP_DIR" -name "*.fa" | wc -l)
        if [[ "$N_GENOMES" -eq 0 ]]; then
            echo "ERROR: no dereplicated genomes found in $DEREP_DIR (run 10.2_checkm2.sh first)"
            exit 1
        fi
        echo "Building combined MAG reference from $N_GENOMES dereplicated genomes"

        mkdir -p "$RENAMED_DIR"
        for f in "$DEREP_DIR"/*.fa; do
            BASENAME=$(basename "$f" .fa)
            awk -v prefix="$BASENAME" '/^>/{sub(/^>/, ">" prefix "__")} {print}' "$f" > "$RENAMED_DIR/$(basename "$f")"
        done

        cat "$RENAMED_DIR"/*.fa > "$COMBINED_REF"
        bowtie2-build --threads "$CPU" "$COMBINED_REF" "$INDEX"
        touch "$INDEX_DONE"
    fi
else
    echo "Waiting for the MAG Bowtie2 index to be built by array task 1..."
    WAIT_ELAPSED=0
    WAIT_TIMEOUT=3600
    while [[ ! -f "$INDEX_DONE" ]]; do
        if [[ "$WAIT_ELAPSED" -ge "$WAIT_TIMEOUT" ]]; then
            echo "ERROR: timed out waiting for the MAG Bowtie2 index"
            exit 1
        fi
        sleep 15
        WAIT_ELAPSED=$((WAIT_ELAPSED+15))
    done
fi

###############################################################################
# 1) Detect the sample (same scheme as 07/08/09)
###############################################################################

cd "$READS_DIR"

SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=${SAMPLE#./}

if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: Sample not found."
    exit 1
fi

cd "$WORKDIR"

R1=$(find "$READS_DIR/$SAMPLE" -maxdepth 1 -name "*_microbiome_R1.fastq.gz" | head -1 || true)
R2=$(find "$READS_DIR/$SAMPLE" -maxdepth 1 -name "*_microbiome_R2.fastq.gz" | head -1 || true)

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: FASTQ files not found for $SAMPLE"
    echo "$READS_DIR/$SAMPLE"
    exit 1
fi

###############################################################################
# 2) Map this sample's reads against the combined MAG catalog
###############################################################################

BAM="$BAM_DIR/${SAMPLE}.sorted.bam"

echo
echo "======================================"
echo "Mapping sample: $SAMPLE"
echo "======================================"
echo

bowtie2 \
    --threads "$CPU" \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    2> "$BAM_DIR/${SAMPLE}.bowtie2.log" \
    | samtools view -@ "$CPU" -b - \
    | samtools sort -@ "$CPU" -o "$BAM"

samtools index -@ "$CPU" "$BAM"

echo
echo "Mapping finished for $SAMPLE"
echo
