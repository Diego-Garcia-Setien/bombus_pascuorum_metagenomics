#!/bin/bash
#SBATCH --job-name=11.2_coverm
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=16000

set -euo pipefail

######################################
# Script: 11.2_coverm.sh
#
# Quantifies the abundance of every dereplicated MAG (10.2_checkm2.sh)
# across all 93 samples, using the BAM files produced by
# 11_bowtie2_mapping.sh (reads mapped against the combined MAG
# catalog) and CoverM's "genome" mode.
#
# This is a SINGLE job, not an array: like GTDB-Tk/CheckM2/dRep, it
# needs to see every sample's BAM at once to build one abundance
# table (MAG x sample); it's not a per-sample calculation.
#
# NOTE: --genome-fasta-directory below points at
# mag_catalog/renamed_genomes/, NOT dRep's original
# dereplicated_genomes/. 11_bowtie2_mapping.sh renames contigs to be
# globally unique before mapping (original megahit contig IDs collide
# across samples), so CoverM needs the same renamed contig IDs to
# match what's in the BAM headers.
#
# Input:
#   data/09.CoverMResults/bam/<sample>.sorted.bam
#   data/09.CoverMResults/mag_catalog/renamed_genomes/*.fa
#
# Output:
#   data/09.CoverMResults/coverm_relative_abundance.tsv
######################################

WORKDIR=$(pwd)

OUTPUT_DIR="$WORKDIR/data/09.CoverMResults"
RENAMED_DIR="$OUTPUT_DIR/mag_catalog/renamed_genomes"
BAM_DIR="$OUTPUT_DIR/bam"

mkdir -p logs

CPU=$SLURM_CPUS_PER_TASK

N_BAMS=$(find "$BAM_DIR" -name "*.sorted.bam" | wc -l)
if [[ "$N_BAMS" -eq 0 ]]; then
    echo "ERROR: no BAM files found in $BAM_DIR (run 11_bowtie2_mapping.sh first)"
    exit 1
fi
echo "Found $N_BAMS mapped samples in $BAM_DIR"

if [[ ! -d "$RENAMED_DIR" ]]; then
    echo "ERROR: $RENAMED_DIR not found (run 11_bowtie2_mapping.sh first)"
    exit 1
fi

source /home/lchueca/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/lchueca/conda-env/coverm

coverm genome \
    --bam-files "$BAM_DIR"/*.sorted.bam \
    --genome-fasta-directory "$RENAMED_DIR" \
    --genome-fasta-extension fa \
    --threads "$CPU" \
    --output-file "$OUTPUT_DIR/coverm_relative_abundance.tsv"

conda deactivate

echo
echo "CoverM finished. Abundance table in $OUTPUT_DIR/coverm_relative_abundance.tsv"
echo
