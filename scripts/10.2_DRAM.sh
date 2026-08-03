#!/bin/bash
#SBATCH --job-name=10.2_DRAM
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=64000

set -euo pipefail

######################################
# Script: 10.2_DRAM.sh
#
# Functional annotation (KOfam/Pfam/dbCAN/VOGDB/MEROPS) of the
# dereplicated MAG catalog from 10.1_checkm2.sh, using DRAM.
#
# REQUIRES DRAM_setup.sh to have completed successfully first
# (check with `DRAM-setup.py print_config` - every database location
# should be populated, not None).
#
# This is a SINGLE job, not an array: with ~20 dereplicated genomes,
# annotating them together in one `DRAM.py annotate` call (it accepts
# a wildcard covering multiple FASTA files) avoids reloading DRAM's
# several-GB reference databases once per genome.
#
# Input:
#   data/08.CheckM2Results/drep_output/dereplicated_genomes/*.fa
#
# Output:
#   data/10.DRAMResults/annotation/annotations.tsv (+ genes.faa, genbank/, trnas.tsv, rrnas.tsv)
#   data/10.DRAMResults/distillation/genome_stats.tsv, metabolism_summary.xlsx, product.html
######################################

WORKDIR=$(pwd)

DEREP_DIR="$WORKDIR/data/08.CheckM2Results/drep_output/dereplicated_genomes"
OUTPUT_DIR="$WORKDIR/data/10.DRAMResults"
ANNOTATE_DIR="$OUTPUT_DIR/annotation"
DISTILL_DIR="$OUTPUT_DIR/distillation"

mkdir -p logs "$OUTPUT_DIR"

CPU=$SLURM_CPUS_PER_TASK

N_GENOMES=$(find "$DEREP_DIR" -name "*.fa" | wc -l)
if [[ "$N_GENOMES" -eq 0 ]]; then
    echo "ERROR: no dereplicated genomes found in $DEREP_DIR (run 10.1_checkm2.sh first)"
    exit 1
fi
echo "Annotating $N_GENOMES dereplicated genomes with DRAM"

if [[ -d "$ANNOTATE_DIR" ]]; then
    echo "ERROR: $ANNOTATE_DIR already exists. DRAM.py annotate refuses to write into an"
    echo "existing output directory. Remove it (or rename it) before re-running."
    exit 1
fi

source /home/lchueca/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/lchueca/conda-env/DRAM

###############################################################################
# 1) Annotate: gene calling (prodigal) + search against KOfam/Pfam/dbCAN/
#    VOGDB/MEROPS/RefSeq-viral
###############################################################################

DRAM.py annotate \
    -i "$DEREP_DIR/*.fa" \
    -o "$ANNOTATE_DIR" \
    --threads "$CPU"

###############################################################################
# 2) Distill: summarize the annotation table into genome-level metabolism
#    stats and the interactive product.html
###############################################################################

DISTILL_ARGS=(-i "$ANNOTATE_DIR/annotations.tsv" -o "$DISTILL_DIR")
[[ -f "$ANNOTATE_DIR/rrnas.tsv" ]] && DISTILL_ARGS+=(--rrna_path "$ANNOTATE_DIR/rrnas.tsv")
[[ -f "$ANNOTATE_DIR/trnas.tsv" ]] && DISTILL_ARGS+=(--trna_path "$ANNOTATE_DIR/trnas.tsv")

DRAM.py distill "${DISTILL_ARGS[@]}"

conda deactivate

echo
echo "DRAM finished."
echo "Annotations: $ANNOTATE_DIR/annotations.tsv"
echo "Distillation: $DISTILL_DIR/"
echo
