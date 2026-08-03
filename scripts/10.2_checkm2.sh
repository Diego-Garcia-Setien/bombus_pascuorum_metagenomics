#!/bin/bash
#SBATCH --job-name=10.2_checkm2
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mem=32000

set -euo pipefail

######################################
# Script: 10.2_checkm2.sh
#
# Quality assessment (CheckM2) and dereplication (dRep) of the refined
# MAGs from 09_MAGScoT.sh.
#
# 1) CheckM2 (completeness/contamination) against the locally
#    downloaded database (DIAMOND db version 3, matches checkm2 1.1.0).
# 2) dRep dereplicate, using CheckM2's quality_report.tsv as the
#    genome quality input (--genomeInfo) instead of letting dRep run
#    its own CheckM, and skani (dRep v4 default) for ANI comparisons.
#
# Both steps are SINGLE jobs, not per-sample arrays: CheckM2 has no
# large fixed cost per run (DIAMOND + a pretrained ML model, no
# reference tree), and dRep needs to compare genomes across ALL
# samples against each other anyway to dereplicate properly, so both
# only make sense run once over every sample's bins combined. Bin
# filenames are already sample-prefixed (set by 09_MAGScoT.sh via
# MAGScoT's -o option), so results stay traceable back to their
# sample of origin.
#
# Input:
#   data/07.BinningResults/<sample>/magscot/refined_bins/*.fa
#
# Output:
#   data/08.CheckM2Results/combined_refined_bins/     (symlinks)
#   data/08.CheckM2Results/checkm2_output/quality_report.tsv
#   data/08.CheckM2Results/genomeInfo.csv
#   data/08.CheckM2Results/drep_output/dereplicated_genomes/
#   data/08.CheckM2Results/drep_output/data_tables/Widb.csv
######################################

WORKDIR=$(pwd)

BINNING_DIR="$WORKDIR/data/07.BinningResults"
OUTPUT_DIR="$WORKDIR/data/08.CheckM2Results"
COMBINED_BINS_DIR="$OUTPUT_DIR/combined_refined_bins"
CHECKM2_OUT_DIR="$OUTPUT_DIR/checkm2_output"
CHECKM2_DB="/data/lchueca/databases/checkm2/202500220/CheckM2_database/uniref100.KO.1.dmnd"
GENOME_INFO="$OUTPUT_DIR/genomeInfo.csv"
DREP_OUT_DIR="$OUTPUT_DIR/drep_output"

mkdir -p logs
rm -rf "$COMBINED_BINS_DIR"
mkdir -p "$COMBINED_BINS_DIR"

CPU=$SLURM_CPUS_PER_TASK

###############################################################################
# 1) Collect refined bins from every sample into one directory (symlinks,
#    filenames are already unique and sample-prefixed from 09_MAGScoT.sh)
###############################################################################

find "$BINNING_DIR" -path "*/magscot/refined_bins/*.fa" -exec ln -s {} "$COMBINED_BINS_DIR/" \;

N_BINS=$(find "$COMBINED_BINS_DIR" -name "*.fa" | wc -l)
echo "Collected $N_BINS refined bins from all samples into $COMBINED_BINS_DIR"

if [[ "$N_BINS" -eq 0 ]]; then
    echo "ERROR: no refined bins found, aborting (run 09_MAGScoT.sh first)."
    exit 1
fi

###############################################################################
# 2) Run CheckM2 on all collected bins in a single call
###############################################################################

source /home/lchueca/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/lchueca/conda-env/checkm2

checkm2 predict \
    --threads "$CPU" \
    --input "$COMBINED_BINS_DIR" \
    --extension .fa \
    --output-directory "$CHECKM2_OUT_DIR" \
    --database_path "$CHECKM2_DB" \
    --force

conda deactivate

echo
echo "CheckM2 finished. Quality report in $CHECKM2_OUT_DIR/quality_report.tsv"
echo

###############################################################################
# 3) Build dRep's genomeInfo table from CheckM2's quality report
#    (genome,completeness,contamination -- genome must include the
#    .fa extension to match the filenames in COMBINED_BINS_DIR)
###############################################################################

awk -F'\t' 'NR==1 {print "genome,completeness,contamination"; next} {print $1".fa,"$2","$3}' \
    "$CHECKM2_OUT_DIR/quality_report.tsv" > "$GENOME_INFO"

###############################################################################
# 4) Dereplicate all collected bins with dRep (skani-based, v4 defaults)
###############################################################################

conda activate /scratch/lchueca/conda-env/drep

dRep dereplicate \
    "$DREP_OUT_DIR" \
    -g "$COMBINED_BINS_DIR"/*.fa \
    --genomeInfo "$GENOME_INFO" \
    -p "$CPU"

conda deactivate

echo
echo "dRep finished. Dereplicated genomes in $DREP_OUT_DIR/dereplicated_genomes/"
echo
