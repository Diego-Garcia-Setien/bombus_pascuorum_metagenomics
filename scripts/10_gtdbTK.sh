#!/bin/bash
#SBATCH --job-name=10_gtdbTK
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=preemption
#SBATCH --qos=regular
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=160000

set -euo pipefail

######################################
# Script: 10_gtdbTK.sh
#
# Taxonomic classification of the refined MAGs from 09_MAGScoT.sh
# (MAGScoT score >= 0.5) using GTDB-Tk against the locally downloaded
# release232 reference package.
#
# This is a SINGLE job, not a per-sample array: GTDB-Tk's classify_wf
# has a large fixed cost (MSA loading, tree placement with pplacer)
# that is nearly independent of the number of genomes, and GTDB-Tk
# r232 needs >=140 GB RAM regardless of how many genomes are
# classified. Running it once for all 93 samples avoids paying that
# fixed cost 93 times. Refined bins from every sample are collected
# into one directory first (filenames prefixed with the sample name,
# so the per-genome rows in GTDB-Tk's output tables can still be
# traced back to their sample).
#
# NOTE: "preemption" is the only partition this account (bc3) has
# access to with enough RAM for GTDB-Tk r232 (general/quantum top out
# around 94 GB). Jobs on this partition can be requeued if preempted;
# GTDB-Tk classify_wf cannot resume mid-run, so a preempted job starts
# over from the beginning.
#
# Input:
#   data/06.MegahitResults/<sample>/final.contigs.fa
#   data/07.BinningResults/<sample>/magscot/<sample>.MAGScoT.refined.contig_to_bin.out
#
# Output:
#   data/08.GTDBtkResults/combined_refined_bins/<sample>__<bin>.fa
#   data/08.GTDBtkResults/gtdbtk_output/  (gtdbtk.bac120/ar53.summary.tsv, etc.)
######################################

WORKDIR=$(pwd)

CONTIGS_DIR="$WORKDIR/data/06.MegahitResults"
BINNING_DIR="$WORKDIR/data/07.BinningResults"
OUTPUT_DIR="$WORKDIR/data/08.GTDBtkResults"
COMBINED_BINS_DIR="$OUTPUT_DIR/combined_refined_bins"
GTDBTK_OUT_DIR="$OUTPUT_DIR/gtdbtk_output"

mkdir -p logs
rm -rf "$COMBINED_BINS_DIR"
mkdir -p "$COMBINED_BINS_DIR"

CPU=$SLURM_CPUS_PER_TASK

###############################################################################
# 1) Collect refined bins from every sample into one directory.
#    Filenames are "<sample>__<bin>.fa" so bins never collide across
#    samples and the sample of origin stays identifiable downstream.
###############################################################################

cd "$CONTIGS_DIR"

for SAMPLE_PATH in $(find . -mindepth 1 -maxdepth 1 -type d | sort); do
    SAMPLE=${SAMPLE_PATH#./}

    CONTIGS="$CONTIGS_DIR/$SAMPLE/final.contigs.fa"
    MAP="$BINNING_DIR/$SAMPLE/magscot/${SAMPLE}.MAGScoT.refined.contig_to_bin.out"

    if [[ ! -f "$MAP" ]]; then
        echo "WARNING: no refined MAGScoT output for $SAMPLE, skipping (run 09_MAGScoT.sh first)"
        continue
    fi

    awk -v mapfile="$MAP" -v outdir="$COMBINED_BINS_DIR" -v sample="$SAMPLE" '
        BEGIN {
            while ((getline line < mapfile) > 0) {
                if (line ~ /^binnew\t/) continue
                split(line, a, "\t")
                n = split(a[1], parts, "/")
                binname = parts[n]
                contigbin[a[2]] = binname
            }
        }
        /^>/ {
            cid = substr($1, 2)
            bin = contigbin[cid]
            printing = (bin != "")
            if (printing) outfile = outdir "/" sample "__" bin ".fa"
            if (printing) print ">" cid >> outfile
            next
        }
        {
            if (printing) print $0 >> outfile
        }
    ' "$CONTIGS"
done

cd "$WORKDIR"

N_BINS=$(find "$COMBINED_BINS_DIR" -name "*.fa" | wc -l)
echo "Collected $N_BINS refined bins from all samples into $COMBINED_BINS_DIR"

if [[ "$N_BINS" -eq 0 ]]; then
    echo "ERROR: no refined bins found, aborting."
    exit 1
fi

###############################################################################
# 2) Classify all collected bins in a single GTDB-Tk run
###############################################################################

source /home/lchueca/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/lchueca/conda-env/gtdbtk

export GTDBTK_DATA_PATH="/data/lchueca/databases/gtdbtk/release232"

mkdir -p "$GTDBTK_OUT_DIR/pplacer_scratch"

gtdbtk classify_wf \
    --genome_dir "$COMBINED_BINS_DIR" \
    --out_dir "$GTDBTK_OUT_DIR" \
    --extension fa \
    --cpus "$CPU" \
    --pplacer_cpus "$CPU" \
    --scratch_dir "$GTDBTK_OUT_DIR/pplacer_scratch"

conda deactivate

echo
echo "GTDB-Tk classification finished. Results in $GTDBTK_OUT_DIR"
echo
