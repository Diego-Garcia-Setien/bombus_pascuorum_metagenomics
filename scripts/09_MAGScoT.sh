#!/bin/bash
#SBATCH --job-name=09_MAGScoT
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --mem=16000
#SBATCH --array=1-93%25

set -euo pipefail

######################################
# Script: 09_MAGScoT.sh
#
# Runs MAGScoT to score and refine the bins produced by the three
# binners in 08_binning.sh (MetaBAT2, MaxBin2, CONCOCT), for every
# sample in data/06.MegahitResults.
#
# MAGScoT needs two inputs per sample:
#   1) A "contigs_to_bin" table (bin <TAB> contig <TAB> binner, no
#      header), built here from the bin FASTA files of the three
#      binners.
#   2) A marker-gene hit table (gene_id <TAB> marker <TAB> e-value),
#      built here by running prodigal on the assembly contigs and
#      hmmsearch against the GTDB r207 marker HMMs bundled with
#      MAGScoT.
#
# Input:
#   data/06.MegahitResults/<sample>/final.contigs.fa
#   data/07.BinningResults/<sample>/metabat2/bin.*.fa
#   data/07.BinningResults/<sample>/maxbin2/bin.*.fasta
#   data/07.BinningResults/<sample>/concoct/fasta_bins/*.fa
#
# Output:
#   data/07.BinningResults/<sample>/magscot/
#       <sample>.contigs_to_bin.tsv
#       <sample>.hmm
#       <sample>.MAGScoT.*  (MAGScoT scoring/refinement results)
######################################

WORKDIR=$(pwd)

CONTIGS_DIR="$WORKDIR/data/06.MegahitResults"
BINNING_DIR="$WORKDIR/data/07.BinningResults"
MAGSCOT_REPO="/scratch/lchueca/conda-env/magscot/MAGScoT"

mkdir -p logs

CPU=$SLURM_CPUS_PER_TASK

###############################################################################
# Automatically detect the sample (one subfolder per sample in
# 06.MegahitResults, same scheme as 07_genome_assembly.sh / 08_binning.sh)
###############################################################################

cd "$CONTIGS_DIR"

SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=${SAMPLE#./}

if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: Sample not found."
    exit 1
fi

cd "$WORKDIR"

CONTIGS="$CONTIGS_DIR/$SAMPLE/final.contigs.fa"
SAMPLE_BINNING_DIR="$BINNING_DIR/$SAMPLE"
METABAT_DIR="$SAMPLE_BINNING_DIR/metabat2"
MAXBIN_DIR="$SAMPLE_BINNING_DIR/maxbin2"
CONCOCT_BINS_DIR="$SAMPLE_BINNING_DIR/concoct/fasta_bins"

if [[ ! -f "$CONTIGS" ]]; then
    echo "ERROR: Contigs not found for $SAMPLE"
    echo "$CONTIGS"
    exit 1
fi

if [[ ! -d "$METABAT_DIR" || ! -d "$MAXBIN_DIR" || ! -d "$CONCOCT_BINS_DIR" ]]; then
    echo "ERROR: Missing binning results for $SAMPLE (run 08_binning.sh first)"
    echo "$SAMPLE_BINNING_DIR"
    exit 1
fi

MAGSCOT_DIR="$SAMPLE_BINNING_DIR/magscot"
mkdir -p "$MAGSCOT_DIR"

echo
echo "======================================"
echo "MAGScoT for sample: $SAMPLE"
echo "======================================"
echo

source /home/lchueca/miniforge3/etc/profile.d/conda.sh
conda activate magscot

###############################################################################
# 1) Build the contigs_to_bin table from the three binners.
#    Bin names are prefixed with the binner name so IDs never collide
#    across binners even if two binners happen to number bins alike.
###############################################################################

CONTIGS_TO_BIN="$MAGSCOT_DIR/${SAMPLE}.contigs_to_bin.tsv"
> "$CONTIGS_TO_BIN"

for f in "$METABAT_DIR"/bin.*.fa; do
    [[ -e "$f" ]] || continue
    BIN_NAME="metabat2_$(basename "$f" .fa)"
    grep ">" "$f" | sed 's/^>//' | awk -v b="$BIN_NAME" -v s="metabat2" '{print b"\t"$1"\t"s}' >> "$CONTIGS_TO_BIN"
done

for f in "$MAXBIN_DIR"/bin.*.fasta; do
    [[ -e "$f" ]] || continue
    BIN_NAME="maxbin2_$(basename "$f" .fasta)"
    grep ">" "$f" | sed 's/^>//' | awk -v b="$BIN_NAME" -v s="maxbin2" '{print b"\t"$1"\t"s}' >> "$CONTIGS_TO_BIN"
done

for f in "$CONCOCT_BINS_DIR"/*.fa; do
    [[ -e "$f" ]] || continue
    BIN_NAME="concoct_$(basename "$f" .fa)"
    grep ">" "$f" | sed 's/^>//' | awk -v b="$BIN_NAME" -v s="concoct" '{print b"\t"$1"\t"s}' >> "$CONTIGS_TO_BIN"
done

###############################################################################
# 2) Marker gene identification (GTDBtk rel 207), same recipe as the
#    MAGScoT README: prodigal for ORF calling, hmmsearch against the
#    bundled TIGRFAM/Pfam marker HMMs.
###############################################################################

PRODIGAL_FAA="$MAGSCOT_DIR/${SAMPLE}.prodigal.faa"
PRODIGAL_FFN="$MAGSCOT_DIR/${SAMPLE}.prodigal.ffn"

prodigal -p meta -i "$CONTIGS" -a "$PRODIGAL_FAA" -d "$PRODIGAL_FFN" -o "$MAGSCOT_DIR/${SAMPLE}.prodigal.tmp"

TIGR_OUT="$MAGSCOT_DIR/${SAMPLE}.hmm.tigr.hit.out"
PFAM_OUT="$MAGSCOT_DIR/${SAMPLE}.hmm.pfam.hit.out"

hmmsearch -o "$MAGSCOT_DIR/${SAMPLE}.hmm.tigr.out" --tblout "$TIGR_OUT" --noali --notextw --cut_nc --cpu "$CPU" \
    "$MAGSCOT_REPO/hmm/gtdbtk_rel207_tigrfam.hmm" "$PRODIGAL_FAA"

hmmsearch -o "$MAGSCOT_DIR/${SAMPLE}.hmm.pfam.out" --tblout "$PFAM_OUT" --noali --notextw --cut_nc --cpu "$CPU" \
    "$MAGSCOT_REPO/hmm/gtdbtk_rel207_Pfam-A.hmm" "$PRODIGAL_FAA"

TIGR_HITS="$MAGSCOT_DIR/${SAMPLE}.tigr"
PFAM_HITS="$MAGSCOT_DIR/${SAMPLE}.pfam"
HMM_FILE="$MAGSCOT_DIR/${SAMPLE}.hmm"

grep -v "^#" "$TIGR_OUT" | awk '{print $1"\t"$3"\t"$5}' > "$TIGR_HITS"
grep -v "^#" "$PFAM_OUT" | awk '{print $1"\t"$4"\t"$5}' > "$PFAM_HITS"
cat "$PFAM_HITS" "$TIGR_HITS" > "$HMM_FILE"

###############################################################################
# 3) Run MAGScoT: scoring and refinement of the combined bin set
#
#    MAGScoT's bin-merging step has a known bug that crashes with
#    "Error in `left_join()`: ... Problem with `contig`" on some
#    samples (contig/bin combinations that trigger an empty
#    intermediate data frame during merging). If that happens, retry
#    once with --skip_merge_bins, which avoids the buggy code path
#    (bins are still scored and refined, just never merged across
#    binners).
###############################################################################

if ! Rscript "$MAGSCOT_REPO/MAGScoT.R" \
    -i "$CONTIGS_TO_BIN" \
    --hmm "$HMM_FILE" \
    -o "$MAGSCOT_DIR/${SAMPLE}.MAGScoT"; then
    echo "WARNING: MAGScoT failed for $SAMPLE with bin-merging enabled, retrying with --skip_merge_bins"
    Rscript "$MAGSCOT_REPO/MAGScoT.R" \
        -i "$CONTIGS_TO_BIN" \
        --hmm "$HMM_FILE" \
        -o "$MAGSCOT_DIR/${SAMPLE}.MAGScoT" \
        --skip_merge_bins
fi

conda deactivate

###############################################################################
# 4) Extract a FASTA file for each refined bin.
#    MAGScoT only outputs a contig-to-bin table (refined.contig_to_bin.out),
#    not the bin sequences themselves, so they are cut out of the original
#    assembly contigs here.
###############################################################################

REFINED_MAP="$MAGSCOT_DIR/${SAMPLE}.MAGScoT.refined.contig_to_bin.out"
REFINED_BINS_DIR="$MAGSCOT_DIR/refined_bins"

rm -rf "$REFINED_BINS_DIR"
mkdir -p "$REFINED_BINS_DIR"

awk -v mapfile="$REFINED_MAP" -v outdir="$REFINED_BINS_DIR" '
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
        if (printing) outfile = outdir "/" bin ".fa"
        if (printing) print ">" cid >> outfile
        next
    }
    {
        if (printing) print $0 >> outfile
    }
' "$CONTIGS"

echo
echo "MAGScoT finished for $SAMPLE"
echo
