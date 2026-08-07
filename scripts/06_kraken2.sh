#!/bin/bash
#SBATCH --job-name=06_kraken2
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=120000
#SBATCH --array=1-93%25

set -euo pipefail

#################################
# Cargar software
#################################
module load Miniforge3/24.11.3-2
conda activate /scratch/lchueca/conda-env/kraken2

CPU=$SLURM_CPUS_PER_TASK

######################################
# Script: kraken2_taxonomy.sh
#
# Obtener la taxonomía de las secuencias no alineadas
# con el genoma del hospedador obtenidas mediante bowtie2
######################################

WORKDIR=$(pwd)
INPUT_DIR="$WORKDIR/data/03.MicrobiomeReads"
KRAKEN_OUT="$WORKDIR/data/05.MicrobiotaTaxonomy"
BRACKEN_OUT="$WORKDIR/data/05.BrackenTaxonomy"
DATABASE="/data/lchueca/databases/kraken_std"	# Here place kraken2 database

mkdir -p "$KRAKEN_OUT"
mkdir -p "$BRACKEN_OUT"
mkdir -p logs

###############################################################################
# Detect sample automatically (one subdirectory per sample, same layout
# as 01_quality_check.sh / 02_fastp.sh / 05_host_depletion.sh)
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

echo
echo "=========================================="
echo "Sample: $SAMPLE"
echo "=========================================="
echo

# Procesando la muestra con kraken2

# Las secuencias clasificadas o no clasificadas se pueden 
# enviar a un archivo para su posterior procesamiento, utilizando los interruptores --classified-out y, 
#respectivamente.--unclassified-out


kraken2 --db "$DATABASE" \
      --threads "$CPU" --paired --minimum-hit-groups 2 \
      --output "$KRAKEN_OUT/${SAMPLE}.kraken2.out" \
      --report "$KRAKEN_OUT/${SAMPLE}.kraken2.report" \
      --gzip-compressed "$R1" "$R2"


echo "Clasificación taxonómica de $SAMPLE terminada"

# Vamos a usar Bracken, que es un programa complementario de Kraken2, 
# Sirve para estimar la abundancia en un solo nivel taxonómico

bracken -d "$DATABASE" -i "$KRAKEN_OUT/${SAMPLE}.kraken2.report" \
      -o "$BRACKEN_OUT/${SAMPLE}.bracken_output" -w "$BRACKEN_OUT/${SAMPLE}.bracken.kreport" -l S \
      -t "$CPU"

###############################################################################
# Combine every sample's kreport into the family/genus/species table and
# the per-sample bacterial fraction summary, once they're all done.
#
# This is a SLURM array (one task per sample), so there's no single task
# that "runs last" on purpose - instead, every task checks after finishing
# its own sample whether all of them are now present, and (re)writes the
# combined files if so. Safe to do redundantly: the logic is deterministic
# (same as the standalone combine_bracken.sh), so it doesn't matter if
# several tasks race to regenerate it at nearly the same time.
#
# "cellular organisms" (taxid 131567) and "Bacteria" (taxid 2) are matched
# by their fixed NCBI taxonomy IDs rather than rank code: in this Kraken2
# database Bacteria sits under a no-rank R1/R2 lineage (cellular organisms
# -> Bacteria), not the "D" (domain) rank code you might expect.
###############################################################################

TOTAL_SAMPLES=$(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
N_KREPORTS=$(find "$BRACKEN_OUT" -name "*.bracken.kreport" | wc -l)

if [[ "$N_KREPORTS" -eq "$TOTAL_SAMPLES" ]]; then
    TAXONOMY_OUT="$BRACKEN_OUT/combined_bracken_family_genus_species.tsv"
    BACTERIA_OUT="$BRACKEN_OUT/bacteria_fraction_per_sample.tsv"

    echo -e "Sample\tRank\tTaxID\tName\tPercentage\tReads_clade" > "$TAXONOMY_OUT"
    echo -e "Sample\tCellular_organisms_pct\tBacteria_pct\tNonBacteria_pct" > "$BACTERIA_OUT"

    for KREPORT in "$BRACKEN_OUT"/*.bracken.kreport; do
        S=$(basename "$KREPORT" .bracken.kreport)

        awk -F'\t' -v sample="$S" '$4=="F" || $4=="G" || $4=="S" {
            name=$6
            sub(/^[ \t]+/, "", name)
            print sample"\t"$4"\t"$5"\t"name"\t"$1"\t"$2
        }' "$KREPORT" >> "$TAXONOMY_OUT"

        awk -F'\t' -v sample="$S" '
            $5=="131567" { cellular=$1 }
            $5=="2"      { bacteria=$1 }
            END {
                nonbacteria = cellular - bacteria
                printf "%s\t%.2f\t%.2f\t%.2f\n", sample, cellular, bacteria, nonbacteria
            }
        ' "$KREPORT" >> "$BACTERIA_OUT"
    done

    echo "All $TOTAL_SAMPLES sample kreports present - wrote $TAXONOMY_OUT and $BACTERIA_OUT"
fi

