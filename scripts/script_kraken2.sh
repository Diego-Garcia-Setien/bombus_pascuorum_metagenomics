#!/bin/bash
#SBATCH --job-name=bp_kraken2_tfm
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=120000
#SBATCH --array=1-3%3

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
INPUT_DIR="$WORKDIR/data/microbiome_reads"
OUTPUT_DATA="$WORKDIR/data/microbiota_taxonomy"
DATABASE="/data/lchueca/databases/kraken_std"

mkdir -p "$OUTPUT_DATA"
mkdir -p logs

###############################################################################
# Detect sample automatically
###############################################################################

R1=$(find "$INPUT_DIR" -maxdepth 1 -name "*_microbiome_R1.fastq.gz" \
      | sort \
      | sed -n "${SLURM_ARRAY_TASK_ID}p")

SAMPLE=$(basename "$R1")
SAMPLE=${SAMPLE%_microbiome_R1.fastq.gz}

R1="$INPUT_DIR/${SAMPLE}_microbiome_R1.fastq.gz"
R2="$INPUT_DIR/${SAMPLE}_microbiome_R2.fastq.gz"

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
      --output "$OUTPUT_DATA/${SAMPLE}.kraken2.out" \
      --report "$OUTPUT_DATA/${SAMPLE}.kraken2.report" \
      --gzip-compressed "$R1" "$R2"


echo "Clasificación taxonómica de $SAMPLE terminada"

# Vamos a usar Bracken, que es un programa complementario de Kraken2, 
# Sirve para estimar la abundancia en un solo nivel taxonómico


BRACKEN_DATA="$WORKDIR/data/bracken_taxonomy_results"

mkdir -p "$BRACKEN_DATA"

bracken -d "$DATABASE" -i "$OUTPUT_DATA/${SAMPLE}.kraken2.report" \
      -o "$BRACKEN_DATA/${SAMPLE}.bracken_output" -w "$BRACKEN_DATA/${SAMPLE}.bracken.kreport" -l S \
      -t "$CPU"

