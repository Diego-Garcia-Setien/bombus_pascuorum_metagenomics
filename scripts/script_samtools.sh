#!/bin/bash
#SBATCH --job-name=bp_extract_unmapped
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=12000
#SBATCH --array=1-93%93

#################################
# Cargar software
#################################
module load Miniforge3/24.11.3-2
conda activate /scratch/lchueca/conda-env/samtools

CPU=8

######################################
# Script: script_samtools_extract_unmapped.sh
#
# Versión simple: coge el SAM de bowtie2 de cada muestra, se queda
# con las lecturas que NO alinean con el genoma del abejorro
# (candidatas a ser microbiota) y las guarda en fastq.
######################################

WORKDIR=$(pwd)
SAM_DIR="$WORKDIR/data/bowtie2_results"
OUTPUT_DIR="$WORKDIR/data/microbiota_reads"
SAMPLES_LIST="$WORKDIR/data/samples.txt"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Selecciona la muestra correspondiente a esta tarea del array
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")

SAM_IN="$SAM_DIR/${SAMPLE}.sam"
BAM_UNMAPPED_SORTED="$OUTPUT_DIR/${SAMPLE}.unmapped.sorted.bam"
FASTQ_R1="$OUTPUT_DIR/${SAMPLE}_microbiota_1.fastq.gz"
FASTQ_R2="$OUTPUT_DIR/${SAMPLE}_microbiota_2.fastq.gz"

echo "Procesando muestra: $SAMPLE"

# SAM -> BAM, quedándonos solo con lecturas donde ambas parejas no alinean (-f 12),
# y ordenado por nombre (-n) para poder separar R1/R2 después
samtools view -b -f 12 -@ "$CPU" "$SAM_IN" | samtools sort -n -@ "$CPU" -o "$BAM_UNMAPPED_SORTED"

# BAM -> fastq (R1 y R2)
samtools fastq -@ "$CPU" -1 "$FASTQ_R1" -2 "$FASTQ_R2" "$BAM_UNMAPPED_SORTED"

echo "Terminado: $SAMPLE"

