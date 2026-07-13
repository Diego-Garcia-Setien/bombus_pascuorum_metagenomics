#!/bin/bash
#SBATCH --job-name=bp_mapping
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
conda activate /scratch/lchueca/conda-env/bowtie2

CPU=8

######################################
# Script: script_bowtie2_align.sh
#
# Alinea las lecturas pareadas de cada muestra contra el genoma
# ya indexado (ver script_bowtie2_index.sh).
#
# Requiere un archivo "samples.txt" en el directorio de trabajo,
# con un nombre de muestra por línea (93 líneas), sin extensión.
# Ejemplo de samples.txt:
# BP01
# BP02
# BP03
# BP04
#	....
#
# Se espera que los fastq se llamen:
#   $INPUT_DIR/<muestra>_1.fastq.gz
#   $INPUT_DIR/<muestra>_2.fastq.gz
######################################

WORKDIR=$(pwd)
INPUT_DIR="$WORKDIR/data/fastp_results"
OUTPUT_DIR="$WORKDIR/data/bowtie2_results"
FAILED_DIR="$WORKDIR/data/bowtie2_failed"
GENOME_DIR="$WORKDIR/data/BP_GENOME"
INDEX="$GENOME_DIR/bp_index"
SAMPLES_LIST="$WORKDIR/data/samples.txt"


mkdir -p "$OUTPUT_DIR"
mkdir -p "$FAILED_DIR"
mkdir -p logs

# Selecciona la muestra correspondiente a esta tarea del array
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: no se encontró muestra para la tarea $SLURM_ARRAY_TASK_ID en $SAMPLES_LIST"
    exit 1
fi

R1="$INPUT_DIR/${SAMPLE}_1.fastq.gz"
R2="$INPUT_DIR/${SAMPLE}_2.fastq.gz"
OUT_SAM="$OUTPUT_DIR/${SAMPLE}.sam"

echo "Tarea $SLURM_ARRAY_TASK_ID -> muestra: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: no se encuentran los archivos fastq de $SAMPLE"
    exit 1
fi

#-p     número de hilos a usar
#-x     índice del genoma de referencia
#-1/-2  archivos de lecturas pareadas 1 y 2
#-S     archivo de salida en formato SAM
#--no-unal descarta del SAM las lecturas que no alinean, no lo ponemos porque precisamente estas lecturas son las que nos interesan

bowtie2 -p "$CPU" -x "$INDEX" -1 "$R1" -2 "$R2" -S "$OUT_SAM"

if [ $? -eq 0 ]; then
    echo "Alineamiento completado correctamente para $SAMPLE"
else
    echo "ERROR en el alineamiento de $SAMPLE, moviendo a bowtie2_failed"
    mv "$OUT_SAM" "$FAILED_DIR/" 2>/dev/null
fi




