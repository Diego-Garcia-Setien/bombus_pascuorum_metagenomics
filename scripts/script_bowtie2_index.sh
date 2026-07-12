#!/bin/bash
#SBATCH --job-name=bp_index
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=12000

#################################
# Cargar software
#################################
module load Miniforge3/24.11.3-2
conda activate /scratch/lchueca/conda-env/bowtie2

######################################
# Script: script_bowtie2_index.sh
#
# 1) Construye el índice de bowtie2 a partir del genoma de referencia.
# 2) Genera "samples.txt", la lista de muestras que usará el array
#    de mapeo (script_bowtie2_align.sh)
# Se ejecuta UNA SOLA VEZ, antes de lanzar el array de mapeo.
######################################

WORKDIR=$(pwd)
GENOME_DIR="$WORKDIR/data/BP_GENOME"
GENOME_FASTA="bombus_pascuorum_genome.fasta"
INDEX_NAME="bp_index"
SAMPLES_LIST="$WORKDIR/data/samples.txt"
INPUT_DIR="$WORKDIR/data/fastp_results"

mkdir -p logs

######################################
# 1) Indexación del genoma
######################################

cd "$GENOME_DIR" || { echo "No se pudo acceder a $GENOME_DIR"; exit 1; }

# Solo construye el índice si no existe ya
if [ ! -f "${INDEX_NAME}.1.bt2" ]; then
    echo "Construyendo índice de bowtie2..."
    bowtie2-build "$GENOME_FASTA" "$INDEX_NAME"
else
    echo "El índice ya existe, no se reconstruye."
fi

######################################
# 2) Generar samples.txt para el array de mapeo
######################################
cd "$WORKDIR" || { echo "No se pudo volver a $WORKDIR"; exit 1; }

echo "Generando lista de muestras en $SAMPLES_LIST..."
ls "$INPUT_DIR"/*_1.fastq.gz | sed 's/_1.fastq.gz//; s#.*/##' > "$SAMPLES_LIST"

N_SAMPLES=$(wc -l < "$SAMPLES_LIST")
echo "Se han detectado $N_SAMPLES muestras."
echo "Recuerda ajustar '#SBATCH --array=1-$N_SAMPLES' en script_bowtie2_align.sh si el número no coincide."

