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
# Construye el índice de bowtie2 a partir del genoma de referencia.
# Se ejecuta UNA SOLA VEZ, antes de lanzar el array de mapeo.
######################################

WORKDIR=$(pwd)
GENOME_DIR="$WORKDIR/data/BP_GENOME"
GENOME_FASTA="bombus_pascuorum_genome.fasta"
INDEX_NAME="bp_index"

mkdir -p logs

cd "$GENOME_DIR" || { echo "No se pudo acceder a $GENOME_DIR"; exit 1; }

# Solo construye el índice si no existe ya
if [ ! -f "${INDEX_NAME}.1.bt2" ]; then
    echo "Construyendo índice de bowtie2..."
    bowtie2-build "$GENOME_FASTA" "$INDEX_NAME"
else
    echo "El índice ya existe, no se reconstruye."
fi

