#!/bin/bash

#SBATCH --job-name=bp_mapping
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntask-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=12000
#SBATCH --array=1-93%93

#################################
# Load software
#################################

module load Miniforge3/24.11.3-2

conda activate /scratch/lchueca/conda-env/fastp

CPU=8




WORKDIR=$(pwd)

IMPUT_DIR="$WORKDIR/data/fastp_results"
OUTPUT_DIR="$WORKDIR/data/bowtie2_results"
FAILED_DIR="$WORKDIR/data/bowtie2_failed"
GENOME= "WORKDIR/data/BP_GENOME"

mkdir -p ./data/bowtie2_results

mkdir -p ./data/bowtie2_failed

#Indexación de un genoma/ secuencia de referencia

cd "$GENOME"

$ bowtie2-build bombus_pascuorum_genome.fasta bp_index

#Alineación de un genoma/ secuancia indexada

$ bowtie2 --no-unal -p n -x bp_index -1 reads_1.fastq -2 reads_2.fastq -S output.sam

#-- no-unal es opcional, es para indicar que las lecturas que no se alineen con el genoma de referencia no se escribiram en la sam salida.

#-p es el número (n) de procesadores/hilos utilizados.

#-x es el índice del genoma.

#-1 es el/los archivos que contiene(n) lectura de pareja 1 
#-2 es el/los archivos que contiene(n) lectura de pareja 2

#indicamos que el formato de salida es .sam
