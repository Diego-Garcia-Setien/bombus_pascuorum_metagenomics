#!/bin/bash

#SBATCH --job-name=bp_multiqc_results_tfm
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=2000

#######################################
# Load software
#######################################

module load Miniforge3/24.11.3-2

conda activate /scratch/lchueca/conda-env/MultiQC

CPU=8

#Data

WORKDIR=$(pwd)

INPUT_DIR="$WORKDIR/data/fastp_results"
OUTPUT_DIR="$WORKDIR/data/multiqc_fastp_report"

#Creamos la carpeta para guardar el report de multiqc

mkdir -p "$WORKDIR"/data/multiqc_fastp_report

#Vamos a utilizar MultiQC, hay que trabajar con los archivos .json de fastp

multiqc "$INPUT_DIR" --outdir "$OUTPUT_DIR" --force  
