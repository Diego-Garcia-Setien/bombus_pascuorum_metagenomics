#!/bin/bash

#SBATCH --job-name=04_build_host_index
#SBATCH --error=logs/%x-%j.err
#SBATCH --output=logs/%x-%j.out

#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=24000

###############################################################################
# Script: 04_build_host_index.sh
#
# Description:
#   Build the Bowtie2 index for the Bombus pascuorum reference genome.
#
# Input:
#   data/references/GCF_905332965.1_iyBomPasc1.1_genomic.fna
#
# Output:
#   data/references/BombusPasc*.bt2
#
###############################################################################

set -euo pipefail

#######################################
# Load software
#######################################

module load Bowtie2/2.5.5-GCC-14.2.0

#######################################
# Settings
#######################################

CPU=16

#######################################
# Directories
#######################################

WORKDIR=$(pwd)

REFERENCE="$WORKDIR/data/references/GCF_905332965.1_iyBomPasc1.1_genomic.fna"

INDEX_PREFIX="$WORKDIR/data/references/BombusPasc"

#######################################
# Check reference genome
#######################################

if [[ ! -f "$REFERENCE" ]]; then
    echo "ERROR: Reference genome not found:"
    echo "$REFERENCE"
    exit 1
fi

#######################################
# Build index
#######################################

echo
echo "======================================="
echo "Building Bowtie2 index"
echo "======================================="
echo

bowtie2-build \
    --threads "$CPU" \
    "$REFERENCE" \
    "$INDEX_PREFIX"

echo
echo "Bowtie2 index created successfully."
echo
