#!/bin/bash
#SBATCH --job-name=08_binning
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=32000
#SBATCH --array=1-93%25

set -euo pipefail

######################################
# Script: 08_binning.sh
#
# Binning de los contigs generados por 07_genome_assembly.sh
# (data/06.MegahitResults) usando tres programas distintos:
#   - MetaBAT2
#   - MaxBin2
#   - CONCOCT
#
# Para poder binear hace falta la cobertura de cada contig, así que
# primero se mapean las lecturas "microbioma" de cada muestra (las
# mismas usadas por megahit, ver 07_genome_assembly.sh) contra sus
# propios contigs con bowtie2. A partir de ese BAM se calcula:
#   - depth.txt (jgi_summarize_bam_contig_depths) -> usado por MetaBAT2
#     y, tras reformatear a 2 columnas, también por MaxBin2 (así no
#     hace falta mapear dos veces).
#   - la tabla de cobertura de subcontigs que necesita CONCOCT
#     (cut_up_fasta.py + concoct_coverage_table.py).
#
# Input:
#   data/06.MegahitResults/<muestra>/final.contigs.fa
#   data/03.MicrobiomeReads/<muestra>/<muestra>_microbiome_R{1,2}.fastq.gz
#
# Output:
#   data/07.BinningResults/<muestra>/
#       mapping/<muestra>.sorted.bam(.bai)
#       depth.txt
#       metabat2/
#       maxbin2/
#       concoct/
######################################

WORKDIR=$(pwd)

CONTIGS_DIR="$WORKDIR/data/06.MegahitResults"
READS_DIR="$WORKDIR/data/03.MicrobiomeReads"
OUTPUT_DIR="$WORKDIR/data/07.BinningResults"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

CPU=$SLURM_CPUS_PER_TASK

###############################################################################
# Detectar la muestra automáticamente (una subcarpeta por muestra en
# 06.MegahitResults, mismo esquema que 07_genome_assembly.sh)
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
R1=$(find "$READS_DIR/$SAMPLE" -maxdepth 1 -name "*_microbiome_R1.fastq.gz" | head -1 || true)
R2=$(find "$READS_DIR/$SAMPLE" -maxdepth 1 -name "*_microbiome_R2.fastq.gz" | head -1 || true)

if [[ ! -f "$CONTIGS" ]]; then
    echo "ERROR: Contigs not found for $SAMPLE"
    echo "$CONTIGS"
    exit 1
fi

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: FASTQ files not found for $SAMPLE"
    echo "$READS_DIR/$SAMPLE"
    exit 1
fi

SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE"
MAP_DIR="$SAMPLE_OUT/mapping"
METABAT_DIR="$SAMPLE_OUT/metabat2"
MAXBIN_DIR="$SAMPLE_OUT/maxbin2"
CONCOCT_DIR="$SAMPLE_OUT/concoct"

mkdir -p "$MAP_DIR" "$METABAT_DIR" "$MAXBIN_DIR" "$CONCOCT_DIR"

echo
echo "======================================"
echo "Binning de la muestra: $SAMPLE"
echo "======================================"
echo

###############################################################################
# 1) Mapear las lecturas de la muestra contra sus propios contigs
#    (necesario para calcular la cobertura que usan los 3 binners)
###############################################################################

module load Bowtie2/2.5.5-GCC-14.2.0
module load SAMtools/1.18-GCC-12.3.0

BAM="$MAP_DIR/${SAMPLE}.sorted.bam"

if [[ ! -f "$BAM" ]]; then
    INDEX="$MAP_DIR/${SAMPLE}_contigs"
    bowtie2-build --threads "$CPU" "$CONTIGS" "$INDEX"

    bowtie2 \
        --threads "$CPU" \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        2> "$MAP_DIR/${SAMPLE}.bowtie2.log" \
        | samtools view -@ "$CPU" -b - \
        | samtools sort -@ "$CPU" -o "$BAM"

    samtools index -@ "$CPU" "$BAM"
else
    echo "BAM ya existe para $SAMPLE, no se remapea."
fi

module purge

###############################################################################
# 2) MetaBAT2
###############################################################################

source /home/lchueca/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/lchueca/conda-env/metabat2

DEPTH_FILE="$SAMPLE_OUT/depth.txt"

jgi_summarize_bam_contig_depths --outputDepth "$DEPTH_FILE" "$BAM"

metabat2 \
    -i "$CONTIGS" \
    -a "$DEPTH_FILE" \
    -o "$METABAT_DIR/bin" \
    -t "$CPU"

conda deactivate

###############################################################################
# 3) MaxBin2
#    Reutiliza el depth.txt de MetaBAT2 (contigName + totalAvgDepth)
#    en vez de volver a mapear las lecturas.
###############################################################################

conda activate /scratch/lchueca/conda-env/maxbin2

MAXBIN_ABUND="$SAMPLE_OUT/maxbin_abundance.txt"
tail -n +2 "$DEPTH_FILE" | awk -F'\t' '{print $1"\t"$3}' > "$MAXBIN_ABUND"

run_MaxBin.pl \
    -contig "$CONTIGS" \
    -abund "$MAXBIN_ABUND" \
    -out "$MAXBIN_DIR/bin" \
    -thread "$CPU"

conda deactivate

###############################################################################
# 4) CONCOCT
###############################################################################

conda activate /scratch/lchueca/conda-env/concoct

CONTIGS_CUT="$CONCOCT_DIR/contigs_10K.fa"
BED_CUT="$CONCOCT_DIR/contigs_10K.bed"
COVERAGE_TABLE="$CONCOCT_DIR/coverage_table.tsv"

cut_up_fasta.py "$CONTIGS" -c 10000 -o 0 --merge_last -b "$BED_CUT" > "$CONTIGS_CUT"

concoct_coverage_table.py "$BED_CUT" "$BAM" > "$COVERAGE_TABLE"

concoct \
    --composition_file "$CONTIGS_CUT" \
    --coverage_file "$COVERAGE_TABLE" \
    -b "$CONCOCT_DIR/" \
    -t "$CPU"

merge_cutup_clustering.py "$CONCOCT_DIR/clustering_gt1000.csv" > "$CONCOCT_DIR/clustering_merged.csv"

mkdir -p "$CONCOCT_DIR/fasta_bins"
extract_fasta_bins.py "$CONTIGS" "$CONCOCT_DIR/clustering_merged.csv" --output_path "$CONCOCT_DIR/fasta_bins"

conda deactivate

echo
echo "Binning completado para $SAMPLE (MetaBAT2, MaxBin2, CONCOCT)"
echo
