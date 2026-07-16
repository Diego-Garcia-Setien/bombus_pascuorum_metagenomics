#!/bin/bash

#SBATCH --job-name=05_host_depletion
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out

#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mem=40000
#SBATCH --array=1-93%25

###############################################################################
# Remove host reads while keeping BOTH host and microbiome reads
#
# FIX vs previous version:
#   - Old script used -f 2 (host) and -f 12 (microbiome). These two filters
#     are NOT complementary: discordant pairs and pairs with only one mate
#     mapped fell into neither file (a "dead zone" of lost reads).
#   - New script uses -G 12 (host) and -f 12 (microbiome), which ARE
#     complementary and exhaustive:
#       -f 12  -> both mates unmapped               => microbiome
#       -G 12  -> NOT (both mates unmapped), i.e.
#                 at least one mate mapped to host   => host
#     Every read pair goes to exactly one of the two files. This is also the
#     conservative standard used in metagenomic decontamination pipelines:
#     if either mate shows similarity to the host genome, the pair is
#     removed from the microbiome set.
#   - Mapping summary is now computed by directly counting reads in the
#     output FASTQ files, not by re-deriving numbers from the bowtie2 log
#     (the old derivation didn't match what ended up in the files).
#   - IMPORTANT: samtools fastq needs mates to be adjacent (name-sorted) to
#     pair them correctly. The BAM used for extraction is coordinate-sorted,
#     so discordant / mixed-mapping pairs can end up far apart in the file.
#     When that happens, samtools fastq treats them as orphan singletons and
#     silently sends them to -s (which we point at /dev/null) -> reads
#     vanish from both host_reads/ and microbiome_reads/ with no error.
#     Fix: run samtools collate on the sorted BAM to group mates together
#     before extracting FASTQ. The coordinate-sorted+indexed BAM is still
#     kept for other downstream uses; only the fastq extraction step uses
#     the collated version.
#
# Outputs:
#
#   data/
#      03.bamFiles/*.bai
#      03.HostReads
#      03.MicrobiomeReads
#
#      03_Mapping.Stats/logs_alignment/
#      03_Mapping.Stats/mapping_summary/
#
###############################################################################

set -euo pipefail

module load Bowtie2/2.5.5-GCC-14.2.0
module load SAMtools/1.18-GCC-12.3.0

CPU_BOWTIE=10
CPU_SAMTOOLS=6

WORKDIR=$(pwd)

INPUT_DIR="$WORKDIR/data/02.CleanReads"
INDEX="$WORKDIR/data/references/BombusPasc"
BAM_DIR="$WORKDIR/data/03.bamFiles"
HOST_DIR="$WORKDIR/data/03.HostReads"
MICRO_DIR="$WORKDIR/data/03.MicrobiomeReads"
ALIGN_DIR="$WORKDIR/data/03_Mapping.Stats/logs_alignment"
SUMMARY_DIR="$WORKDIR/data/03_Mapping.Stats/mapping_summary"

########################################
# Detect sample automatically (one subdirectory per sample,
# same layout as 01_quality_check.sh / 02_fastp.sh)
########################################

cd "$INPUT_DIR"

SAMPLE=$(find . -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=${SAMPLE#./}

if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: Sample not found."
    exit 1
fi

SAMPLE_DIR="$INPUT_DIR/$SAMPLE"

echo
echo "======================================"
echo "Processing sample: $SAMPLE"
echo "======================================"
echo

R1=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*_1.fq.gz" | head -1 || true)
R2=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*_2.fq.gz" | head -1 || true)

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: FASTQ files not found."
    echo "$SAMPLE_DIR"
    exit 1
fi

########################################
# Per-sample output directories
########################################

SAMPLE_BAM_DIR="$BAM_DIR/$SAMPLE"
SAMPLE_HOST_DIR="$HOST_DIR/$SAMPLE"
SAMPLE_MICRO_DIR="$MICRO_DIR/$SAMPLE"
SAMPLE_ALIGN_DIR="$ALIGN_DIR/$SAMPLE"
SAMPLE_SUMMARY_DIR="$SUMMARY_DIR/$SAMPLE"

mkdir -p "$SAMPLE_BAM_DIR"
mkdir -p "$SAMPLE_HOST_DIR"
mkdir -p "$SAMPLE_MICRO_DIR"
mkdir -p "$SAMPLE_ALIGN_DIR"
mkdir -p "$SAMPLE_SUMMARY_DIR"

########################################
# Alignment
########################################

bowtie2 \
    --very-sensitive \
    --threads "$CPU_BOWTIE" \
    --reorder \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    2> "$SAMPLE_ALIGN_DIR/${SAMPLE}.bowtie2.log" \
    | samtools view \
        -@ "$CPU_SAMTOOLS" \
        -b - \
    | samtools sort \
        -@ "$CPU_SAMTOOLS" \
        -o "$SAMPLE_BAM_DIR/${SAMPLE}.sorted.bam"

########################################
# Index BAM
########################################

samtools index \
    -@ "$CPU_SAMTOOLS" \
    "$SAMPLE_BAM_DIR/${SAMPLE}.sorted.bam"

########################################
# Collate (group mates together) before FASTQ extraction
# samtools fastq requires mates to be adjacent to pair them correctly.
# collate is much faster than a full name-sort for this purpose.
########################################

COLLATE_TMP="$SAMPLE_BAM_DIR/${SAMPLE}.collate_tmp"

samtools collate \
    -@ "$CPU_SAMTOOLS" \
    -o "$SAMPLE_BAM_DIR/${SAMPLE}.collated.bam" \
    "$SAMPLE_BAM_DIR/${SAMPLE}.sorted.bam" \
    "$COLLATE_TMP"

########################################
# Host reads
# -G 12 = exclude pairs where BOTH mates are unmapped
#       = keep pairs where AT LEAST ONE mate mapped to host
########################################

samtools fastq \
    -@ "$CPU_SAMTOOLS" \
    -G 12 \
    "$SAMPLE_BAM_DIR/${SAMPLE}.collated.bam" \
    -1 "$SAMPLE_HOST_DIR/${SAMPLE}_host_R1.fastq.gz" \
    -2 "$SAMPLE_HOST_DIR/${SAMPLE}_host_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n

########################################
# Microbiome reads
# -f 12 = keep pairs where BOTH mates are unmapped
# (exact complement of the host filter above)
########################################

samtools fastq \
    -@ "$CPU_SAMTOOLS" \
    -f 12 \
    "$SAMPLE_BAM_DIR/${SAMPLE}.collated.bam" \
    -1 "$SAMPLE_MICRO_DIR/${SAMPLE}_microbiome_R1.fastq.gz" \
    -2 "$SAMPLE_MICRO_DIR/${SAMPLE}_microbiome_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n

########################################
# Clean up the collated BAM (it's a temporary reordering just for
# extraction; the coordinate-sorted+indexed BAM above is kept)
########################################

rm -f "$SAMPLE_BAM_DIR/${SAMPLE}.collated.bam"

########################################
# Mapping summary
# Count reads directly from the output FASTQ files, so the summary
# always reflects what is actually in host_reads/ and microbiome_reads/.
########################################

TOTAL_READS=$(grep "reads; of these:" \
    "$SAMPLE_ALIGN_DIR/${SAMPLE}.bowtie2.log" | awk '{print $1}')

HOST_PAIRS=$(( $(zcat "$SAMPLE_HOST_DIR/${SAMPLE}_host_R1.fastq.gz" | wc -l) / 4 ))
MICRO_PAIRS=$(( $(zcat "$SAMPLE_MICRO_DIR/${SAMPLE}_microbiome_R1.fastq.gz" | wc -l) / 4 ))

HOST_PERCENT=$(awk "BEGIN {printf \"%.2f\",100*$HOST_PAIRS/$TOTAL_READS}")
MICRO_PERCENT=$(awk "BEGIN {printf \"%.2f\",100*$MICRO_PAIRS/$TOTAL_READS}")

CHECK_SUM=$((HOST_PAIRS+MICRO_PAIRS))

echo -e "Sample\tTotal_pairs\tHost_pairs\tHost_percent\tMicrobiome_pairs\tMicrobiome_percent\tSum_check" \
> "$SAMPLE_SUMMARY_DIR/${SAMPLE}.summary.tsv"

echo -e "${SAMPLE}\t${TOTAL_READS}\t${HOST_PAIRS}\t${HOST_PERCENT}\t${MICRO_PAIRS}\t${MICRO_PERCENT}\t${CHECK_SUM}" \
>> "$SAMPLE_SUMMARY_DIR/${SAMPLE}.summary.tsv"

if [[ "$CHECK_SUM" -ne "$TOTAL_READS" ]]; then
    echo "WARNING: Host_pairs + Microbiome_pairs ($CHECK_SUM) does not equal Total_pairs ($TOTAL_READS) for $SAMPLE"
fi

echo
echo "Finished $SAMPLE"
echo
