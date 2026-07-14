#!/bin/bash

#SBATCH --job-name=06_host_depletion
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --output=logs/%x-%A_%a.out

#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=32000
#SBATCH --array=1-3%3

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
#      bam/
#      bam/*.bai
#
#      host_reads/
#
#      microbiome_reads/
#
#      logs_alignment/
#
#      mapping_summary/
#
###############################################################################

set -euo pipefail

module load Bowtie2/2.5.5-GCC-14.2.0
module load SAMtools/1.18-GCC-12.3.0

CPU_BOWTIE=10
CPU_SAMTOOLS=6

WORKDIR=$(pwd)

INPUT_DIR="$WORKDIR/data/fastp_results"

INDEX="$WORKDIR/data/references/BombusPasc"

BAM_DIR="$WORKDIR/data/bam"

HOST_DIR="$WORKDIR/data/host_reads"

MICRO_DIR="$WORKDIR/data/microbiome_reads"

ALIGN_DIR="$WORKDIR/data/logs_alignment"

SUMMARY_DIR="$WORKDIR/data/mapping_summary"

mkdir -p "$BAM_DIR"
mkdir -p "$HOST_DIR"
mkdir -p "$MICRO_DIR"
mkdir -p "$ALIGN_DIR"
mkdir -p "$SUMMARY_DIR"

########################################
# Detect sample automatically
########################################

SAMPLE=$(find "$INPUT_DIR" -maxdepth 1 -name "*_R1.fq.gz" \
        | sort \
        | sed -n "${SLURM_ARRAY_TASK_ID}p")

SAMPLE=$(basename "$SAMPLE")
SAMPLE=${SAMPLE%_R1.fq.gz}

echo
echo "======================================"
echo "Processing sample: $SAMPLE"
echo "======================================"
echo

R1="$INPUT_DIR/${SAMPLE}_R1.fq.gz"
R2="$INPUT_DIR/${SAMPLE}_R2.fq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "FASTQ files not found"
    echo "$R1"
    echo "$R2"
    exit 1
fi

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
    2> "$ALIGN_DIR/${SAMPLE}.bowtie2.log" \
    | samtools view \
        -@ "$CPU_SAMTOOLS" \
        -b - \
    | samtools sort \
        -@ "$CPU_SAMTOOLS" \
        -o "$BAM_DIR/${SAMPLE}.sorted.bam"

########################################
# Index BAM
########################################

samtools index \
    -@ "$CPU_SAMTOOLS" \
    "$BAM_DIR/${SAMPLE}.sorted.bam"

########################################
# Collate (group mates together) before FASTQ extraction
# samtools fastq requires mates to be adjacent to pair them correctly.
# collate is much faster than a full name-sort for this purpose.
########################################

COLLATE_TMP="$BAM_DIR/${SAMPLE}.collate_tmp"

samtools collate \
    -@ "$CPU_SAMTOOLS" \
    -o "$BAM_DIR/${SAMPLE}.collated.bam" \
    "$BAM_DIR/${SAMPLE}.sorted.bam" \
    "$COLLATE_TMP"

########################################
# Host reads
# -G 12 = exclude pairs where BOTH mates are unmapped
#       = keep pairs where AT LEAST ONE mate mapped to host
########################################

samtools fastq \
    -@ "$CPU_SAMTOOLS" \
    -G 12 \
    "$BAM_DIR/${SAMPLE}.collated.bam" \
    -1 "$HOST_DIR/${SAMPLE}_host_R1.fastq.gz" \
    -2 "$HOST_DIR/${SAMPLE}_host_R2.fastq.gz" \
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
    "$BAM_DIR/${SAMPLE}.collated.bam" \
    -1 "$MICRO_DIR/${SAMPLE}_microbiome_R1.fastq.gz" \
    -2 "$MICRO_DIR/${SAMPLE}_microbiome_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n

########################################
# Clean up the collated BAM (it's a temporary reordering just for
# extraction; the coordinate-sorted+indexed BAM above is kept)
########################################

rm -f "$BAM_DIR/${SAMPLE}.collated.bam"

########################################
# Mapping summary
# Count reads directly from the output FASTQ files, so the summary
# always reflects what is actually in host_reads/ and microbiome_reads/.
########################################

TOTAL_READS=$(grep "reads; of these:" \
    "$ALIGN_DIR/${SAMPLE}.bowtie2.log" | awk '{print $1}')

HOST_PAIRS=$(( $(zcat "$HOST_DIR/${SAMPLE}_host_R1.fastq.gz" | wc -l) / 4 ))
MICRO_PAIRS=$(( $(zcat "$MICRO_DIR/${SAMPLE}_microbiome_R1.fastq.gz" | wc -l) / 4 ))

HOST_PERCENT=$(awk "BEGIN {printf \"%.2f\",100*$HOST_PAIRS/$TOTAL_READS}")
MICRO_PERCENT=$(awk "BEGIN {printf \"%.2f\",100*$MICRO_PAIRS/$TOTAL_READS}")

CHECK_SUM=$((HOST_PAIRS+MICRO_PAIRS))

echo -e "Sample\tTotal_pairs\tHost_pairs\tHost_percent\tMicrobiome_pairs\tMicrobiome_percent\tSum_check" \
> "$SUMMARY_DIR/${SAMPLE}.summary.tsv"

echo -e "${SAMPLE}\t${TOTAL_READS}\t${HOST_PAIRS}\t${HOST_PERCENT}\t${MICRO_PAIRS}\t${MICRO_PERCENT}\t${CHECK_SUM}" \
>> "$SUMMARY_DIR/${SAMPLE}.summary.tsv"

if [[ "$CHECK_SUM" -ne "$TOTAL_READS" ]]; then
    echo "WARNING: Host_pairs + Microbiome_pairs ($CHECK_SUM) does not equal Total_pairs ($TOTAL_READS) for $SAMPLE"
fi

echo
echo "Finished $SAMPLE"
echo
