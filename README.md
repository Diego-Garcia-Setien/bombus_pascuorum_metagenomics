# Bombus_pascuorum_metagenomics

A workflow for Genome Resolved Metagenomics from the gut of *Bombus pascuorum* bumblebee 🐝

Metagenomic profiling pipeline for gut microbiome samples of *Bombus pascuorum*. The pipeline performs raw read quality control, adapter/quality trimming, host read depletion, and taxonomic classification of the microbiome fraction.

## Repository structure

```
bombus_pascuorum_metagenomics/
├── data/
│   ├── 01.RawReads/            # Raw fastq files (one subdirectory per sample)
│   ├── 02.CleanReads/          # Trimmed reads (fastp output)
│   ├── 03.bamFiles/            # Alignment bam files
│   ├── 03.HostReads/           # Reads mapped to the host genome
│   ├── 03.MicrobiomeReads/     # Non-host reads (microbiome fraction)
│   ├── 03.Mapping.Stats/       # Host-mapping statistics
│   │   ├── logs_alignment/
│   │   └── mapping_summary/
│   ├── 06.MicrobiotaTaxonomy/  # Kraken2 taxonomic classification
│   ├── 06.BrackenTaxonomy/     # Bracken abundance estimation
│   ├── QC/
│   │   ├── 01.FastQC_MultiQC/  # QC report on raw reads
│   │   ├── 03.FastQC_MultiQC/  # QC report on clean reads
│   │   └── 03.Fastp_MultiQC/   # MultiQC report on fastp output
│   └── references/             # Host genome and Bowtie2 index
├── scripts/                    # Pipeline scripts (submitted via SLURM)
└── logs/                       # SLURM .out / .err files
```

## Pipeline overview

| Step | Script | Description | Input | Output |
|------|--------|-------------|-------|--------|
| 1 | `01_quality_check.sh` | Run FastQC on raw reads, then aggregate with MultiQC | `data/01.RawReads/` | `data/QC/01.FastQC_MultiQC/` |
| 2 | `02_fastp.sh` | Adapter and quality trimming with fastp | `data/01.RawReads/` | `data/02.CleanReads/` |
| 3 | `03_qulity_check_fastp.sh` | Run FastQC on clean reads and MultiQC on FastQC + fastp reports | `data/02.CleanReads/` | `data/QC/03_FastQC_MultiQC/`, `data/QC/03_Fastp_MultiQC/` |
| 4 | `04_build_host_index.sh` | Build Bowtie2 index of the *Bombus pascuorum* reference genome | Reference genome (FASTA) | `data/references/` |
| 5 | `05_host_depletion.sh` | Map clean reads to the host genome with Bowtie2 and split host vs. non-host reads | `data/02.CleanReads/` | `data/03.HostReads/`, `data/03.MicrobiomeReads/`, `data/03_Mapping.Stats/` |
| 6 | `06_kraken2.sh` | Taxonomic classification (Kraken2) and abundance estimation (Bracken) on the microbiome fraction | `data/03.HostReads/` | `data/06.MicrobiotaTaxonomy/`, `data/06.BrackenTaxonomy/` |

## Requirements

- FastQC
- MultiQC
- fastp
- Bowtie2
- Kraken2
- Bracken
- SLURM workload manager

## Usage

All scripts are designed to be submitted as SLURM jobs from the `scripts/` directory:

```bash
sbatch scripts/01_quality_check.sh
sbatch scripts/02_fastp.sh
sbatch scripts/03_qulity_check_fastp.sh
sbatch scripts/04_build_host_index.sh
sbatch scripts/05_host_depletion.sh
sbatch scripts/06_kraken2.sh
```

Job logs (`.out` / `.err`) are written to the `logs/` directory.

## Sample layout

Each sample directory under `data/01.RawReads/` (and correspondingly under `data/02.CleanReads/`) contains paired-end reads:

```
data/01.RawReads/<sample>/
├── <sample>_fw.fastq.gz
└── <sample>_rv.fastq.gz
```
