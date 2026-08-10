######################################
# Script: visualize_gtdbtk.R
#
# Exploratory visualization of the GTDB-Tk taxonomic classification
# produced by 10_gtdbTK.sh
# (data/08.GTDBtkResults/gtdbtk_output/gtdbtk.bac120.summary.tsv).
#
#
# Structure:
#   1) Setup
#   2) Load + tidy (split Sample/Bin, split taxonomy into ranks)
#   3) Per-sample summary (MAGs recovered + their taxonomy)
#   4) Overall composition (counts per rank, bar plots)
#   5) Metadata integration (TEMPLATE - same as visualize_coverm.R)
######################################

## ---- 1) Setup ----------------------------------------------------------

required_packages <- c("tidyverse")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(tidyverse)

# Adjust to your local path
gtdbtk_file <- "gtdbtk.bac120.summary.tsv"

# If you also have an archaeal run (gtdbtk.ar53.summary.tsv), read and
# bind_rows() it with the bacterial table below - same column layout.

## ---- 2) Load + tidy ------------------------------------------------------

# GTDB-Tk writes literal "N/A" (not a blank cell) when a value doesn't
# apply (e.g. msa_percent/closest_genome_ani/warnings for genomes
# classified by ANI alone). Without na = "N/A" these columns get read
# in as character instead of numeric, and warnings looks non-empty
# for every row.
gtdbtk_raw <- read_tsv(gtdbtk_file, na = c("", "N/A"), show_col_types = FALSE)

# user_genome looks like "BPGAS240901__BPGAS240901.MAGScoT_cleanbin_000001"
# (sample prefix added by 10_gtdbTK.sh so genomes from every sample stay
# unique in the combined classify_wf run). Split it back into Sample + Bin.
#
# classification looks like "d__Bacteria;p__Pseudomonadota;c__...;s__Genus species"
# Split on ";" into the 7 standard ranks and strip the "x__" prefixes.
gtdbtk <- gtdbtk_raw |>
    select(user_genome, classification, closest_genome_ani, msa_percent, warnings) |>
    separate(user_genome, into = c("Sample", "Bin"), sep = "__", extra = "merge") |>
    separate(classification,
        into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
        sep = ";", fill = "right"
    ) |>
    mutate(across(Domain:Species, ~ str_remove(., "^[a-z]__")))

## ---- 3) Per-sample summary ------------------------------------------------

# How many MAGs were recovered per sample, and what are they
mags_per_sample <- gtdbtk |>
    count(Sample, name = "n_MAGs") |>
    arrange(desc(n_MAGs))

print(mags_per_sample)

# Convenience: taxonomy of all MAGs for one specific sample
# (change "BPGAS240901" for any sample you want to look up)
gtdbtk |>
    filter(Sample == "BPGAS240901") |>
    select(Bin, Genus, Species, closest_genome_ani, msa_percent)

## ---- 4) Overall composition ----------------------------------------------

# 4a. How many MAGs recovered per sample (bar plot, all samples)
ggplot(mags_per_sample, aes(x = reorder(Sample, n_MAGs), y = n_MAGs)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "MAGs recovered per sample", x = NULL, y = "Number of MAGs") +
    theme_minimal(base_size = 7)

# 4b. Genus-level composition across the whole dataset (how many MAGs
#     of each genus were recovered, regardless of sample)
genus_counts <- gtdbtk |>
    count(Genus, sort = TRUE) |>
    mutate(Genus = replace_na(Genus, "Unclassified"))

ggplot(genus_counts, aes(x = reorder(Genus, n), y = n)) +
    geom_col(fill = "darkgreen") +
    coord_flip() +
    labs(title = "MAG count by genus (all samples combined)", x = NULL, y = "Number of MAGs") +
    theme_minimal(base_size = 9)

# 4c. Per-sample composition, stacked bar by genus (same style as the
#     CoverM abundance plot in visualize_coverm.R, but counting
#     MAGs instead of relative abundance)
ggplot(gtdbtk, aes(x = Sample, fill = Genus)) +
    geom_bar() +
    labs(title = "MAG taxonomic composition per sample (genus)", x = NULL, y = "Number of MAGs") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

# 4d. Quality sanity check: MSA percent and ANI to closest reference
#     genome (low msa_percent or low ANI are worth double-checking)
ggplot(gtdbtk, aes(x = msa_percent)) +
    geom_histogram(binwidth = 5, fill = "grey50") +
    labs(title = "MSA percent (marker-gene coverage) across all MAGs", x = "MSA %", y = "Count") +
    theme_minimal()

gtdbtk |>
    filter(!is.na(warnings)) |>
    select(Sample, Bin, Genus, Species, warnings) |>
    print(n = Inf)

## ---- 5) Metadata integration (TEMPLATE) ----------------------------------
#
# Same idea as in visualize_coverm.R: once the metadata table
# (locality, habitat_type, collection_year, host_size) is available,
# join it here on "Sample" to look at whether taxonomic composition
# (e.g. genus-level richness, presence of a given genus) relates to
# habitat/locality/altitude/host size.
#
# metadata <- read_csv("metadata.csv", show_col_types = FALSE)
#
# mags_meta <- mags_per_sample |> left_join(metadata, by = "Sample")
#
# ggplot(mags_meta, aes(x = habitat_type, y = n_MAGs, fill = habitat_type)) +
#     geom_boxplot() +
#     labs(title = "Number of MAGs recovered by habitat type") +
#     theme_minimal()
