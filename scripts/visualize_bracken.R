######################################
# Script: visualize_bracken.R
#
# Exploratory visualization of the Bracken family/genus/species
# relative abundance table and the per-sample bacterial fraction
# produced by combine_bracken.sh
# (data/05.BrackenTaxonomy/combined_bracken_family_genus_species.tsv,
#  data/05.BrackenTaxonomy/bacteria_fraction_per_sample.tsv).
#
# Structure:
#   1) Setup
#   2) Load + split into family/genus/species tables
#   3) Per-sample summary + composition plots (genus, family)
#   4) Diversity and prevalence
#   5) Bacterial vs. non-bacterial fraction per sample
#   6) Metadata integration (TEMPLATE - same as the other visualize_*.R scripts)
######################################

## ---- 1) Setup ----------------------------------------------------------

required_packages <- c("tidyverse", "vegan", "pheatmap")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(tidyverse)
library(vegan)
library(pheatmap)

# Adjust to your local path
bracken_file <- "visualize/combined_bracken_family_genus_species.tsv"
bacteria_fraction_file <- "bacteria_fraction_per_sample.tsv"

## ---- 2) Load + split ------------------------------------------------------

bracken <- read_tsv(bracken_file, show_col_types = FALSE)

# Non-bacterial contamination filter: Kraken2/Bracken databases built
# against NCBI RefSeq/GTDB can pick up host or reagent-contaminant
# reads even after host-depletion (05_host_depletion.sh already
# removes most Bombus host DNA, but some non-bacterial reads can still
# slip through, e.g. human handling contamination). Checked against
# the full dataset (not just the top 15): "Homo" (genus) / "Hominidae"
# (family) / "Homo sapiens" (species) are the only clearly
# non-bacterial names present among the 2179 genera and 620 families
# recovered here. This list is NOT exhaustive by design (it's not
# derived from the actual domain-level lineage, just from what was
# observed) - re-check top_genus/top_family (section 3b/3e) after
# re-running combine_bracken.sh on new data and extend it if needed.
EXCLUDED_GENERA <- c("Homo")
EXCLUDED_FAMILIES <- c("Hominidae")

bracken <- bracken |>
    filter(
        !(Rank == "G" & Name %in% EXCLUDED_GENERA),
        !(Rank == "F" & Name %in% EXCLUDED_FAMILIES),
        !(Rank == "S" & str_starts(Name, paste0(EXCLUDED_GENERA, " ", collapse = "|")))
    )

family <- bracken |> filter(Rank == "F") |> select(Sample, Name, Percentage)
genus <- bracken |> filter(Rank == "G") |> select(Sample, Name, Percentage)
species <- bracken |> filter(Rank == "S") |> select(Sample, Name, Percentage)

# Wide matrices (Taxon x Sample), used for diversity/ordination below.
# values_fill = 0 because a missing Sample/Name row means that taxon
# had 0% in that sample (not NA).
family_matrix <- family |>
    pivot_wider(names_from = Sample, values_from = Percentage, values_fill = 0) |>
    column_to_rownames("Name") |>
    as.matrix()

genus_matrix <- genus |>
    pivot_wider(names_from = Sample, values_from = Percentage, values_fill = 0) |>
    column_to_rownames("Name") |>
    as.matrix()

species_matrix <- species |>
    pivot_wider(names_from = Sample, values_from = Percentage, values_fill = 0) |>
    column_to_rownames("Name") |>
    as.matrix()

## ---- 3) Per-sample summary + composition ----------------------------------

# Convenience: full family/genus/species breakdown for one specific sample
bracken |>
    filter(Sample == "BPZEB250704", Rank == "G") |>
    arrange(desc(Percentage)) |>
    select(Name, Percentage, Reads_clade)

# 3a. Richness: how many families/genera/species detected (>0%) per sample
richness <- tibble(
    Sample = colnames(genus_matrix),
    Family_richness = colSums(family_matrix > 0),
    Genus_richness = colSums(genus_matrix > 0),
    Species_richness = colSums(species_matrix > 0)
)
print(richness)

# 3b. Top genera overall (Homo already excluded in section 2 - check
# here for any other host/contaminant genus that might need adding to
# EXCLUDED_GENERA)
top_genus <- genus |>
    group_by(Name) |>
    summarise(mean_pct = mean(Percentage)) |>
    slice_max(mean_pct, n = 15)

print(top_genus)

# 3c. Stacked bar of genus composition per sample (top 15 + Other)
genus_grouped <- genus |>
    mutate(Name_grp = if_else(Name %in% top_genus$Name, Name, "Other"))

ggplot(genus_grouped, aes(x = Sample, y = Percentage, fill = Name_grp)) +
    geom_col() +
    labs(
        title = "Genus-level community composition per sample (top 15 + Other)",
        x = NULL, y = "Relative abundance (%)", fill = "Genus"
    ) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

# 3d. Heatmap of the top 30 genera across all samples
top30_genus <- genus |>
    group_by(Name) |>
    summarise(mean_pct = mean(Percentage)) |>
    slice_max(mean_pct, n = 30) |>
    pull(Name)

pheatmap(
    log1p(genus_matrix[top30_genus, ]),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 5,
    main = "Top 30 genera, relative abundance (log1p-scaled)"
)

# 3e. Same idea at family level (top 15 + Other, stacked bar)
top_family <- family |>
    group_by(Name) |>
    summarise(mean_pct = mean(Percentage)) |>
    slice_max(mean_pct, n = 15)

print(top_family)

family_grouped <- family |>
    mutate(Name_grp = if_else(Name %in% top_family$Name, Name, "Other"))

ggplot(family_grouped, aes(x = Sample, y = Percentage, fill = Name_grp)) +
    geom_col() +
    labs(
        title = "Family-level community composition per sample (top 15 + Other)",
        x = NULL, y = "Relative abundance (%)", fill = "Family"
    ) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

## ---- 4) Diversity and prevalence -------------------------------------------

diversity_df <- tibble(
    Sample = colnames(genus_matrix),
    Shannon_genus = diversity(t(genus_matrix), index = "shannon"),
    Simpson_genus = diversity(t(genus_matrix), index = "simpson")
)

ggplot(diversity_df, aes(x = reorder(Sample, Shannon_genus), y = Shannon_genus)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Shannon diversity (genus level) per sample", x = NULL) +
    theme_minimal(base_size = 7)

# Prevalence: in how many samples is each genus detected (>0%)?
prevalence_genus <- genus |>
    group_by(Name) |>
    summarise(
        n_samples_present = sum(Percentage > 0),
        pct_samples_present = 100 * n_samples_present / n_distinct(genus$Sample),
        mean_pct_when_present = mean(Percentage[Percentage > 0])
    ) |>
    arrange(desc(n_samples_present))

print(head(prevalence_genus, 20))

## ---- 5) Bacterial vs. non-bacterial fraction per sample -------------------
#
# From bacteria_fraction_per_sample.tsv: Cellular_organisms_pct is
# essentially "everything except unclassified/artificial reads",
# Bacteria_pct is the slice of that which Kraken2 assigned to the
# Bacteria domain, and NonBacteria_pct = Cellular_organisms_pct -
# Bacteria_pct (host DNA, Archaea, or anything else non-bacterial
# still counted as a cellular organism).

bacteria_fraction_file <- "visualize/bacteria_fraction_per_sample.tsv"
bacteria_fraction <- read_tsv(bacteria_fraction_file, show_col_types = FALSE)

print(bacteria_fraction |> arrange(desc(NonBacteria_pct)))

ggplot(bacteria_fraction, aes(x = reorder(Sample, NonBacteria_pct), y = NonBacteria_pct)) +
    geom_col(fill = "firebrick") +
    coord_flip() +
    labs(
        title = "Non-bacterial fraction per sample",
        subtitle = "Cellular_organisms_pct - Bacteria_pct (host DNA, Archaea, etc.)",
        x = NULL, y = "Non-bacteria (%)"
    ) +
    theme_minimal(base_size = 7)

# Stacked view: Bacteria vs. non-bacteria per sample, in one bar each
bacteria_fraction |>
    select(Sample, Bacteria_pct, NonBacteria_pct) |>
    pivot_longer(-Sample, names_to = "Fraction", values_to = "Percentage") |>
    ggplot(aes(x = Sample, y = Percentage, fill = Fraction)) +
    geom_col() +
    labs(title = "Bacterial vs. non-bacterial fraction per sample", x = NULL, y = "%") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

## ---- 6) Metadata integration (TEMPLATE) ----------------------------------
#
# Same idea as visualize_coverm.R / visualize_gtdbtk.R: once the
# metadata table (locality, habitat_type, collection_year, host_size)
# is available, join it on "Sample" to relate community composition,
# diversity, or the non-bacterial fraction to those variables.
#
# metadata <- read_csv("metadata.csv", show_col_types = FALSE)
#
# diversity_meta <- diversity_df |> left_join(metadata, by = "Sample")
#
# ggplot(diversity_meta, aes(x = habitat_type, y = Shannon_genus, fill = habitat_type)) +
#     geom_boxplot() +
#     labs(title = "Genus-level Shannon diversity by habitat type") +
#     theme_minimal()
