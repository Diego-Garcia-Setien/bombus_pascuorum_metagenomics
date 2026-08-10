######################################
# Script: visualize_coverm.R
#
# Exploratory visualization of the MAG relative-abundance table
# produced by 11.2_coverm.sh (data/09.CoverMResults/coverm_relative_abundance.tsv).
#
#
# Structure:
#   1) Setup
#   2) Load + tidy the CoverM table
#   3) Community composition (stacked bars, heatmap, prevalence)
#   4) Diversity (Shannon/Simpson) and ordination (NMDS, Bray-Curtis)
#   5) Metadata integration (TEMPLATE - needs the real metadata file)
#
# Section 5 is a skeleton: it needs the metadata table (locality,
# habitat type, collection year, host size) to be finished. Once
# that file is available, join it in step 5 and re-run steps 3-4
# with color/facet by the relevant metadata column.
######################################

## ---- 1) Setup ----------------------------------------------------------

required_packages <- c("tidyverse", "vegan", "pheatmap", "ggrepel")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(tidyverse)
library(vegan)
library(pheatmap)

# Adjust to your local path
coverm_file <- "coverm_relative_abundance.tsv"

## ---- 2) Load + tidy the CoverM table -----------------------------------

coverm_wide <- read_tsv(coverm_file, na = c("", "NaN"), show_col_types = FALSE)

# Column names look like "BPART241001.sorted Relative Abundance (%)";
# strip that down to just the sample ID.
names(coverm_wide) <- names(coverm_wide) |>
    str_remove("\\.sorted Relative Abundance \\(%\\)$")
names(coverm_wide)[1] <- "Genome"

# Long format: one row per Genome x Sample
coverm_long <- coverm_wide |>
    pivot_longer(-Genome, names_to = "Sample", values_to = "RelAbundance") |>
    mutate(RelAbundance = replace_na(RelAbundance, 0))

# Split off the "unmapped" fraction (not a genome) as its own table
unmapped <- coverm_long |> filter(Genome == "unmapped")
mag_long <- coverm_long |> filter(Genome != "unmapped")

# Wide matrix of MAGs x samples (genomes as rows), used for
# ordination/diversity below
mag_matrix <- mag_long |>
    pivot_wider(names_from = Sample, values_from = RelAbundance) |>
    column_to_rownames("Genome") |>
    as.matrix()

## ---- 3) Community composition -------------------------------------------

# 3a. Mean mapped fraction across samples (how much of each sample's
#     reads are explained by the MAG catalog)
mean_unmapped <- mean(unmapped$RelAbundance, na.rm = TRUE)
cat(sprintf("Mean unmapped fraction across samples: %.1f%%\n", mean_unmapped))

ggplot(unmapped, aes(x = reorder(Sample, RelAbundance), y = RelAbundance)) +
    geom_col(fill = "grey60") +
    coord_flip() +
    labs(
        title = "Unmapped read percentage per sample",
        x = NULL, y = "Unmapped (%)"
    ) +
    theme_minimal(base_size = 8)

# 3b. Stacked bar plot of MAG composition per sample (top 15 MAGs by
#     mean relative abundance, rest grouped as "Other")
top_mags <- mag_long |>
    group_by(Genome) |>
    summarise(mean_abund = mean(RelAbundance)) |>
    slice_max(mean_abund, n = 15) |>
    pull(Genome)

mag_long_grouped <- mag_long |>
    mutate(Genome_grp = if_else(Genome %in% top_mags, Genome, "Other"))

ggplot(mag_long_grouped, aes(x = Sample, y = RelAbundance, fill = Genome_grp)) +
    geom_col() +
    labs(
        title = "MAG relative abundance per sample (top 15 + Other)",
        x = NULL, y = "Relative abundance (%)", fill = "MAG"
    ) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

# 3c. Heatmap of MAG relative abundance across samples
pheatmap(
    log1p(mag_matrix),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 5,
    main = "MAG relative abundance (log1p-scaled)"
)

# 3d. Prevalence: in how many samples is each MAG detected (>0%)?
prevalence <- mag_long |>
    group_by(Genome) |>
    summarise(
        n_samples_present = sum(RelAbundance > 0),
        pct_samples_present = 100 * n_samples_present / n_distinct(mag_long$Sample),
        mean_abund_when_present = mean(RelAbundance[RelAbundance > 0])
    ) |>
    arrange(desc(n_samples_present))

print(prevalence)

## ---- 4) Diversity and ordination ----------------------------------------

# Shannon/Simpson diversity per sample, from MAG relative abundances
# (vegan's diversity() expects samples as rows)
sample_matrix <- t(mag_matrix)
sample_matrix[is.na(sample_matrix)] <- 0

diversity_df <- tibble(
    Sample = rownames(sample_matrix),
    Shannon = diversity(sample_matrix, index = "shannon"),
    Simpson = diversity(sample_matrix, index = "simpson"),
    Richness = rowSums(sample_matrix > 0)
)

ggplot(diversity_df, aes(x = reorder(Sample, Shannon), y = Shannon)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Shannon diversity of MAG community per sample", x = NULL) +
    theme_minimal(base_size = 8)

# NMDS ordination (Bray-Curtis) of samples based on MAG composition.
# Samples with 0% reads mapped to any MAG in the catalog (all-zero
# row) have to be dropped first: Bray-Curtis is undefined between two
# empty samples and metaMDS() errors out otherwise. This can
# legitimately happen for samples whose own bins never passed
# MAGScoT's quality threshold (see 09_MAGScoT.sh), so their community
# isn't represented in the dereplicated MAG catalog at all.
zero_samples <- rownames(sample_matrix)[rowSums(sample_matrix) == 0]
if (length(zero_samples) > 0) {
    cat(sprintf(
        "Dropping %d sample(s) with 0%% reads mapped to any MAG (excluded from NMDS): %s\n",
        length(zero_samples), paste(zero_samples, collapse = ", ")
    ))
}
ordination_matrix <- sample_matrix[rowSums(sample_matrix) > 0, , drop = FALSE]

set.seed(1)
nmds <- metaMDS(ordination_matrix, distance = "bray", trymax = 100)

nmds_scores <- as_tibble(scores(nmds, display = "sites"), rownames = "Sample")

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, label = Sample)) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(size = 2, max.overlaps = 15) +
    labs(title = sprintf("NMDS of MAG community (stress = %.3f)", nmds$stress)) +
    theme_minimal()

## ---- 5) Metadata integration (TEMPLATE) ----------------------------------
#
# TODO once the metadata file is available:
#   - Expected: one row per sample, with a column that matches the
#     "Sample" values used above (e.g. "BPART241001").
#   - Expected columns: locality, habitat_type, collection_year, host_size
#     (adjust names below to match the real file).
#
# metadata <- read_csv("metadata.csv", show_col_types = FALSE)
#
# diversity_meta <- diversity_df |>
#     left_join(metadata, by = "Sample")
#
# ggplot(diversity_meta, aes(x = habitat_type, y = Shannon, fill = habitat_type)) +
#     geom_boxplot() +
#     labs(title = "MAG Shannon diversity by habitat type") +
#     theme_minimal()
#
# nmds_meta <- nmds_scores |>
#     left_join(metadata, by = "Sample")
#
# ggplot(nmds_meta, aes(x = NMDS1, y = NMDS2, color = habitat_type, shape = as.factor(collection_year))) +
#     geom_point(size = 3) +
#     labs(title = "NMDS of MAG community, coloured by habitat, shaped by year") +
#     theme_minimal()
#
# PERMANOVA: does community composition differ by habitat/locality/year?
# adonis2(sample_matrix ~ habitat_type + locality + collection_year,
#         data = metadata[match(rownames(sample_matrix), metadata$Sample), ],
#         method = "bray")
