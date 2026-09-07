######################################
# Script: visualize_coverm_glm.R
#
# Standalone script (separate from visualize_coverm.R): basic GLMs to
# look at MAG richness and Shannon diversity (from CoverM relative
# abundances, data/09.CoverMResults/coverm_relative_abundance.tsv) by
# Habitat and Year.
#
# Structure:
#   1) Setup
#   2) Load CoverM + compute MAG richness and Shannon diversity per sample
#   3) Load metadata + join (saved as metadata_diversity, also written
#      to visualize/b_pascuorum_samples_metadata_con_diversidad.csv)
#   4) GLM: Richness ~ Habitat (Poisson)
#   5) GLM: Richness ~ Habitat * Year (Poisson) + interaction plot
#   6) GLM: Shannon ~ Habitat (Gamma - Shannon isn't a count)
#   7) GLM: Shannon ~ Habitat * Year (Gamma) + interaction plot
######################################

## ---- 1) Setup ----------------------------------------------------------

required_packages <- c("tidyverse", "vegan")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(tidyverse)
library(vegan)
library(fitdistrplus)

# Adjust to your local path
coverm_file <- "visualize/coverm_relative_abundance.tsv"
metadata_file <- "visualize/b_pascuorum_samples_metadata_con_bodysize.csv"

## ---- 2) Load CoverM + compute MAG richness and Shannon diversity ----------

coverm_wide <- read_tsv(coverm_file, na = c("", "NaN"), show_col_types = FALSE)

# Column names look like "BPART241001.sorted Relative Abundance (%)";
# strip that down to just the sample ID.
names(coverm_wide) <- names(coverm_wide) |>
    str_remove("\\.sorted Relative Abundance \\(%\\)$")
names(coverm_wide)[1] <- "Genome"

mag_matrix <- coverm_wide |>
    filter(Genome != "unmapped") |>
    pivot_longer(-Genome, names_to = "Sample", values_to = "RelAbundance") |>
    mutate(RelAbundance = replace_na(RelAbundance, 0)) |>
    pivot_wider(names_from = Sample, values_from = RelAbundance) |>
    column_to_rownames("Genome") |>
    as.matrix()

sample_matrix <- t(mag_matrix)
sample_matrix[is.na(sample_matrix)] <- 0

# Richness = number of distinct MAGs (from the shared, dereplicated
# catalog) with nonzero relative abundance in a sample's own reads;
# Shannon = diversity index (vegan::diversity(), expects samples as
# rows) - same definitions as visualize_coverm.R.
diversity_df <- tibble(
    Sample = rownames(sample_matrix),
    Richness = rowSums(sample_matrix > 0),
    Shannon = diversity(sample_matrix, index = "shannon")
)

## ---- 3) Load metadata + join -----------------------------------------------

# Same conventions as visualize_coverm.R: Code lacks the "BP" prefix
# used in CoverM sample names, and a trailing "M" is the project's
# legacy metagenomics-vs-metatranscriptomics naming (same individual),
# stripped only for the join key.
metadata <- read_csv(metadata_file, show_col_types = FALSE) |>
    mutate(Sample = paste0("BP", Code))

# Year is numeric (2024/2025) as read from the CSV; only 2 distinct
# values exist, so it's used as a fixed-effect factor below (not a
# random effect/mixed model, which needs more levels than that to
# estimate a variance component reliably). Converted here, before the
# join, so richness_meta actually inherits it as a factor.
metadata$Year <- as.factor(metadata$Year)

strip_m_suffix <- function(x) str_remove(x, "M$")

# metadata_diversity IS the original metadata (every column from
# `metadata` comes along via left_join) with Richness and Shannon
# added - the "save into the original metadata" step. It only lives in
# this R session unless written out; see the write_csv() call below to
# persist it as an actual file.
metadata_diversity <- diversity_df |>
    mutate(Metadata_key = strip_m_suffix(Sample)) |>
    left_join(metadata, by = c("Metadata_key" = "Sample")) |>
    dplyr::select(-Metadata_key)

# Persist metadata_diversity to disk, so Richness/Shannon are saved
# alongside the original metadata columns as a real file, not just an
# in-session object. Comment out if you don't want this written.
write_csv(metadata_diversity, "visualize/b_pascuorum_samples_metadata_con_diversidad.csv")

## ---- 4) GLM: Richness ~ Habitat --------------------------------------------
#
descdist(metadata_diversity$Richness, discrete = TRUE)
# Richness is a count (number of MAGs detected) -> Poisson family
# (log link).

model_richness_habitat <- glm(Richness ~ Habitat, data = metadata_diversity, family = poisson(link = "log"))

cat("\nGLM summary (Richness ~ Habitat):\n")
print(summary(model_richness_habitat))

cat("\nAnalysis of deviance (sequential):\n")
print(anova(model_richness_habitat, test = "Chisq"))

## ---- 5) GLM: Richness ~ Habitat * Year -------------------------------------
#
# Interaction term: does the Habitat effect on Richness differ by
# Year (2024 vs. 2025)?

model_richness_habitat_year <- glm(Richness ~ Habitat * Year, data = metadata_diversity, family = poisson(link = "log"))

cat("\nGLM summary (Richness ~ Habitat * Year):\n")
print(summary(model_richness_habitat_year))

cat("\nAnalysis of deviance (sequential):\n")
print(anova(model_richness_habitat_year, test = "Chisq"))

# Interaction plot: model-predicted mean Richness per Habitat x Year
# combination, back-transformed from the log link (exp()), with a 95%
# CI. Roughly parallel lines = no interaction (consistent with the
# non-significant Habitat:Year term above); lines that cross or
# diverge sharply would indicate one.
#
# sort(unique(...)) rather than levels() for Habitat - it's a
# character column here, not a factor, and levels() only works on an
# actual factor (unlike Year, which already is one).
newdata <- expand.grid(
    Habitat = sort(unique(metadata_diversity$Habitat)),
    Year = levels(metadata_diversity$Year)
)

pred <- predict(model_richness_habitat_year, newdata = newdata, type = "link", se.fit = TRUE)
newdata$fit <- exp(pred$fit)
newdata$lower <- exp(pred$fit - 1.96 * pred$se.fit)
newdata$upper <- exp(pred$fit + 1.96 * pred$se.fit)

print(
    ggplot(newdata, aes(x = Habitat, y = fit, color = Year, group = Year)) +
        geom_point(position = position_dodge(width = 0.3), size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(width = 0.3)) +
        geom_line(position = position_dodge(width = 0.3)) +
        labs(title = "Predicted MAG richness by Habitat and Year", y = "Predicted Richness", x = "Habitat") +
        theme_minimal()
)

## ---- 6) GLM: Shannon ~ Habitat ----------------------------------------------
#
# Same idea as sections 4-5, but for Shannon diversity instead of
# Richness. Shannon is a continuous, strictly positive index - NOT a
# count, so Poisson doesn't apply to it (see the earlier chat
# discussion on this) - Gamma family (log link) is the right choice
# here, same as visualize_diversity_size.R.
#
# Gamma requires y > 0 strictly: 3 samples have Shannon == 0 (2 with 0
# MAGs detected at all, 1 with exactly 1 MAG - Shannon = 0 is the
# correct value for a single-taxon community, not a data error) and
# have to be excluded before fitting, same as visualize_coverm.R
# section 9 / visualize_diversity_size.R.
n_before <- nrow(metadata_diversity)
metadata_diversity_shannon <- metadata_diversity |> filter(Shannon > 0)
cat(sprintf(
    "\nExcluding %d / %d sample(s) with Shannon == 0 (required for the Gamma family).\n",
    n_before - nrow(metadata_diversity_shannon), n_before
))

descdist(metadata_diversity_shannon$Shannon, discrete = FALSE)

model_shannon_habitat <- glm(Shannon ~ Habitat, data = metadata_diversity_shannon, family = Gamma(link = "log"))

cat("\nGLM summary (Shannon ~ Habitat):\n")
print(summary(model_shannon_habitat))

cat("\nAnalysis of deviance (sequential, F-test - Gamma estimates its own dispersion):\n")
print(anova(model_shannon_habitat, test = "F"))

## ---- 7) GLM: Shannon ~ Habitat * Year --------------------------------------
#
# Interaction term: does the Habitat effect on Shannon diversity
# differ by Year (2024 vs. 2025)?

model_shannon_habitat_year <- glm(Shannon ~ Habitat * Year, data = metadata_diversity_shannon, family = Gamma(link = "log"))

cat("\nGLM summary (Shannon ~ Habitat * Year):\n")
print(summary(model_shannon_habitat_year))

cat("\nAnalysis of deviance (sequential, F-test):\n")
print(anova(model_shannon_habitat_year, test = "F"))

# Interaction plot - same recipe as section 5's, just Shannon instead
# of Richness (still exp() to undo the log link; still sort(unique())
# for Habitat since it's a character column, not a factor).
newdata_shannon <- expand.grid(
    Habitat = sort(unique(metadata_diversity_shannon$Habitat)),
    Year = levels(metadata_diversity_shannon$Year)
)

pred_shannon <- predict(model_shannon_habitat_year, newdata = newdata_shannon, type = "link", se.fit = TRUE)
newdata_shannon$fit <- exp(pred_shannon$fit)
newdata_shannon$lower <- exp(pred_shannon$fit - 1.96 * pred_shannon$se.fit)
newdata_shannon$upper <- exp(pred_shannon$fit + 1.96 * pred_shannon$se.fit)

print(
    ggplot(newdata_shannon, aes(x = Habitat, y = fit, color = Year, group = Year)) +
        geom_point(position = position_dodge(width = 0.3), size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(width = 0.3)) +
        geom_line(position = position_dodge(width = 0.3)) +
        labs(title = "Predicted Shannon diversity by Habitat and Year", y = "Predicted Shannon diversity", x = "Habitat") +
        theme_minimal()
)
