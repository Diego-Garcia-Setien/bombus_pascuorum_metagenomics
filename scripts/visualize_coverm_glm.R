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
#   8-11) Sections 4-7 duplicated for Provincia_grouped instead of
#      Habitat (Navarra folded into Gipuzkoa - see metadata$Provincia_grouped).
#      8 and 10 each split into .a (all samples) and .b (robustness
#      check excluding Shannon < 0.75 or Richness < 2.5); 10 also has
#      a .c (same as .b, with significance brackets/stars and
#      Spanish title/axis on the combined violin+GLM plot)
#   12) GLM: Bracken Shannon ~ Provincia_grouped - same as section 10,
#      but Shannon from Kraken2/Bracken read classification (raw reads
#      vs. reference DB) instead of the CoverM MAG catalog
######################################

## ---- 1) Setup ----------------------------------------------------------

required_packages <- c("tidyverse", "vegan", "emmeans", "ggsignif")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(tidyverse)
library(vegan)
library(fitdistrplus)
library(emmeans)

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

metadata <- metadata |>
  mutate(Provincia_grouped = if_else(Provincias == "Navarra", "Gipuzkoa", Provincias))

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

## ---- 8.a) GLM: Richness ~ Provincia_grouped (all samples) -------------------
#
# Same as section 4, but grouped by Provincia_grouped (Navarra folded
# into Gipuzkoa, added manually above) instead of Habitat.
descdist(metadata_diversity$Richness, discrete = TRUE)
# Richness is a count (number of MAGs detected) -> Poisson family
# (log link). (Same distribution check as section 4 - Richness itself
# doesn't change with the grouping variable, only the model does.)

model_richness_provincia <- glm(Richness ~ Provincia_grouped, data = metadata_diversity, family = poisson(link = "log"))

cat("\nGLM summary (Richness ~ Provincia_grouped):\n")
print(summary(model_richness_provincia))

cat("\nAnalysis of deviance (sequential):\n")
print(anova(model_richness_provincia, test = "Chisq"))

# Custom color per province (Araba = green, Bizkaia = pink, Gipuzkoa =
# dark purple) - reused for sections 8 and 10's plots.
provincia_colors_custom <- c("Araba" = "#2E8B57", "Bizkaia" = "#E75480", "Gipuzkoa" = "#4B0082")

# Predicted-mean plot: no Year here (this model doesn't include it),
# so just one point + 95% CI per Provincia, back-transformed from the
# log link.
newdata_richness_provincia <- data.frame(
    Provincia_grouped = sort(unique(metadata_diversity$Provincia_grouped))
)
pred_richness_provincia <- predict(model_richness_provincia, newdata = newdata_richness_provincia, type = "link", se.fit = TRUE)
newdata_richness_provincia$fit <- exp(pred_richness_provincia$fit)
newdata_richness_provincia$lower <- exp(pred_richness_provincia$fit - 1.96 * pred_richness_provincia$se.fit)
newdata_richness_provincia$upper <- exp(pred_richness_provincia$fit + 1.96 * pred_richness_provincia$se.fit)

print(
    ggplot(newdata_richness_provincia, aes(x = Provincia_grouped, y = fit, color = Provincia_grouped)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
        scale_color_manual(values = provincia_colors_custom) +
        labs(title = "Predicted MAG richness by Provincia", y = "Predicted Richness", x = "Provincia", color = "Provincia") +
        theme_minimal()
)

# Violin of the raw Richness values by Provincia.
print(
    ggplot(metadata_diversity, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        scale_fill_manual(values = provincia_colors_custom) +
        labs(title = "MAG richness by Provincia", y = "Richness", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

# Combined: raw-data violin + the GLM's predicted mean/CI overlaid
# (kept a neutral dark red/black for the model estimate, so it stands
# out against the colored violins instead of blending in).
print(
    ggplot(metadata_diversity, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.3) +
        scale_fill_manual(values = provincia_colors_custom) +
        geom_pointrange(
            data = newdata_richness_provincia,
            aes(x = Provincia_grouped, y = fit, ymin = lower, ymax = upper),
            color = "black", size = 0.8, inherit.aes = FALSE
        ) +
        labs(title = "MAG richness by Provincia (violin + GLM estimate)", y = "Richness", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

# Pairwise comparison between the 3 provinces, Tukey-adjusted.
# emmeans() rather than plain Tukey HSD (aov()'s TukeyHSD only works
# on lm/aov objects, not GLMs) - it computes the estimated marginal
# means correctly on the Poisson model's scale, back-transforms them
# to counts (type = "response"), and adjusts the 3 pairwise p-values
# for multiple comparisons using the Tukey method by default.
cat("\nPairwise comparison between provinces (Tukey-adjusted):\n")
print(emmeans(model_richness_provincia, pairwise ~ Provincia_grouped, type = "response"))

## ---- 8.b) GLM: Richness ~ Provincia_grouped (low-diversity samples excluded) ----
#
# Robustness check: same model as 8.a, but excluding samples with
# Shannon < 0.75 OR Richness < 2.5 - low-diversity/low-richness
# samples that could be noisy or unreliable (e.g. very shallow
# sequencing for that individual). All 5 samples this removes happen
# to be from Bizkaia - worth keeping in mind when reading this
# section, since it's specifically testing whether Bizkaia's effect
# was being driven by its most extreme low-diversity samples.
metadata_diversity_filtered <- metadata_diversity |>
    filter(Shannon >= 0.75, Richness >= 2.5)

cat(sprintf(
    "\nExcluding %d / %d sample(s) with Shannon < 0.75 or Richness < 2.5:\n",
    nrow(metadata_diversity) - nrow(metadata_diversity_filtered), nrow(metadata_diversity)
))
print(metadata_diversity |> filter(Shannon < 0.75 | Richness < 2.5) |> dplyr::select(Sample, Shannon, Richness, Provincia_grouped))

descdist(metadata_diversity_filtered$Richness, discrete = TRUE)

model_richness_provincia_filtered <- glm(Richness ~ Provincia_grouped, data = metadata_diversity_filtered, family = poisson(link = "log"))

cat("\nGLM summary (Richness ~ Provincia_grouped, filtered):\n")
print(summary(model_richness_provincia_filtered))

cat("\nAnalysis of deviance (sequential):\n")
print(anova(model_richness_provincia_filtered, test = "Chisq"))

newdata_richness_provincia_filtered <- data.frame(
    Provincia_grouped = sort(unique(metadata_diversity_filtered$Provincia_grouped))
)
pred_richness_provincia_filtered <- predict(model_richness_provincia_filtered, newdata = newdata_richness_provincia_filtered, type = "link", se.fit = TRUE)
newdata_richness_provincia_filtered$fit <- exp(pred_richness_provincia_filtered$fit)
newdata_richness_provincia_filtered$lower <- exp(pred_richness_provincia_filtered$fit - 1.96 * pred_richness_provincia_filtered$se.fit)
newdata_richness_provincia_filtered$upper <- exp(pred_richness_provincia_filtered$fit + 1.96 * pred_richness_provincia_filtered$se.fit)

print(
    ggplot(newdata_richness_provincia_filtered, aes(x = Provincia_grouped, y = fit, color = Provincia_grouped)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
        scale_color_manual(values = provincia_colors_custom) +
        labs(title = "Predicted MAG richness by Provincia (filtered)", y = "Predicted Richness", x = "Provincia", color = "Provincia") +
        theme_minimal()
)

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        scale_fill_manual(values = provincia_colors_custom) +
        labs(title = "MAG richness by Provincia (filtered)", y = "Richness", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.3) +
        scale_fill_manual(values = provincia_colors_custom) +
        geom_pointrange(
            data = newdata_richness_provincia_filtered,
            aes(x = Provincia_grouped, y = fit, ymin = lower, ymax = upper),
            color = "black", size = 0.8, inherit.aes = FALSE
        ) +
        labs(title = "MAG richness by Provincia (filtered, violin + GLM estimate)", y = "Richness", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

cat("\nPairwise comparison between provinces, filtered (Tukey-adjusted):\n")
print(emmeans(model_richness_provincia_filtered, pairwise ~ Provincia_grouped, type = "response"))

## ---- 8.c) Same as 8.b, with significance brackets + Spanish labels ---------
#
# Duplicate of all 3 plots from 8.b - reuses model_richness_provincia_filtered
# and metadata_diversity_filtered as-is (same model/data, no need to
# refit). Every plot here gets: (1) horizontal brackets between
# provinces with a significant Tukey-adjusted pairwise difference,
# with asterisks for the usual significance thresholds (*** < 0.001,
# ** < 0.01, * < 0.05); (2) title/y-axis label translated to Spanish.
# Mirrors 10.c exactly (see the comments there for the plotmath/
# single-line-title rationale) - p_to_stars() is already defined above.

# Single-line plotmath title (species name italicized), sized down so
# it fits the plot width without wrapping - same approach as 10.c's.
plot_title_richness_es <- bquote(
    "Riqueza de MAGs de la comunidad de bacterias del intestino de "*italic("Bombus pascuorum")*" entre provincias."
)

# Pairwise comparison (same numbers as 8.b's) - computed once here,
# shared by all 3 plots' brackets below.
emm_richness_provincia_c <- emmeans(model_richness_provincia_filtered, pairwise ~ Provincia_grouped, type = "response")
print(emm_richness_provincia_c)

# p-value -> significance stars, then keep only the pairs that reach
# p < 0.05 - those are the ones that get a bracket drawn. Defined
# here (first use) and reused as-is in section 10.c below.
p_to_stars <- function(p) {
    dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ ""
    )
}

sig_contrasts_richness_provincia <- as.data.frame(emm_richness_provincia_c$contrasts) |>
    mutate(stars = p_to_stars(p.value)) |>
    separate(contrast, into = c("group1", "group2"), sep = " / ") |>
    filter(stars != "")

# 8.c-i. Predicted-mean plot - brackets stacked above the tallest
# upper CI (this plot's own y-range is much narrower than the raw
# Richness values, so it needs its own y_position column).
sig_contrasts_pred_richness <- sig_contrasts_richness_provincia
sig_contrasts_pred_richness$y_position <- max(newdata_richness_provincia_filtered$upper) +
    seq_len(nrow(sig_contrasts_pred_richness)) * 0.3

print(
    ggplot(newdata_richness_provincia_filtered, aes(x = Provincia_grouped, y = fit, color = Provincia_grouped)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
        scale_color_manual(values = provincia_colors_custom) +
        ggsignif::geom_signif(
            data = sig_contrasts_pred_richness,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, color = "black"
        ) +
        labs(title = plot_title_richness_es, y = "Riqueza de MAGs", x = "Provincia", color = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

# 8.c-ii. Violin of the raw Richness values - brackets stacked above
# the tallest violin.
sig_contrasts_violin_richness <- sig_contrasts_richness_provincia
sig_contrasts_violin_richness$y_position <- max(metadata_diversity_filtered$Richness) +
    seq_len(nrow(sig_contrasts_violin_richness)) * 1.5

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        scale_fill_manual(values = provincia_colors_custom) +
        ggsignif::geom_signif(
            data = sig_contrasts_violin_richness,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01
        ) +
        labs(title = plot_title_richness_es, y = "Riqueza de MAGs", x = "Provincia", fill = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

# 8.c-iii. Combined: raw-data violin + GLM-estimate overlaid -
# brackets stacked above the tallest violin (same base as 8.c-ii's).
sig_contrasts_combined_richness <- sig_contrasts_violin_richness

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.3) +
        scale_fill_manual(values = provincia_colors_custom) +
        geom_pointrange(
            data = newdata_richness_provincia_filtered,
            aes(x = Provincia_grouped, y = fit, ymin = lower, ymax = upper),
            color = "black", size = 0.8, inherit.aes = FALSE
        ) +
        ggsignif::geom_signif(
            data = sig_contrasts_combined_richness,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01
        ) +
        labs(title = plot_title_richness_es, y = "Riqueza de MAGs", x = "Provincia", fill = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

## ---- 8.d) Same as 8.c-ii, jitter points shaped by Habitat ------------------
#
# Duplicate of 8.c-ii's violin plot (raw Richness values, brackets,
# Spanish title) - reuses sig_contrasts_violin_richness and
# plot_title_richness_es as-is. Only change: geom_jitter's points get
# a shape per Habitat instead of all being the same shape. Habitat has
# 3 actual values in the data - "Natural", "Agrícola", "Urbano" -
# mapped here to square/triangle/circle respectively. Defined here
# (first use) and reused as-is in section 10.d below.
#
# Agrícola (triangle) and Natural (square) are drawn at double the
# size of Urbano (circle), per request - a separate scale_size_manual
# mapped to the same Habitat variable, with the same name/labels as
# scale_shape_manual so ggplot merges both into a single legend
# instead of drawing two.
#
# The legend LABEL for "Natural" is displayed as "Seminatural" (the
# write-up's preferred term - see earlier discussion), via a labels
# function on both scales; the underlying data/column is untouched,
# still "Natural".
#
# Named vector's names (and the legend title string) are re-tagged as
# UTF-8 (Encoding<-, NOT enc2utf8()) for the same reason as
# habitat_colors in visualize_coverm.R: a literal accented string
# typed in the source under a non-UTF-8 session locale renders as
# garbled dots (e.g. "H..bitat") or fails to match a UTF-8-tagged
# string from a CSV, instead of erroring.
habitat_shapes <- c("Natural" = 15, "Agrícola" = 17, "Urbano" = 16)
Encoding(names(habitat_shapes)) <- "UTF-8"

habitat_sizes <- c("Natural" = 3, "Agrícola" = 3, "Urbano" = 1.5)
Encoding(names(habitat_sizes)) <- "UTF-8"

habitat_legend_title <- "Hábitat"
Encoding(habitat_legend_title) <- "UTF-8"

habitat_legend_labels <- function(x) if_else(x == "Natural", "Seminatural", x)

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Richness, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(aes(shape = Habitat, size = Habitat), width = 0.1, alpha = 0.6) +
        scale_fill_manual(values = provincia_colors_custom) +
        scale_shape_manual(values = habitat_shapes, name = habitat_legend_title, labels = habitat_legend_labels) +
        scale_size_manual(values = habitat_sizes, name = habitat_legend_title, labels = habitat_legend_labels) +
        ggsignif::geom_signif(
            data = sig_contrasts_violin_richness,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01
        ) +
        labs(title = plot_title_richness_es, y = "Riqueza de MAGs", x = "Provincia", fill = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

## ---- 9) GLM: Richness ~ Provincia_grouped * Year ---------------------------
#
# Interaction term: does the Provincia effect on Richness differ by
# Year (2024 vs. 2025)?

model_richness_provincia_year <- glm(Richness ~ Provincia_grouped * Year, data = metadata_diversity, family = poisson(link = "log"))

cat("\nGLM summary (Richness ~ Provincia_grouped * Year):\n")
print(summary(model_richness_provincia_year))

cat("\nAnalysis of deviance (sequential):\n")
print(anova(model_richness_provincia_year, test = "Chisq"))

# Interaction plot - same recipe as section 5's, just Provincia_grouped
# instead of Habitat.
newdata_provincia <- expand.grid(
    Provincia_grouped = sort(unique(metadata_diversity$Provincia_grouped)),
    Year = levels(metadata_diversity$Year)
)

pred_provincia <- predict(model_richness_provincia_year, newdata = newdata_provincia, type = "link", se.fit = TRUE)
newdata_provincia$fit <- exp(pred_provincia$fit)
newdata_provincia$lower <- exp(pred_provincia$fit - 1.96 * pred_provincia$se.fit)
newdata_provincia$upper <- exp(pred_provincia$fit + 1.96 * pred_provincia$se.fit)

print(
    ggplot(newdata_provincia, aes(x = Provincia_grouped, y = fit, color = Year, group = Year)) +
        geom_point(position = position_dodge(width = 0.3), size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(width = 0.3)) +
        geom_line(position = position_dodge(width = 0.3)) +
        labs(title = "Predicted MAG richness by Provincia and Year", y = "Predicted Richness", x = "Provincia") +
        theme_minimal()
)

## ---- 10.a) GLM: Shannon ~ Provincia_grouped (all samples) -------------------
#
# Same as section 6, but grouped by Provincia_grouped instead of
# Habitat. Reuses metadata_diversity_shannon from section 6 - the
# Shannon == 0 exclusion doesn't depend on the grouping variable.
descdist(metadata_diversity_shannon$Shannon, discrete = FALSE)

model_shannon_provincia <- glm(Shannon ~ Provincia_grouped, data = metadata_diversity_shannon, family = Gamma(link = "log"))

cat("\nGLM summary (Shannon ~ Provincia_grouped):\n")
print(summary(model_shannon_provincia))

cat("\nAnalysis of deviance (sequential, F-test - Gamma estimates its own dispersion):\n")
print(anova(model_shannon_provincia, test = "F"))

# Predicted-mean plot: no Year here (this model doesn't include it),
# so just one point + 95% CI per Provincia, back-transformed from the
# log link. Named "_only" to avoid clashing with section 11's
# newdata_shannon_provincia (that one includes Year).
newdata_shannon_provincia_only <- data.frame(
    Provincia_grouped = sort(unique(metadata_diversity_shannon$Provincia_grouped))
)
pred_shannon_provincia_only <- predict(model_shannon_provincia, newdata = newdata_shannon_provincia_only, type = "link", se.fit = TRUE)
newdata_shannon_provincia_only$fit <- exp(pred_shannon_provincia_only$fit)
newdata_shannon_provincia_only$lower <- exp(pred_shannon_provincia_only$fit - 1.96 * pred_shannon_provincia_only$se.fit)
newdata_shannon_provincia_only$upper <- exp(pred_shannon_provincia_only$fit + 1.96 * pred_shannon_provincia_only$se.fit)

print(
    ggplot(newdata_shannon_provincia_only, aes(x = Provincia_grouped, y = fit, color = Provincia_grouped)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
        scale_color_manual(values = provincia_colors_custom) +
        labs(title = "Predicted Shannon diversity by Provincia", y = "Predicted Shannon diversity", x = "Provincia", color = "Provincia") +
        theme_minimal()
)

# Violin of the raw Shannon values by Provincia.
print(
    ggplot(metadata_diversity_shannon, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        scale_fill_manual(values = provincia_colors_custom) +
        labs(title = "Shannon diversity by Provincia", y = "Shannon diversity", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

# Combined: raw-data violin + the GLM's predicted mean/CI overlaid
# (kept a neutral black for the model estimate, so it stands out
# against the colored violins instead of blending in).
print(
    ggplot(metadata_diversity_shannon, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.3) +
        scale_fill_manual(values = provincia_colors_custom) +
        geom_pointrange(
            data = newdata_shannon_provincia_only,
            aes(x = Provincia_grouped, y = fit, ymin = lower, ymax = upper),
            color = "black", size = 0.8, inherit.aes = FALSE
        ) +
        labs(title = "Shannon diversity by Provincia (violin + GLM estimate)", y = "Shannon diversity", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

# Pairwise comparison between the 3 provinces, Tukey-adjusted - same
# approach as section 8's, applied to the Gamma model.
cat("\nPairwise comparison between provinces (Tukey-adjusted):\n")
print(emmeans(model_shannon_provincia, pairwise ~ Provincia_grouped, type = "response"))

## ---- 10.b) GLM: Shannon ~ Provincia_grouped (low-diversity samples excluded) ----
#
# Robustness check: same model as 10.a, but on metadata_diversity_filtered
# (from section 8.b - Shannon < 0.75 or Richness < 2.5 excluded, all 5
# such samples happen to be from Bizkaia). Note this is a stricter
# filter than 10.a's own Shannon > 0 requirement, so it's a subset of
# metadata_diversity_shannon, not built from it separately.

descdist(metadata_diversity_filtered$Shannon, discrete = FALSE)

model_shannon_provincia_filtered <- glm(Shannon ~ Provincia_grouped, data = metadata_diversity_filtered, family = Gamma(link = "log"))

cat("\nGLM summary (Shannon ~ Provincia_grouped, filtered):\n")
print(summary(model_shannon_provincia_filtered))

cat("\nAnalysis of deviance (sequential, F-test - Gamma estimates its own dispersion):\n")
print(anova(model_shannon_provincia_filtered, test = "F"))

newdata_shannon_provincia_filtered <- data.frame(
    Provincia_grouped = sort(unique(metadata_diversity_filtered$Provincia_grouped))
)
pred_shannon_provincia_filtered <- predict(model_shannon_provincia_filtered, newdata = newdata_shannon_provincia_filtered, type = "link", se.fit = TRUE)
newdata_shannon_provincia_filtered$fit <- exp(pred_shannon_provincia_filtered$fit)
newdata_shannon_provincia_filtered$lower <- exp(pred_shannon_provincia_filtered$fit - 1.96 * pred_shannon_provincia_filtered$se.fit)
newdata_shannon_provincia_filtered$upper <- exp(pred_shannon_provincia_filtered$fit + 1.96 * pred_shannon_provincia_filtered$se.fit)

print(
    ggplot(newdata_shannon_provincia_filtered, aes(x = Provincia_grouped, y = fit, color = Provincia_grouped)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
        scale_color_manual(values = provincia_colors_custom) +
        labs(title = "Predicted Shannon diversity by Provincia (filtered)", y = "Predicted Shannon diversity", x = "Provincia", color = "Provincia") +
        theme_minimal()
)

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        scale_fill_manual(values = provincia_colors_custom) +
        labs(title = "Shannon diversity by Provincia (filtered)", y = "Shannon diversity", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.3) +
        scale_fill_manual(values = provincia_colors_custom) +
        geom_pointrange(
            data = newdata_shannon_provincia_filtered,
            aes(x = Provincia_grouped, y = fit, ymin = lower, ymax = upper),
            color = "black", size = 0.8, inherit.aes = FALSE
        ) +
        labs(title = "Shannon diversity by Provincia (filtered, violin + GLM estimate)", y = "Shannon diversity", x = "Provincia", fill = "Provincia") +
        theme_minimal()
)

cat("\nPairwise comparison between provinces, filtered (Tukey-adjusted):\n")
print(emmeans(model_shannon_provincia_filtered, pairwise ~ Provincia_grouped, type = "response"))

## ---- 10.c) Same as 10.b, with significance brackets + Spanish labels -------
#
# Duplicate of all 3 plots from 10.b - reuses model_shannon_provincia_filtered
# and metadata_diversity_filtered as-is (same model/data, no need to
# refit). Every plot here gets: (1) horizontal brackets between
# provinces with a significant Tukey-adjusted pairwise difference,
# with asterisks for the usual significance thresholds (*** < 0.001,
# ** < 0.01, * < 0.05); (2) title/y-axis label translated to Spanish.

# A plotmath expression (bquote()) instead of a plain string, so
# "Bombus pascuorum" renders in italics (species names conventionally
# are) while the rest of the title stays upright. ggtext (which would
# let this be plain markdown, e.g. "*Bombus pascuorum*") isn't
# installed, so this uses base R's plotmath instead - no extra
# package needed. "*" between the quoted strings and italic(...) means
# "concatenate", not multiplication, in this context. Kept on a single
# line (no atop()) per request - element_text(size = ...) below is
# shrunk so it still fits the plot width instead of wrapping/overflowing.
# element_text(hjust = 0.5) centers it (left-aligned otherwise).
plot_title_es <- bquote(
    "Diversidad de Shannon de la comunidad de bacterias del intestino de "*italic("Bombus pascuorum")*" entre provincias."
)

# Pairwise comparison (same numbers as 10.b's) - computed once here,
# shared by all 3 plots' brackets below.
emm_shannon_provincia_c <- emmeans(model_shannon_provincia_filtered, pairwise ~ Provincia_grouped, type = "response")
print(emm_shannon_provincia_c)

# p_to_stars() (p-value -> significance stars) is already defined
# above, in section 8.c - reused here as-is.

sig_contrasts_shannon_provincia <- as.data.frame(emm_shannon_provincia_c$contrasts) |>
    mutate(stars = p_to_stars(p.value)) |>
    separate(contrast, into = c("group1", "group2"), sep = " / ") |>
    filter(stars != "")

# 10.c-i. Predicted-mean plot - brackets stacked above the tallest
# upper CI (this plot's own y-range is much narrower than the raw
# Shannon values, so it needs its own y_position column).
sig_contrasts_pred <- sig_contrasts_shannon_provincia
sig_contrasts_pred$y_position <- max(newdata_shannon_provincia_filtered$upper) +
    seq_len(nrow(sig_contrasts_pred)) * 0.05

print(
    ggplot(newdata_shannon_provincia_filtered, aes(x = Provincia_grouped, y = fit, color = Provincia_grouped)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
        scale_color_manual(values = provincia_colors_custom) +
        ggsignif::geom_signif(
            data = sig_contrasts_pred,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, color = "black"
        ) +
        labs(title = plot_title_es, y = "Diversidad de Shannon", x = "Provincia", color = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

# 10.c-ii. Violin of the raw Shannon values - brackets stacked above
# the tallest violin.
sig_contrasts_violin <- sig_contrasts_shannon_provincia
sig_contrasts_violin$y_position <- max(metadata_diversity_filtered$Shannon) +
    seq_len(nrow(sig_contrasts_violin)) * 0.15

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        scale_fill_manual(values = provincia_colors_custom) +
        ggsignif::geom_signif(
            data = sig_contrasts_violin,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01
        ) +
        labs(title = plot_title_es, y = "Diversidad de Shannon", x = "Provincia", fill = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

# 10.c-iii. Combined: raw-data violin + GLM-estimate overlaid -
# brackets stacked above the tallest violin (same base as 10.c-ii's).
sig_contrasts_combined <- sig_contrasts_violin

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.3) +
        scale_fill_manual(values = provincia_colors_custom) +
        geom_pointrange(
            data = newdata_shannon_provincia_filtered,
            aes(x = Provincia_grouped, y = fit, ymin = lower, ymax = upper),
            color = "black", size = 0.8, inherit.aes = FALSE
        ) +
        ggsignif::geom_signif(
            data = sig_contrasts_combined,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01
        ) +
        labs(title = plot_title_es, y = "Diversidad de Shannon", x = "Provincia", fill = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

## ---- 10.d) Same as 10.c-ii, jitter points shaped by Habitat ----------------
#
# Duplicate of 10.c-ii's violin plot (raw Shannon values, brackets,
# Spanish title) - reuses sig_contrasts_violin and plot_title_es as-is.
# Only change: geom_jitter's points get a shape per Habitat instead of
# all being the same shape, so a reader can see how the 3 habitats are
# distributed within each province's violin, with Agrícola/Natural
# drawn at double the size of Urbano and "Natural" displayed as
# "Seminatural" in the legend. habitat_shapes, habitat_sizes,
# habitat_legend_title and habitat_legend_labels are all already
# defined above, in section 8.d (same rationale there) - reused here
# as-is.

print(
    ggplot(metadata_diversity_filtered, aes(x = Provincia_grouped, y = Shannon, fill = Provincia_grouped)) +
        geom_violin(alpha = 0.5) +
        geom_jitter(aes(shape = Habitat, size = Habitat), width = 0.1, alpha = 0.6) +
        scale_fill_manual(values = provincia_colors_custom) +
        scale_shape_manual(values = habitat_shapes, name = habitat_legend_title, labels = habitat_legend_labels) +
        scale_size_manual(values = habitat_sizes, name = habitat_legend_title, labels = habitat_legend_labels) +
        ggsignif::geom_signif(
            data = sig_contrasts_violin,
            aes(xmin = group1, xmax = group2, y_position = y_position, annotations = stars),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01
        ) +
        labs(title = plot_title_es, y = "Diversidad de Shannon", x = "Provincia", fill = "Provincia") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
)

## ---- 11) GLM: Shannon ~ Provincia_grouped * Year ---------------------------
#
# Interaction term: does the Provincia effect on Shannon diversity
# differ by Year (2024 vs. 2025)?

model_shannon_provincia_year <- glm(Shannon ~ Provincia_grouped * Year, data = metadata_diversity_shannon, family = Gamma(link = "log"))

cat("\nGLM summary (Shannon ~ Provincia_grouped * Year):\n")
print(summary(model_shannon_provincia_year))

cat("\nAnalysis of deviance (sequential, F-test):\n")
print(anova(model_shannon_provincia_year, test = "F"))

# Interaction plot - same recipe as section 7's, just Provincia_grouped
# instead of Habitat.
newdata_shannon_provincia <- expand.grid(
    Provincia_grouped = sort(unique(metadata_diversity_shannon$Provincia_grouped)),
    Year = levels(metadata_diversity_shannon$Year)
)

pred_shannon_provincia <- predict(model_shannon_provincia_year, newdata = newdata_shannon_provincia, type = "link", se.fit = TRUE)
newdata_shannon_provincia$fit <- exp(pred_shannon_provincia$fit)
newdata_shannon_provincia$lower <- exp(pred_shannon_provincia$fit - 1.96 * pred_shannon_provincia$se.fit)
newdata_shannon_provincia$upper <- exp(pred_shannon_provincia$fit + 1.96 * pred_shannon_provincia$se.fit)

print(
    ggplot(newdata_shannon_provincia, aes(x = Provincia_grouped, y = fit, color = Year, group = Year)) +
        geom_point(position = position_dodge(width = 0.3), size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(width = 0.3)) +
        geom_line(position = position_dodge(width = 0.3)) +
        labs(title = "Predicted Shannon diversity by Provincia and Year", y = "Predicted Shannon diversity", x = "Provincia") +
        theme_minimal()
)

print(
    ggplot(metadata_diversity_shannon, aes(x = Provincia_grouped, y = Shannon, fill = Year)) +
        geom_violin(alpha = 0.5, position = position_dodge(width = 0.8)) +
        geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1), alpha = 0.4, size = 1) +
        labs(title = "Shannon diversity by Provincia and Year", y = "Shannon diversity", x = "Provincia") +
        theme_minimal()
)

print(
    ggplot(metadata_diversity_shannon, aes(x = Provincia_grouped, y = Shannon)) +
        geom_violin(aes(fill = Year), alpha = 0.4, position = position_dodge(width = 0.8)) +
        geom_pointrange(
          data = newdata_shannon_provincia,
          aes(y = fit, ymin = lower, ymax = upper, color = Year),
          position = position_dodge(width = 0.8), size = 0.8
    ) +
    labs(title = "Shannon diversity by Provincia and Year (violin + GLM estimate)", y = "Shannon diversity", x = "Provincia") +
    theme_minimal()
)

## ---- 12) GLM: Bracken Shannon ~ Provincia_grouped (raw reads vs. reference DB) ----
#
# Same question as section 10 (Shannon ~ Provincia_grouped, Gamma
# family), but Shannon computed from Kraken2/Bracken read
# classification (genus level) instead of the CoverM MAG catalog -
# i.e. every read classified against the full reference database,
# with no assembly/binning step in between. Same data source as
# visualize_diversity_size_bracken.R / visualize_bracken.R; loaded
# independently here so this script doesn't depend on sourcing those.
#
# This is a genuinely different (not redundant) view from sections
# 4-11: the MAG catalog only reflects organisms that assembled/binned
# well enough to pass QC, so it can miss low-abundance taxa that
# Bracken still picks up directly from the reads.

bracken_file <- "visualize/combined_bracken_family_genus_species.tsv"

# Same non-bacterial contamination filter as visualize_bracken.R -
# Kraken2/Bracken databases can pick up host/reagent-contaminant reads
# even after host-depletion.
bracken_raw <- read_tsv(bracken_file, show_col_types = FALSE) |>
    filter(
        !(Rank == "G" & Name == "Homo"),
        !(Rank == "F" & Name == "Hominidae"),
        !(Rank == "S" & str_starts(Name, "Homo "))
    )

genus_matrix_bracken <- bracken_raw |>
    filter(Rank == "G") |>
    dplyr::select(Sample, Name, Percentage) |>
    pivot_wider(names_from = Sample, values_from = Percentage, values_fill = 0) |>
    column_to_rownames("Name") |>
    as.matrix()

# Bracken's Sample column already carries the "BP" prefix (unlike the
# raw metadata Code, but same as CoverM's Genome column) - no prefix
# fix needed here, just the usual "M" suffix join-key handling.
diversity_df_bracken <- tibble(
    Sample = colnames(genus_matrix_bracken),
    Shannon_bracken = diversity(t(genus_matrix_bracken), index = "shannon")
)

metadata_diversity_bracken <- diversity_df_bracken |>
    mutate(Metadata_key = strip_m_suffix(Sample)) |>
    left_join(metadata, by = c("Metadata_key" = "Sample")) |>
    dplyr::select(-Metadata_key)

# No Shannon == 0 samples expected here (genus-level Bracken richness
# runs in the hundreds per sample, unlike the sparse MAG catalog), but
# checked rather than assumed, same Gamma requirement as section 10.
n_before <- nrow(metadata_diversity_bracken)
metadata_diversity_bracken <- metadata_diversity_bracken |> filter(Shannon_bracken > 0)
cat(sprintf(
    "\nExcluding %d / %d sample(s) with Shannon_bracken == 0 (required for the Gamma family).\n",
    n_before - nrow(metadata_diversity_bracken), n_before
))

descdist(metadata_diversity_bracken$Shannon_bracken, discrete = FALSE)

model_shannon_bracken_provincia <- glm(
    Shannon_bracken ~ Provincia_grouped,
    data = metadata_diversity_bracken, family = Gamma(link = "log")
)

cat("\nGLM summary (Shannon_bracken ~ Provincia_grouped):\n")
print(summary(model_shannon_bracken_provincia))

cat("\nAnalysis of deviance (sequential, F-test - Gamma estimates its own dispersion):\n")
print(anova(model_shannon_bracken_provincia, test = "F"))

# Predicted-mean plot - same recipe as sections 5/7/9/11's interaction
# plots, but without Year (this model doesn't include it, so there's
# no second group to color/dodge by - just one point + 95% CI per
# Provincia).
newdata_shannon_bracken <- data.frame(
    Provincia_grouped = sort(unique(metadata_diversity_bracken$Provincia_grouped))
)

pred_shannon_bracken <- predict(model_shannon_bracken_provincia, newdata = newdata_shannon_bracken, type = "link", se.fit = TRUE)
newdata_shannon_bracken$fit <- exp(pred_shannon_bracken$fit)
newdata_shannon_bracken$lower <- exp(pred_shannon_bracken$fit - 1.96 * pred_shannon_bracken$se.fit)
newdata_shannon_bracken$upper <- exp(pred_shannon_bracken$fit + 1.96 * pred_shannon_bracken$se.fit)

print(
    ggplot(newdata_shannon_bracken, aes(x = Provincia_grouped, y = fit)) +
        geom_jitter(
            data = metadata_diversity_bracken, aes(x = Provincia_grouped, y = Shannon_bracken),
            width = 0.1, height = 0, alpha = 0.3, inherit.aes = FALSE
        ) +
        geom_pointrange(aes(ymin = lower, ymax = upper), color = "steelblue", linewidth = 1, size = 0.8) +
        labs(
            title = "Predicted Bracken Shannon diversity by Provincia (Gamma family, 95% CI)",
            y = "Shannon diversity (Bracken genus)", x = "Provincia"
        ) +
        theme_minimal()
)

