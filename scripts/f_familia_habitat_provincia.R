# ============================================================
# Gráfico combinado de composición a nivel de Familia:
# Bloque 1 -> Hábitat (Urbano, Agricola, Semi-natural)
# Bloque 2 -> Provincia (Bizkaia, Gipuzkoa, Araba)
# ============================================================

# 1. Librerías ----------------------------------------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(ragg)

# 2. Leer datos -----------------------------------------------------------
gtdbtk_sin_rejects <- readRDS("visualize/gtdbtk_sin_rejects.rds")

# Salvaguarda por si los nombres no se corrigieron en el origen
gtdbtk_sin_rejects <- gtdbtk_sin_rejects %>%
  mutate(
    Provincia = recode(Provincia, "Guipuzcoa" = "Gipuzkoa", "Vizcaya" = "Bizkaia"),
    Habitat   = recode(Habitat,   "Natural" = "Semi-natural", "Agricola" = "Agrícola")
  )

# 3. Ranking de Familia por Hábitat ---------------------------------------
ranking_familia_habitat <- gtdbtk_sin_rejects %>%
  filter(!is.na(Habitat)) %>%
  mutate(Familia = if_else(is.na(Familia) | Familia == "", "Sin determinar", Familia)) %>%
  count(Habitat, Familia, name = "Cantidad") %>%
  rename(Grupo = Habitat) %>%
  mutate(Bloque = "Hábitat")

# 4. Ranking de Familia por Provincia --------------------------------------
ranking_familia_provincia <- gtdbtk_sin_rejects %>%
  filter(!is.na(Provincia)) %>%
  mutate(Familia = if_else(is.na(Familia) | Familia == "", "Sin determinar", Familia)) %>%
  count(Provincia, Familia, name = "Cantidad") %>%
  rename(Grupo = Provincia) %>%
  mutate(Bloque = "Provincia")

# 5. Unir ambos bloques -----------------------------------------------------
datos_familia <- bind_rows(ranking_familia_habitat, ranking_familia_provincia)

# Orden del eje X dentro de cada bloque
orden_grupo <- c("Urbano", "Agrícola", "Semi-natural", "Bizkaia", "Gipuzkoa", "Araba")
datos_familia$Grupo <- factor(datos_familia$Grupo, levels = orden_grupo)

# Orden de los bloques (paneles), de izquierda a derecha
datos_familia$Bloque <- factor(datos_familia$Bloque, levels = c("Hábitat", "Provincia"))

# Orden de las familias por abundancia total (Sin determinar al final)
orden_familias <- datos_familia %>%
  group_by(Familia) %>%
  summarise(total = sum(Cantidad), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Familia)
orden_familias <- c(setdiff(orden_familias, "Sin determinar"), "Sin determinar")
datos_familia$Familia <- factor(datos_familia$Familia, levels = orden_familias)

# 6. Gráfico ------------------------------------------------------------
p_familia <- ggplot(datos_familia, aes(x = Grupo, y = Cantidad, fill = Familia)) +
  geom_col(position = "fill", color = "white") +
  facet_wrap(~ Bloque, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Abundancia relativa", fill = "Familia") +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y     = element_text(size = 11, color = "black"),
    legend.text     = element_text(size = 10, color = "black"),
    legend.position = "right",
    strip.text      = element_text(face = "bold", size = 12),
    panel.spacing   = unit(1, "lines")
  )
ragg::agg_png("plots/composicion_familia_habitat_provincia.png", width = 10, height = 6, units = "in", res = 300)
print(p_familia)
dev.off()
print(p_familia)

ggsave(
  filename = "plots/composicion_familia_habitat_provincia.png",
  plot     = p_familia,
  device   = ragg::agg_png,
  width    = 10,
  height   = 6,
  units    = "in",
  res      = 300
)
#ggsave(filename = "plots/composicion_familia_habitat_provincia.png", plot = p_familia, width = 10, height = 6, dpi = 300)
#ggsave(filename = "plots/composicion_familia_habitat_provincia.pdf", plot = p_familia, width = 10, height = 6)
#ggsave(filename = "plots/composicion_familia_habitat_provincia.svg", plot = p_familia, width = 10, height = 6)