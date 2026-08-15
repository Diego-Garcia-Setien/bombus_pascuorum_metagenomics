## ============================================================
## Este script consta de dos partes 
##
## Las partes 1ª y 2ª (lectura del TSV, separación de taxonomía,
## Especie_final, y el cruce con localización/provincia/hábitat) son
## IDÉNTICAS a las del script de Hábitat -- no hace falta repetirlas.
##
## Tampoco hace falta repetir el filtrado de Rejects: el script de
## Hábitat ya genera "visualize/gtdbtk_sin_rejects.rds" (que también
## tiene la columna Provincia, porque viene del mismo cruce con
## tabla_localizaciones). Este script lee ese archivo directamente.
##
## Requiere haber corrido antes:
##   - a_checkm2_quality_classification.R  -> checkm2_clasificado.rds
##   - b_gtdbtk_habitat.R (1ª a 3ª parte)  -> gtdbtk_sin_rejects.rds
##
## *****1ª parte*****
##
## Se hace un análisis de taxonomía por PROVINCIA (Araba / Vizcaya /
## Guipuzcoa), sobre los MAGs que ya vienen sin Rejects.
##
## *****2ª parte*****
##
## Se hace el mismo análisis de la 3ª parte pero excluyendo los 5 taxones
## más abundantes (los mismos 5 en las 3 provincias), para poder ver mejor
## las diferencias entre la microbiota "secundaria", no core.
## ============================================================

# 1. Cargar librerías -------------------------------------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(scales)
library(patchwork)

# 2. Leer el archivo intermedio ya filtrado (sin Rejects) --------------------
gtdbtk_sin_rejects <- readRDS("visualize/gtdbtk_sin_rejects.rds")

# Chequeo rápido: confirmar que la columna Provincia está presente y sin NA
table(gtdbtk_sin_rejects$Provincia, useNA = "always")

# 3. Funciones reutilizables para los gráficos de barras apiladas -----------
preparar_datos_stack <- function(ranking_tabla, columna_taxon, columna_grupo) {
  ranking_tabla %>%
    rename(TaxonPlot = all_of(columna_taxon), GrupoPlot = all_of(columna_grupo)) %>%
    group_by(GrupoPlot) %>%
    mutate(Porcentaje = Cantidad / sum(Cantidad) * 100) %>%
    ungroup()
}

graficar_stack_grupo <- function(datos, titulo) {
  
  orden_taxones <- datos %>%
    group_by(TaxonPlot) %>%
    summarise(total = sum(Cantidad)) %>%
    arrange(desc(total)) %>%
    pull(TaxonPlot)
  orden_taxones <- c(setdiff(orden_taxones, "Sin determinar"), "Sin determinar")
  datos$TaxonPlot <- factor(datos$TaxonPlot, levels = orden_taxones)
  
  p_abs <- ggplot(datos, aes(x = GrupoPlot, y = Cantidad, fill = TaxonPlot)) +
    geom_col(position = "stack", color = "white") +
    coord_flip() +
    labs(x = NULL, y = "Cantidad de MAGs", fill = NULL) +
    theme_minimal()
  
  p_pct <- ggplot(datos, aes(x = GrupoPlot, y = Cantidad, fill = TaxonPlot)) +
    geom_col(position = "fill", color = "white") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Porcentaje", fill = NULL) +
    theme_minimal()
  
  (p_abs + p_pct) +
    plot_layout(guides = "collect") +
    plot_annotation(title = titulo,
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
}

##=======================================================
#########################################################
##=======================================================
######################################################||#
##                                                    ||#
## 1ª parte del script                                ||#
##                                                    ||#
######################################################||#
##=======================================================
#########################################################
##=======================================================

# 4. Ranking de Género por Provincia, sin Rejects (cantidad y porcentaje) --
ranking_genero_provincia_without_rejects <- gtdbtk_sin_rejects %>%
  filter(!is.na(Provincia)) %>%
  count(Provincia, Genero, name = "Cantidad") %>%
  group_by(Provincia) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Provincia, desc(Cantidad)) %>%
  ungroup()

ranking_genero_provincia_without_rejects

# 5. Ranking de Especie por Provincia, sin Rejects (cantidad y porcentaje) -
ranking_especie_provincia_without_rejects <- gtdbtk_sin_rejects %>%
  filter(!is.na(Provincia)) %>%
  count(Provincia, Especie_final, name = "Cantidad") %>%
  group_by(Provincia) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Provincia, desc(Cantidad)) %>%
  ungroup()

ranking_especie_provincia_without_rejects

# 6. Exportar taxonomía filtrada + rankings a Excel --------------------------
wb <- createWorkbook()

addWorksheet(wb, "Taxonomy_without_rejects")
writeData(wb, "Taxonomy_without_rejects", gtdbtk_sin_rejects)

addWorksheet(wb, "Ranking_Genus_Provincia")
writeData(wb, "Ranking_Genus_Provincia", ranking_genero_provincia_without_rejects)

addWorksheet(wb, "Ranking_Species_Provincia")
writeData(wb, "Ranking_Species_Provincia", ranking_especie_provincia_without_rejects)

#saveWorkbook(wb, "excels/taxonomy_MAGs_gtdbtk_provincia_without_rejects.xlsx", overwrite = TRUE)

#write_tsv(gtdbtk_sin_rejects, "visualize/taxonomy_without_rejects_provincia.tsv")
#write_tsv(ranking_genero_provincia_without_rejects, "visualize/ranking_genus_provincia_without_rejects.tsv")
#write_tsv(ranking_especie_provincia_without_rejects, "visualize/ranking_species_provincia_without_rejects.tsv")

# 7. Gráficos de barras apiladas horizontales, por provincia, sin Rejects --
datos_stack_genero_sr <- preparar_datos_stack(ranking_genero_provincia_without_rejects, "Genero", "Provincia")
p_genero_sr <- graficar_stack_grupo(datos_stack_genero_sr, "Composición de Géneros por Provincia (sin Rejects)")
print(p_genero_sr)
#ggsave(filename = "plots/stacked_bars_genus_provincia_without_rejects.png", plot = p_genero_sr, width = 12, height = 5, dpi = 300)
#ggsave(filename = "plots/stacked_bars_genus_provincia_without_rejects.pdf", plot = p_genero_sr, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_genus_provincia_without_rejects.svg", plot = p_genero_sr, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_genus_provincia_without_rejects.eps", plot = p_genero_sr, width = 12, height = 5, device = cairo_ps)

datos_stack_especie_sr <- preparar_datos_stack(ranking_especie_provincia_without_rejects, "Especie_final", "Provincia")
p_especie_sr <- graficar_stack_grupo(datos_stack_especie_sr, "Composición de Especies por Provincia (sin Rejects)")
print(p_especie_sr)
#ggsave(filename = "plots/stacked_bars_species_provincia_without_rejects.png", plot = p_especie_sr, width = 12, height = 5, dpi = 300)
#ggsave(filename = "plots/stacked_bars_species_provincia_without_rejects.pdf", plot = p_especie_sr, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_species_provincia_without_rejects.svg", plot = p_especie_sr, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_species_provincia_without_rejects.eps", plot = p_especie_sr, width = 12, height = 5, device = cairo_ps)

# 8. Excels con cantidad y porcentaje de Género/Especie por Provincia, sin Rejects
# --- Género -----------------------------------------------------------------
cantidad_genero <- ranking_genero_provincia_without_rejects %>%
  select(Provincia, Genero, Cantidad) %>%
  pivot_wider(names_from = Provincia, values_from = Cantidad, values_fill = 0) %>%
  arrange(Genero)

porcentaje_genero <- ranking_genero_provincia_without_rejects %>%
  select(Provincia, Genero, Porcentaje) %>%
  pivot_wider(names_from = Provincia, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Genero)

wb_genero <- createWorkbook()
addWorksheet(wb_genero, "Cantidad")
writeData(wb_genero, "Cantidad", cantidad_genero)
addWorksheet(wb_genero, "Porcentaje")
writeData(wb_genero, "Porcentaje", porcentaje_genero)
#saveWorkbook(wb_genero, "excels/resumen_genero_por_provincia_sin_rejects.xlsx", overwrite = TRUE)

# --- Especie -----------------------------------------------------------------
cantidad_especie <- ranking_especie_provincia_without_rejects %>%
  select(Provincia, Especie_final, Cantidad) %>%
  pivot_wider(names_from = Provincia, values_from = Cantidad, values_fill = 0) %>%
  arrange(Especie_final)

porcentaje_especie <- ranking_especie_provincia_without_rejects %>%
  select(Provincia, Especie_final, Porcentaje) %>%
  pivot_wider(names_from = Provincia, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Especie_final)

wb_especie <- createWorkbook()
addWorksheet(wb_especie, "Cantidad")
writeData(wb_especie, "Cantidad", cantidad_especie)
addWorksheet(wb_especie, "Porcentaje")
writeData(wb_especie, "Porcentaje", porcentaje_especie)
#saveWorkbook(wb_especie, "excels/resumen_especie_por_provincia_sin_rejects.xlsx", overwrite = TRUE)

# 9. Guardar archivo intermedio (para usar en otros scripts, ej. Familia) --
#saveRDS(gtdbtk_sin_rejects, "visualize/gtdbtk_sin_rejects_provincia.rds")

##=======================================================
#########################################################
##=======================================================
######################################################||#
##                                                    ||#
## 2ª parte del script                                ||#
##                                                    ||#
######################################################||#
##=======================================================
#########################################################
##=======================================================

# 10. Identificar las 5 especies más abundantes en total (3 provincias juntas)
top5_especies <- ranking_especie_provincia_without_rejects %>%
  group_by(Especie_final) %>%
  summarise(Total = sum(Cantidad), .groups = "drop") %>%
  slice_max(Total, n = 5) %>%
  pull(Especie_final)

top5_especies

# 11. Filtrar el dataset, sacando esas 5 especies -----------------------------
gtdbtk_sin_top5 <- gtdbtk_sin_rejects %>%
  filter(!Especie_final %in% top5_especies)

nrow(gtdbtk_sin_rejects)
nrow(gtdbtk_sin_top5)

# 12. Ranking de Especie por Provincia, sin Rejects y sin el top 5 ----------
ranking_especie_provincia_sin_top5 <- gtdbtk_sin_top5 %>%
  filter(!is.na(Provincia)) %>%
  count(Provincia, Especie_final, name = "Cantidad") %>%
  group_by(Provincia) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Provincia, desc(Cantidad)) %>%
  ungroup()

ranking_especie_provincia_sin_top5

# 13. Gráfico de barras apiladas horizontales, sin el top 5 de especies -----
datos_stack_especie_sin_top5 <- preparar_datos_stack(ranking_especie_provincia_sin_top5, "Especie_final", "Provincia")
p_especie_sin_top5 <- graficar_stack_grupo(
  datos_stack_especie_sin_top5,
  "Composición de Especies por Provincia (sin Rejects, sin top 5 especies)"
)
print(p_especie_sin_top5)

#ggsave(filename = "plots/stacked_bars_species_provincia_sin_top5.png", plot = p_especie_sin_top5, width = 12, height = 5, dpi = 300)
#ggsave(filename = "plots/stacked_bars_species_provincia_sin_top5.pdf", plot = p_especie_sin_top5, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_species_provincia_sin_top5.svg", plot = p_especie_sin_top5, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_species_provincia_sin_top5.eps", plot = p_especie_sin_top5, width = 12, height = 5, device = cairo_ps)

# 14. Excel con cantidad y porcentaje de Especie por Provincia, sin top 5 --
cantidad_especie_sin_top5 <- ranking_especie_provincia_sin_top5 %>%
  select(Provincia, Especie_final, Cantidad) %>%
  pivot_wider(names_from = Provincia, values_from = Cantidad, values_fill = 0) %>%
  arrange(Especie_final)

porcentaje_especie_sin_top5 <- ranking_especie_provincia_sin_top5 %>%
  select(Provincia, Especie_final, Porcentaje) %>%
  pivot_wider(names_from = Provincia, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Especie_final)

wb_especie_sin_top5 <- createWorkbook()
addWorksheet(wb_especie_sin_top5, "Cantidad")
writeData(wb_especie_sin_top5, "Cantidad", cantidad_especie_sin_top5)
addWorksheet(wb_especie_sin_top5, "Porcentaje")
writeData(wb_especie_sin_top5, "Porcentaje", porcentaje_especie_sin_top5)
#saveWorkbook(wb_especie_sin_top5, "excels/resumen_especie_por_provincia_sin_top5.xlsx", overwrite = TRUE)

# 15. TSV con el ranking completo, sin top 5 especies -------------------------
#write_tsv(ranking_especie_provincia_sin_top5, "visualize/ranking_species_provincia_sin_top5.tsv")

# 16. Identificar los 5 géneros más abundantes en total ------------------------
top5_generos <- ranking_genero_provincia_without_rejects %>%
  group_by(Genero) %>%
  summarise(Total = sum(Cantidad), .groups = "drop") %>%
  slice_max(Total, n = 5) %>%
  pull(Genero)

top5_generos

# 17. Filtrar el dataset, sacando esos 5 géneros ------------------------------
gtdbtk_sin_top5_genero <- gtdbtk_sin_rejects %>%
  filter(!Genero %in% top5_generos)

nrow(gtdbtk_sin_rejects)
nrow(gtdbtk_sin_top5_genero)

# 18. Ranking de Género por Provincia, sin Rejects y sin el top 5 -----------
ranking_genero_provincia_sin_top5 <- gtdbtk_sin_top5_genero %>%
  filter(!is.na(Provincia)) %>%
  count(Provincia, Genero, name = "Cantidad") %>%
  group_by(Provincia) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Provincia, desc(Cantidad)) %>%
  ungroup()

ranking_genero_provincia_sin_top5

# 19. Gráfico de barras apiladas horizontales, sin el top 5 de géneros ------
datos_stack_genero_sin_top5 <- preparar_datos_stack(ranking_genero_provincia_sin_top5, "Genero", "Provincia")
p_genero_sin_top5 <- graficar_stack_grupo(
  datos_stack_genero_sin_top5,
  "Composición de Géneros por Provincia (sin Rejects, sin top 5 géneros)"
)
print(p_genero_sin_top5)

#ggsave(filename = "plots/stacked_bars_genus_provincia_sin_top5.png", plot = p_genero_sin_top5, width = 12, height = 5, dpi = 300)
#ggsave(filename = "plots/stacked_bars_genus_provincia_sin_top5.pdf", plot = p_genero_sin_top5, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_genus_provincia_sin_top5.svg", plot = p_genero_sin_top5, width = 12, height = 5)
#ggsave(filename = "plots/stacked_bars_genus_provincia_sin_top5.eps", plot = p_genero_sin_top5, width = 12, height = 5, device = cairo_ps)

# 20. Excel con cantidad y porcentaje de Género por Provincia, sin top 5 ----
cantidad_genero_sin_top5 <- ranking_genero_provincia_sin_top5 %>%
  select(Provincia, Genero, Cantidad) %>%
  pivot_wider(names_from = Provincia, values_from = Cantidad, values_fill = 0) %>%
  arrange(Genero)

porcentaje_genero_sin_top5 <- ranking_genero_provincia_sin_top5 %>%
  select(Provincia, Genero, Porcentaje) %>%
  pivot_wider(names_from = Provincia, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Genero)

wb_genero_sin_top5 <- createWorkbook()
addWorksheet(wb_genero_sin_top5, "Cantidad")
writeData(wb_genero_sin_top5, "Cantidad", cantidad_genero_sin_top5)
addWorksheet(wb_genero_sin_top5, "Porcentaje")
writeData(wb_genero_sin_top5, "Porcentaje", porcentaje_genero_sin_top5)
#saveWorkbook(wb_genero_sin_top5, "excels/resumen_genero_por_provincia_sin_top5.xlsx", overwrite = TRUE)

# 21. TSV con el ranking completo, sin top 5 géneros ---------------------------
#write_tsv(ranking_genero_provincia_sin_top5, "visualize/ranking_genus_provincia_sin_top5.tsv")

#save.image("visualize/b_gtdbtk_provincia.RData")