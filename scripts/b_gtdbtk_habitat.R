## ============================================================
## Este script consta de cuatro partes
##
## *****1ª parte*****
##
## Procesamiento de clasificación taxonómica (GTDB-Tk)
##
## A partir de un TSV con columnas "user_genome" y "classification":
##   - separa la clasificación en Dominio/Phylum/Clase/Orden/Familia/Genero/Especie
##   - extrae metadatos de muestra (organismo, localización, año, etc.)
##     a partir de "user_genome", igual que se hizo con el reporte de CheckM2
##
## Formato de "user_genome" asumido:
##   BPART241001__BPART241001.MAGScoT_cleanbin_000001
##   |__________|  |________________________________|
##   prefijo dup.   Name real (mismo formato que en CheckM2)
##
## Formato de "classification" asumido:
##   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;
##   f__Orbaceae;g__Gilliamella;s__Gilliamella mensalis
##   (la especie puede venir vacía: "s__" -> se convierte a NA)
##
## *****2ª parte*****
##
## Se agrega una nueva columna de especie final, para el caso de aquellos
## MAGs donde más de una vez solo se ha identificado a nivel de género a 
## un mismo género, se le añade el término "indeterminado" en la columna de
## especie final, como el caso de Xylocopillactobacillus. Para los MAGs que 
## solo han sido identificados a nivel de género pero el género es único,
## solamente se repite una vez, se le añade la abreviatura de especie "sp"
## ya que se puede asegurar que es una unica especie, en el caso de 
## Xylocopillactobacillus, pueden ser tantas especies como MAGs identificados
## con ese género
##
## *****3ª parte*****
##
## Se hace un análisis de taxonomía por hábitat, EXCLUYENDO los MAGs
## clasificados como "Rejects" en el script de calidad (CheckM2)
##
## A partir de los valores de completitud y contaminación obtenidos con CheckM2,
## los MAGs se clasificaron en cinco categorías de calidad,
## basadas en los criterios MIMAG (Bowers et al. 2017)
## con dos categorías adicionales ("Calidad aceptable" y "Rejected") 
##
## Este script repite exactamente el mismo análisis,
## pero sobre el subconjunto de MAGs con calidad
## aceptable o mejor (se descartan los "Rejects").
##
## *****4ª parte*****
##
## Se hace el mismo análisis de la 3ª parte pero excluyendo los 5 taxones
## más abundantes los cuales son principalmente bacterias más típicas de 
## las microbiota intestinal
##
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
##
##
# 1. Cargar librerías -------------------------------------------------------
# *instalar los paquetes si no están previamente instalados
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(scales)
library(patchwork)

# 2. Leer el archivo TSV de GTDB-Tk ------------------------------------------
gtdbtk <- read_tsv("visualize/gtdbtk.bac120.summary.tsv")

colnames(gtdbtk)

# 3. Extraer el "Name" real desde user_genome --------------------------------
# Quita el prefijo duplicado antes del "__", dejando el mismo formato
# que la columna "Name" del reporte de CheckM2
gtdbtk <- gtdbtk %>%
  mutate(Name = str_remove(user_genome, "^[^_]+__"))

# 4. Separar la clasificación taxonómica en columnas -------------------------
gtdbtk <- gtdbtk %>%
  separate(
    classification,
    into  = c("Dominio", "Phylum", "Clase", "Orden", "Familia", "Genero", "Especie"),
    sep   = ";",
    fill  = "right",   # por si algún registro viniera con menos niveles
    remove = FALSE      # conserva la columna original "classification" también
  ) %>%
  mutate(across(Dominio:Especie, ~ str_remove(., "^[a-z]__"))) %>%  # quita "d__", "p__", etc.
  mutate(across(Dominio:Especie, ~ na_if(., "")))                    # "" -> NA (ej. especie sin determinar)

# 5. Extraer metadatos de muestra desde Name ---------------------------------
regex_mag <- "^([A-Z])P([A-Z]{3})(\\d{2})(\\d{4})([A-Z]?)\\.(.+)$"

matched <- str_match(gtdbtk$Name, regex_mag)
dim(matched)  # debería tener 7 columnas (1 match completo + 6 grupos)

gtdbtk <- gtdbtk %>%
  mutate(
    Organismo     = matched[, 2],
    Localizacion  = matched[, 3],
    Anio          = as.numeric(matched[, 4]) + 2000,
    CodigoMuestra = matched[, 5],
    Sufijo        = matched[, 6],
    BinID         = matched[, 7]
  )

# Nombres que no matchearon el patrón esperado (debería dar 0)
sum(is.na(gtdbtk$Localizacion))

# 6. Resumen rápido ----------------------------------------------------
table(gtdbtk$Dominio)
table(gtdbtk$Phylum)

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

# 7. Agregar columna "Especie_final" con la especie refinada -----------------
# Reglas:
#  - Si la especie SÍ está determinada -> se deja igual.
#  - Si NO está determinada y ese género tiene una única muestra sin
#    especie en todo el dataset -> se puede asumir con seguridad "Genero sp."
#  - Si tiene VARIAS muestras sin especie determinada -> no podemos asumir
#    que sean la misma especie -> "Genero indeterminado"
conteo_sin_especie <- gtdbtk %>%
  filter(is.na(Especie)) %>%
  count(Genero, name = "n_sin_especie")

gtdbtk <- gtdbtk %>%
  left_join(conteo_sin_especie, by = "Genero") %>%
  mutate(
    Especie_final = case_when(
      !is.na(Especie)                     ~ Especie,
      is.na(Especie) & n_sin_especie == 1 ~ paste(Genero, "sp."),
      is.na(Especie) & n_sin_especie > 1  ~ paste(Genero, "indeterminado"),
      TRUE                                 ~ NA_character_
    )
  ) %>%
  select(-n_sin_especie) %>%
  relocate(Especie_final, .after = Especie)

# Ver qué quedó en cada caso sin especie determinada
gtdbtk %>%
  filter(is.na(Especie)) %>%
  distinct(Genero, Especie_final)

# 8. Guardar archivo intermedio (versión SIN cruce de hábitat todavía) ------
saveRDS(gtdbtk, "visualize/gtdbtk_taxonomia.rds")

# 9. Cruzar con la tabla de localización / provincia / hábitat --------------
# (ÚNICO join en todo el script -- hacerlo dos veces duplica las columnas
# Locality/Provincia/Habitat en .x/.y y rompe todo lo que sigue)
tabla_localizaciones <- tribble(
  ~Localizacion, ~Locality,          ~Provincia,   ~Habitat,
  "VAL",         "Valderejo",        "Araba",      "Natural",
  "KAR",         "Karkamo",          "Araba",      "Agricola",
  "GAS",         "Vitoria_Gasteiz",  "Araba",      "Urbano",
  "GOR",         "Gorbea",           "Vizcaya",    "Natural",
  "ZEB",         "Zeberio",          "Vizcaya",    "Agricola",
  "BIL",         "Bilbao",           "Vizcaya",    "Urbano",
  "ART",         "Artikutza",        "Guipuzcoa",  "Natural",
  "ATZ",         "Asteazu",          "Guipuzcoa",  "Agricola",
  "DON",         "Donostia",         "Guipuzcoa",  "Urbano"
)

gtdbtk <- gtdbtk %>%
  left_join(tabla_localizaciones, by = "Localizacion")

# Chequeo: localizaciones que no matchearon con la tabla de referencia (debería dar 0 filas)
gtdbtk %>% filter(is.na(Habitat)) %>% distinct(Localizacion)

# 10. Ranking de Género por Hábitat (cantidad y porcentaje) -------------------
ranking_genero_habitat <- gtdbtk %>%
  filter(!is.na(Habitat)) %>%
  mutate(Genero = if_else(is.na(Genero) | Genero == "", "Sin determinar", Genero)) %>%
  count(Habitat, Genero, name = "Cantidad") %>%
  group_by(Habitat) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Habitat, desc(Cantidad)) %>%
  ungroup()

ranking_genero_habitat

# 11. Ranking de Especie por Hábitat (cantidad y porcentaje) ------------------
# Usa "Especie_final" (con las etiquetas "sp." / "indeterminado" ya aplicadas)
ranking_especie_habitat <- gtdbtk %>%
  filter(!is.na(Habitat)) %>%
  count(Habitat, Especie_final, name = "Cantidad") %>%
  group_by(Habitat) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Habitat, desc(Cantidad)) %>%
  ungroup()

ranking_especie_habitat

# 12. Guardar archivo intermedio (para la 3ª parte: filtrado de Rejects) ----
saveRDS(gtdbtk, "visualize/gtdbtk_taxonomia.rds")

# 13. Exportar a Excel ----------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "Taxonomia_completa")
writeData(wb, "Taxonomia_completa", gtdbtk)

addWorksheet(wb, "Ranking_Genero_Habitat")
writeData(wb, "Ranking_Genero_Habitat", ranking_genero_habitat)

addWorksheet(wb, "Ranking_Especie_Habitat")
writeData(wb, "Ranking_Especie_Habitat", ranking_especie_habitat)

#saveWorkbook(wb, "excels/taxonomia_MAGs_gtdbtk.xlsx", overwrite = TRUE)

# 14. TSV de cada tabla por separado ------------------------------------------
#write_tsv(gtdbtk, "visualize/taxonomy_all_species.tsv")
#write_tsv(ranking_genero_habitat, "visualize/ranking_genus_habitat_all_species.tsv")
#write_tsv(ranking_especie_habitat, "visualize/ranking_species_habitat_all_species.tsv")

# 15. Gráficos de barras apiladas horizontales, por hábitat -------------------
# Uno para Género y otro para Especie. Cada uno combina dos paneles:
# cantidad absoluta (izquierda) y porcentaje (derecha). Se muestran TODOS
# los taxones (sin agrupar en "Otros").

preparar_datos_stack <- function(ranking_tabla, columna_taxon) {
  ranking_tabla %>%
    rename(TaxonPlot = all_of(columna_taxon)) %>%
    group_by(Habitat) %>%
    mutate(Porcentaje = Cantidad / sum(Cantidad) * 100) %>%
    ungroup()
}

graficar_stack_habitat <- function(datos, titulo) {
  
  orden_taxones <- datos %>%
    group_by(TaxonPlot) %>%
    summarise(total = sum(Cantidad)) %>%
    arrange(desc(total)) %>%
    pull(TaxonPlot)
  orden_taxones <- c(setdiff(orden_taxones, "Sin determinar"), "Sin determinar")
  datos$TaxonPlot <- factor(datos$TaxonPlot, levels = orden_taxones)
  
  p_abs <- ggplot(datos, aes(x = Habitat, y = Cantidad, fill = TaxonPlot)) +
    geom_col(position = "stack", color = "white") +
    coord_flip() +
    labs(x = NULL, y = "Cantidad de MAGs", fill = NULL) +
    theme_minimal()
  
  p_pct <- ggplot(datos, aes(x = Habitat, y = Cantidad, fill = TaxonPlot)) +
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

# Género --------------------------------------------------------------------
datos_stack_genero <- preparar_datos_stack(ranking_genero_habitat, "Genero")
p_genero <- graficar_stack_habitat(datos_stack_genero, "Composición de Géneros por Hábitat")
print(p_genero)
ggsave(filename = "plots/barras_apiladas_genero_habitat.png", plot = p_genero, width = 12, height = 5, dpi = 300)
ggsave(filename = "plots/barras_apiladas_genero_habitat.pdf", plot = p_genero, width = 12, height = 5)
ggsave(filename = "plots/barras_apiladas_genero_habitat.svg", plot = p_genero, width = 12, height = 5)
ggsave(filename = "plots/barras_apiladas_genero_habitat.eps", plot = p_genero, width = 12, height = 5, device = cairo_ps)

# Especie ---------------------------------------------------------------------
datos_stack_especie <- preparar_datos_stack(ranking_especie_habitat, "Especie_final")
p_especie <- graficar_stack_habitat(datos_stack_especie, "Composición de Especies por Hábitat")
print(p_especie)
ggsave(filename = "plots/barras_apiladas_especie_habitat.png", plot = p_especie, width = 12, height = 5, dpi = 300)
ggsave(filename = "plots/barras_apiladas_especie_habitat.pdf", plot = p_especie, width = 12, height = 5)
ggsave(filename = "plots/barras_apiladas_especie_habitat.svg", plot = p_especie, width = 12, height = 5)
ggsave(filename = "plots/barras_apiladas_especie_habitat.eps", plot = p_especie, width = 12, height = 5, device = cairo_ps)

##=======================================================
#########################################################
##=======================================================
######################################################||#
##                                                    ||#
## 3ª parte del script                                ||#
##                                                    ||#
######################################################||#
##=======================================================
#########################################################
##=======================================================

# 16. Leer los archivos intermedios necesarios -------------------------------
checkm2_clasificado <- readRDS("visualize/checkm2_clasificado.rds")
#gtdbtk              <- readRDS("visualize/gtdbtk_taxonomia.rds")

# 17. Unir calidad + taxonomía por "Name", y descartar los Rejects -----------
gtdbtk_sin_rejects <- gtdbtk %>%
  left_join(checkm2_clasificado %>% select(Name, Categoria), by = "Name")

# Chequeo: MAGs de la taxonomía sin pareja en checkm2 (deberían ser 0)
sum(is.na(gtdbtk_sin_rejects$Categoria))

# Cuántos MAGs había en total, y cuántos quedan tras descartar Rejects
nrow(gtdbtk_sin_rejects)                     # total original
table(gtdbtk_sin_rejects$Categoria)          # por categoría de calidad

gtdbtk_sin_rejects <- gtdbtk_sin_rejects %>%
  filter(Categoria != "Rejects")

nrow(gtdbtk_sin_rejects)                     # total tras descartar Rejects

# 18. Ranking de Género por Hábitat, sin Rejects (cantidad y porcentaje) ----
ranking_genero_habitat_without_rejects <- gtdbtk_sin_rejects %>%
  filter(!is.na(Habitat)) %>%
  count(Habitat, Genero, name = "Cantidad") %>%
  group_by(Habitat) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Habitat, desc(Cantidad)) %>%
  ungroup()

ranking_genero_habitat_without_rejects

# 19. Ranking de Especie por Hábitat, sin Rejects (cantidad y porcentaje) ---
ranking_especie_habitat_without_rejects <- gtdbtk_sin_rejects %>%
  filter(!is.na(Habitat)) %>%
  count(Habitat, Especie_final, name = "Cantidad") %>%
  group_by(Habitat) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Habitat, desc(Cantidad)) %>%
  ungroup()

ranking_especie_habitat_without_rejects

# 20. Exportar taxonomía filtrada + rankings a Excel --------------------------
wb <- createWorkbook()

addWorksheet(wb, "Taxonomy_without_rejects")
writeData(wb, "Taxonomy_without_rejects", gtdbtk_sin_rejects)

addWorksheet(wb, "Ranking_Genus_Habitat")
writeData(wb, "Ranking_Genus_Habitat", ranking_genero_habitat_without_rejects)

addWorksheet(wb, "Ranking_Species_Habitat")
writeData(wb, "Ranking_Species_Habitat", ranking_especie_habitat_without_rejects)

saveWorkbook(wb, "excels/taxonomy_MAGs_gtdbtk_without_rejects.xlsx", overwrite = TRUE)

# TSV de cada tabla por separado
write_tsv(gtdbtk_sin_rejects, "visualize/taxonomy_without_rejects.tsv")
write_tsv(ranking_genero_habitat_without_rejects, "visualize/ranking_genus_habitat_without_rejects.tsv")
write_tsv(ranking_especie_habitat_without_rejects, "visualize/ranking_species_habitat_without_rejects.tsv")

# 21. Gráficos de barras apiladas horizontales, por hábitat, sin Rejects ----
# Reutiliza las mismas funciones preparar_datos_stack / graficar_stack_habitat
# definidas en la 2ª parte -- no hace falta redefinirlas.

datos_stack_genero_sr <- preparar_datos_stack(ranking_genero_habitat_without_rejects, "Genero")
p_genero_sr <- graficar_stack_habitat(datos_stack_genero_sr, "Composición de Géneros por Hábitat (sin Rejects)")
print(p_genero_sr)
ggsave(filename = "plots/stacked_bars_genus_habitat_without_rejects.png", plot = p_genero_sr, width = 12, height = 5, dpi = 300)
ggsave(filename = "plots/stacked_bars_genus_habitat_without_rejects.pdf", plot = p_genero_sr, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_genus_habitat_without_rejects.svg", plot = p_genero_sr, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_genus_habitat_without_rejects.eps", plot = p_genero_sr, width = 12, height = 5, device = cairo_ps)

datos_stack_especie_sr <- preparar_datos_stack(ranking_especie_habitat_without_rejects, "Especie_final")
p_especie_sr <- graficar_stack_habitat(datos_stack_especie_sr, "Composición de Especies por Hábitat (sin Rejects)")
print(p_especie_sr)
ggsave(filename = "plots/stacked_bars_species_habitat_without_rejects.png", plot = p_especie_sr, width = 12, height = 5, dpi = 300)
ggsave(filename = "plots/stacked_bars_species_habitat_without_rejects.pdf", plot = p_especie_sr, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_species_habitat_without_rejects.svg", plot = p_especie_sr, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_species_habitat_without_rejects.eps", plot = p_especie_sr, width = 12, height = 5, device = cairo_ps)

# 22. Excels con cantidad y porcentaje de Género/Especie por Hábitat, sin Rejects
# --- Género -----------------------------------------------------------------
cantidad_genero <- ranking_genero_habitat_without_rejects %>%
  select(Habitat, Genero, Cantidad) %>%
  pivot_wider(names_from = Habitat, values_from = Cantidad, values_fill = 0) %>%
  arrange(Genero)

porcentaje_genero <- ranking_genero_habitat_without_rejects %>%
  select(Habitat, Genero, Porcentaje) %>%
  pivot_wider(names_from = Habitat, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Genero)

wb_genero <- createWorkbook()
addWorksheet(wb_genero, "Cantidad")
writeData(wb_genero, "Cantidad", cantidad_genero)
addWorksheet(wb_genero, "Porcentaje")
writeData(wb_genero, "Porcentaje", porcentaje_genero)
#saveWorkbook(wb_genero, "excels/resumen_genero_por_habitat_sin_rejects.xlsx", overwrite = TRUE)

# --- Especie -----------------------------------------------------------------
cantidad_especie <- ranking_especie_habitat_without_rejects %>%
  select(Habitat, Especie_final, Cantidad) %>%
  pivot_wider(names_from = Habitat, values_from = Cantidad, values_fill = 0) %>%
  arrange(Especie_final)

porcentaje_especie <- ranking_especie_habitat_without_rejects %>%
  select(Habitat, Especie_final, Porcentaje) %>%
  pivot_wider(names_from = Habitat, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Especie_final)

wb_especie <- createWorkbook()
addWorksheet(wb_especie, "Cantidad")
writeData(wb_especie, "Cantidad", cantidad_especie)
addWorksheet(wb_especie, "Porcentaje")
writeData(wb_especie, "Porcentaje", porcentaje_especie)
saveWorkbook(wb_especie, "excels/resumen_especie_por_habitat_sin_rejects.xlsx", overwrite = TRUE)

# 23. Guardar archivo intermedio (para usar en otros scripts, ej. Familia) --
saveRDS(gtdbtk_sin_rejects, "visualize/gtdbtk_sin_rejects.rds")


##=======================================================
#########################################################
##=======================================================
######################################################||#
##                                                    ||#
## 4ª parte del script                                ||#
##                                                    ||#
######################################################||#
##=======================================================
#########################################################
##=======================================================

# 24. Identificar las 5 especies más abundantes en total ---------------------
top5_especies <- ranking_especie_habitat_without_rejects %>%
  group_by(Especie_final) %>%
  summarise(Total = sum(Cantidad), .groups = "drop") %>%
  slice_max(Total, n = 5) %>%
  pull(Especie_final)

top5_especies  # revisar cuáles son antes de excluirlas

# 25. Filtrar el dataset, sacando esas 5 especies -----------------------------
gtdbtk_sin_top5 <- gtdbtk_sin_rejects %>%
  filter(!Especie_final %in% top5_especies)

nrow(gtdbtk_sin_rejects)   # total antes de sacar las 5 especies
nrow(gtdbtk_sin_top5)      # total después

# 26. Ranking de Especie por Hábitat, sin Rejects y sin el top 5 ------------
ranking_especie_habitat_sin_top5 <- gtdbtk_sin_top5 %>%
  filter(!is.na(Habitat)) %>%
  count(Habitat, Especie_final, name = "Cantidad") %>%
  group_by(Habitat) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Habitat, desc(Cantidad)) %>%
  ungroup()

ranking_especie_habitat_sin_top5

# 27. Gráfico de barras apiladas horizontales, sin el top 5 ------------------
# Reutiliza preparar_datos_stack() / graficar_stack_habitat(), ya definidas
datos_stack_especie_sin_top5 <- preparar_datos_stack(ranking_especie_habitat_sin_top5, "Especie_final")
p_especie_sin_top5 <- graficar_stack_habitat(
  datos_stack_especie_sin_top5,
  "Composición de Especies por Hábitat (sin Rejects, sin top 5 especies)"
)
print(p_especie_sin_top5)

ggsave(filename = "plots/stacked_bars_species_habitat_sin_top5.png", plot = p_especie_sin_top5, width = 12, height = 5, dpi = 300)
ggsave(filename = "plots/stacked_bars_species_habitat_sin_top5.pdf", plot = p_especie_sin_top5, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_species_habitat_sin_top5.svg", plot = p_especie_sin_top5, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_species_habitat_sin_top5.eps", plot = p_especie_sin_top5, width = 12, height = 5, device = cairo_ps)

# 28. Excel con cantidad y porcentaje de Especie por Hábitat, sin top 5 -----
cantidad_especie_sin_top5 <- ranking_especie_habitat_sin_top5 %>%
  select(Habitat, Especie_final, Cantidad) %>%
  pivot_wider(names_from = Habitat, values_from = Cantidad, values_fill = 0) %>%
  arrange(Especie_final)

porcentaje_especie_sin_top5 <- ranking_especie_habitat_sin_top5 %>%
  select(Habitat, Especie_final, Porcentaje) %>%
  pivot_wider(names_from = Habitat, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Especie_final)

wb_especie_sin_top5 <- createWorkbook()
addWorksheet(wb_especie_sin_top5, "Cantidad")
writeData(wb_especie_sin_top5, "Cantidad", cantidad_especie_sin_top5)
addWorksheet(wb_especie_sin_top5, "Porcentaje")
writeData(wb_especie_sin_top5, "Porcentaje", porcentaje_especie_sin_top5)
saveWorkbook(wb_especie_sin_top5, "excels/resumen_especie_por_habitat_sin_top5.xlsx", overwrite = TRUE)

# 29. TSV con el ranking completo, sin top 5 especies -------------------------
write_tsv(ranking_especie_habitat_sin_top5, "visualize/ranking_species_habitat_sin_top5.tsv")

# 30. Identificar los 5 géneros más abundantes en total -----------------------
top5_generos <- ranking_genero_habitat_without_rejects %>%
  group_by(Genero) %>%
  summarise(Total = sum(Cantidad), .groups = "drop") %>%
  slice_max(Total, n = 5) %>%
  pull(Genero)

top5_generos  # revisar cuáles son antes de excluirlos

# 31. Filtrar el dataset, sacando esos 5 géneros ------------------------------
gtdbtk_sin_top5_genero <- gtdbtk_sin_rejects %>%
  filter(!Genero %in% top5_generos)

nrow(gtdbtk_sin_rejects)          # total antes de sacar los 5 géneros
nrow(gtdbtk_sin_top5_genero)      # total después

# 32. Ranking de Género por Hábitat, sin Rejects y sin el top 5 --------------
ranking_genero_habitat_sin_top5 <- gtdbtk_sin_top5_genero %>%
  filter(!is.na(Habitat)) %>%
  count(Habitat, Genero, name = "Cantidad") %>%
  group_by(Habitat) %>%
  mutate(Porcentaje = round(Cantidad / sum(Cantidad) * 100, 1)) %>%
  arrange(Habitat, desc(Cantidad)) %>%
  ungroup()

ranking_genero_habitat_sin_top5

# 33. Gráfico de barras apiladas horizontales, sin el top 5 de géneros ------
datos_stack_genero_sin_top5 <- preparar_datos_stack(ranking_genero_habitat_sin_top5, "Genero")
p_genero_sin_top5 <- graficar_stack_habitat(
  datos_stack_genero_sin_top5,
  "Composición de Géneros por Hábitat (sin Rejects, sin top 5 géneros)"
)
print(p_genero_sin_top5)

ggsave(filename = "plots/stacked_bars_genus_habitat_sin_top5.png", plot = p_genero_sin_top5, width = 12, height = 5, dpi = 300)
ggsave(filename = "plots/stacked_bars_genus_habitat_sin_top5.pdf", plot = p_genero_sin_top5, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_genus_habitat_sin_top5.svg", plot = p_genero_sin_top5, width = 12, height = 5)
ggsave(filename = "plots/stacked_bars_genus_habitat_sin_top5.eps", plot = p_genero_sin_top5, width = 12, height = 5, device = cairo_ps)

# 34. Excel con cantidad y porcentaje de Género por Hábitat, sin top 5 ------
cantidad_genero_sin_top5 <- ranking_genero_habitat_sin_top5 %>%
  select(Habitat, Genero, Cantidad) %>%
  pivot_wider(names_from = Habitat, values_from = Cantidad, values_fill = 0) %>%
  arrange(Genero)

porcentaje_genero_sin_top5 <- ranking_genero_habitat_sin_top5 %>%
  select(Habitat, Genero, Porcentaje) %>%
  pivot_wider(names_from = Habitat, values_from = Porcentaje, values_fill = 0) %>%
  arrange(Genero)

wb_genero_sin_top5 <- createWorkbook()
addWorksheet(wb_genero_sin_top5, "Cantidad")
writeData(wb_genero_sin_top5, "Cantidad", cantidad_genero_sin_top5)
addWorksheet(wb_genero_sin_top5, "Porcentaje")
writeData(wb_genero_sin_top5, "Porcentaje", porcentaje_genero_sin_top5)
saveWorkbook(wb_genero_sin_top5, "excels/resumen_genero_por_habitat_sin_top5.xlsx", overwrite = TRUE)

# 35. TSV con el ranking completo, sin top 5 géneros ---------------------------
write_tsv(ranking_genero_habitat_sin_top5, "visualize/ranking_genus_habitat_sin_top5.tsv")

save.image("visualize/b_gtfbtk_habitat.RData")