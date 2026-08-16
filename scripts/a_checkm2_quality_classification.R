## ============================================================
## 
## 1ª parte del script
##
## Clasificación de MAGs según calidad (CheckM2)
## Basado en criterios MIMAG (Bowers et al. 2017), con dos
## categorías adicionales: "Calidad aceptable" y "Rejects"
##
## Categorías:
##   Alta calidad       -> Completeness >= 90 & Contamination < 5
##   Calidad media       -> Completeness >= 50 & Contamination < 10
##   Calidad aceptable   -> Completeness >= 50 & Contamination < 15
##   Baja calidad        -> Completeness < 50
##   Rejects              -> todas las que no entran en ninguna categoria anterior
##
##
## 2ª parte del script
##
## Metadatos de muestra + cruce con localización / provincia / hábitat
##
##
## Usa el archivo intermedio "checkm2_clasificado.rds"
##
## Formato de nombre de MAG asumido:
##   B  P  ART  24    1001  .MAGScoT_cleanbin_000001
##   |  |  |    |     |      |
##   |  |  |    |     |      +-- resto del ID (bin generado por MAGScoT)
##   |  |  |    |     +--------- código de muestra (4 dígitos)
##   |  |  |    +--------------- año (2 dígitos -> se le suma 2000)
##   |  |  +-------------------- código de localización (3 letras, ej. ART)
##   |  +----------------------- pascuorum
##   +-------------------------- Bombus



## ============================================================
##
## 1º parte del script
##
## Establecemos nuestro entorno de trabajo y nuestros directorios





# 1. Cargar librerías -------------------------------------------------------


library(openxlsx)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)

# 2. Leer el archivo TSV -----------------------------------------------------
checkm2 <- read_tsv("visualize/checkm2_quality_report.tsv")

colnames(checkm2)

# 3. Clasificar los MAGs ------------------------------------------------------
checkm2_clasificado <- checkm2 %>%
  mutate(
    Categoria = case_when(
      Completeness >= 90 & Contamination < 5  ~ "Alta calidad",
      Completeness >= 50 & Contamination < 10 ~ "Calidad media",
      Completeness >= 50 & Contamination < 15 ~ "Calidad aceptable",
      Completeness < 50  & Contamination < 10 ~ "Baja calidad",
      TRUE                                    ~ "Rejects" 
    )
  )

# 4. Ver un resumen de cuántos MAGs caen en cada categoría -------------------
table(checkm2_clasificado$Categoria)

# 5. Guardar el archivo completo con la nueva columna ------------------------
write_tsv(checkm2_clasificado, "visualize/checkm2_quality_report_classified.tsv")

# 6. Guardar un TSV separado por categoría ------------------------
categorias <- split(checkm2_clasificado, checkm2_clasificado$Categoria)

for (cat in names(categorias)) {
  nombre_archivo <- paste0("MAGs_", gsub(" ", "_", cat), ".tsv")
  write_tsv(categorias[[cat]], nombre_archivo)
}
saveRDS(checkm2_clasificado, "visualize/checkm2_clasificado.rds")


# 7. (Opcional) Obtener solo la lista de nombres de MAGs por categoría -------
listas_ids <- lapply(categorias, function(df) df$Name)

# Guardar como .txt, uno por categoría
#for (cat in names(listas_ids)) {
#  nombre_archivo <- paste0("IDs_", gsub(" ", "_", cat), ".txt")
#  writeLines(listas_ids[[cat]], nombre_archivo)
#}

## ============================================================

## 2ª parte del script

# 8. Leer el archivo intermedio generado por clasificar_MAGs_checkm2.R ------
checkm2_clasificado <- readRDS("visualize/checkm2_clasificado.rds")

# 9. Tabla de correspondencias: localización -> provincia -> hábitat --------
# (viene de codigos_localizaciones_habitats.xlsx)
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

# 10. Extraer metadatos del nombre del MAG ------------------------------------

regex_mag <- "^([A-Z])P([A-Z]{3})(\\d{2})(\\d{4})([A-Z]?)\\.(.+)$"

checkm2_metadatos <- checkm2_clasificado %>%
  mutate(
    Organismo     = str_match(Name, regex_mag)[, 2],
    Localizacion  = str_match(Name, regex_mag)[, 3],
    Anio          = as.numeric(str_match(Name, regex_mag)[, 4]) + 2000,
    CodigoMuestra = str_match(Name, regex_mag)[, 5],
    Sufijo        = str_match(Name, regex_mag)[, 6],   
    BinID         = str_match(Name, regex_mag)[, 7]
  )

# Nombres que no coinciden con el patrón esperado 
sum(is.na(checkm2_metadatos$Localizacion))
checkm2_metadatos %>%
  filter(is.na(Localizacion)) %>%
  pull(Name) %>%
  head(46)

# 11. Cruzar con la tabla de localización / provincia / hábitat --------------
checkm2_metadatos <- checkm2_metadatos %>%
  left_join(tabla_localizaciones, by = "Localizacion")

# Localizaciones que no coinciden con la tabla de referencia
checkm2_metadatos %>% filter(is.na(Habitat)) %>% distinct(Localizacion)

# 12. Resumen de calidad por Hábitat (tabla cruzada) --------------------------
resumen_habitat <- checkm2_metadatos %>%
  count(Habitat, Categoria) %>%
  pivot_wider(names_from = Categoria, values_from = n, values_fill = 0)

resumen_habitat

# 13. Exportar todo a un excel y tsv ------------------------------------

wb <- createWorkbook()

addWorksheet(wb, "Todos_metadatos")
writeData(wb, "Todos_metadatos", checkm2_metadatos)

addWorksheet(wb, "Resumen_Habitat")
writeData(wb, "Resumen_Habitat", resumen_habitat)

saveWorkbook(wb, "excels/MAGs_metadatos_habitat.xlsx", overwrite = TRUE)

# TSV de cada tabla por separado
write_tsv(checkm2_metadatos, "visualize/MAGs_metadatos_habitat.tsv")
write_tsv(resumen_habitat, "visualize/Resumen_Habitat.tsv")

# 14. Diagramas de quesitos (pie charts) por Hábitat --------------------------
# Uno por cada hábitat (Natural / Agricola / Urbano), mostrando en cada
# "quesito" la categoría, la cantidad absoluta y el porcentaje.

datos_pie <- checkm2_metadatos %>%
  filter(!is.na(Habitat)) %>%
  count(Habitat, Categoria) %>%
  group_by(Habitat) %>%
  mutate(
    Porcentaje = n / sum(n) * 100,
    Etiqueta   = paste0(Categoria, "\n", n, " (", round(Porcentaje, 1), "%)")
  ) %>%
  ungroup()

habitats <- unique(datos_pie$Habitat)

for (h in habitats) {
  
  datos_h <- datos_pie %>% filter(Habitat == h)
  
  p <- ggplot(datos_h, aes(x = "", y = n, fill = Categoria)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    geom_text(aes(label = Etiqueta),
              position = position_stack(vjust = 0.5), size = 2.5) +
    labs(title = paste("Calidad de MAGs -", h),
         x = NULL, y = NULL, fill = "Categoría") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  
  #ggsave(filename = paste0("plots/pie_calidad_", h, ".png"), plot = p, width = 6, height = 6, dpi = 300)
  #ggsave(filename = paste0("plots/pie_calidad_", h, ".pdf"), plot = p, width = 6, height = 6)
  #ggsave(filename = paste0("plots/pie_calidad_", h, ".svg"), plot = p, width = 6, height = 6)
  #ggsave(filename = paste0("plots/pie_calidad_", h, ".eps"), plot = p, width = 6, height = 6, device = cairo_ps)
}

# 14b. Un único gráfico con los tres quesitos juntos  --------------
p_conjunto <- ggplot(datos_pie, aes(x = "", y = n, fill = Categoria)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  geom_text(aes(label = Etiqueta),
            position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~ Habitat) +
  labs(title = "Calidad de MAGs por Hábitat",
       x = NULL, y = NULL, fill = "Categoría") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_conjunto)

  #ggsave(filename = "plots/pie_calidad_por_habitat.png", plot = p_conjunto, width = 12, height = 5, dpi = 300)
  #ggsave(filename = "plots/pie_calidad_por_habitat.pdf", plot = p_conjunto, width = 12, height = 5)
  #ggsave(filename = "plots/pie_calidad_por_habitat.svg", plot = p_conjunto, width = 12, height = 5)
  #ggsave(filename = "plots/pie_calidad_por_habitat.eps", plot = p_conjunto, width = 12, height = 5, device = cairo_ps)


#save.image("visualize/a_checkm2_classification.RData")
