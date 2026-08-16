## ============================================================
## Nivel funcional/metabólico ponderado por abundancia de especie,
## separado en los DOS TIPOS de datos que trae DRAM:
##
##   TIPO 1 - Rutas metabólicas (32 columnas numéricas 0-1):
##            fracción de completitud de la ruta.
##            El ponderado da la "completitud promedio" de la ruta
##            en la comunidad de ese hábitat/provincia.
##
##   TIPO 2 - CAZy / genes marcadores (66 columnas True/False):
##            presencia/ausencia real del gen.
##            El ponderado da directamente el "% de la comunidad que
##            tiene el gen", pesando por cuán común es cada especie.
##
## Para cada tipo, se hace el análisis por HÁBITAT y por PROVINCIA.
##
## ============================================================

# 1. Cargar librerías -------------------------------------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(openxlsx)
library(ggplot2)

# 2. Leer el archivo de DRAM y el de taxonomía sin Rejects -------------------
dram <- read_tsv("visualize/DRAMR_product.tsv")
gtdbtk_sin_rejects <- readRDS("visualize/gtdbtk_sin_rejects.rds")

dim(dram)  # 20 filas x 99 columnas (genome + 98 funciones)

# 3. Separar las columnas funcionales en los dos tipos -----------------------
# read_tsv ya detecta el tipo de cada columna al leer el archivo:
# "logical" = True/False (CAZy/genes), "numeric" = fracción 0-1 (rutas)
tipos_columna <- sapply(dram, class)

columnas_rutas     <- names(tipos_columna)[tipos_columna == "numeric"]
columnas_booleanas <- names(tipos_columna)[tipos_columna == "logical"]

# Sacar las columnas de CAZy del análisis booleano: en este dataset no
# se anotaron correctamente (quedan todas vacías/FALSE) y solo ocupan
# espacio sin aportar información real.
columnas_cazy <- columnas_booleanas[str_starts(columnas_booleanas, "CAZy:")]
columnas_booleanas <- setdiff(columnas_booleanas, columnas_cazy)

length(columnas_cazy)       # cuántas columnas de CAZy se sacaron
columnas_cazy                # revisar cuáles son, por si alguna sí tenía datos

length(columnas_rutas)      # debería dar 32
length(columnas_booleanas)  # debería dar 66 - 19 (CAZy) = 47

# Convertir booleanas a 0/1 para poder promediarlas igual que las numéricas
dram <- dram %>%
  mutate(across(all_of(columnas_booleanas), ~ as.numeric(.)))

# 4. Averiguar a qué ESPECIE pertenece el MAG representante de DRAM ---------
mapa_especie <- gtdbtk_sin_rejects %>%
  distinct(Name, Especie_final)

dram_especie <- dram %>%
  left_join(mapa_especie, by = c("genome" = "Name"))

# Chequeo: MAGs de DRAM sin especie asignada (no deberían aparecer)
dram_especie %>% filter(is.na(Especie_final)) %>% pull(genome)

##=======================================================
######################################################||#
## Funciones reutilizables                            ||#
######################################################||#
##=======================================================

# Pasa un subconjunto de columnas funcionales a formato largo
pasar_a_largo <- function(dram_especie, columnas) {
  dram_especie %>%
    select(Especie_final, all_of(columnas)) %>%
    pivot_longer(cols = all_of(columnas), names_to = "Funcion", values_to = "Valor")
}

# Calcula el nivel ponderado por abundancia de cada función, por grupo
calcular_nivel_ponderado <- function(gtdbtk_sin_rejects, dram_largo, columna_grupo) {
  
  abundancia_especie <- gtdbtk_sin_rejects %>%
    filter(!is.na(.data[[columna_grupo]])) %>%
    count(.data[[columna_grupo]], Especie_final, name = "Cantidad") %>%
    rename(Grupo = all_of(columna_grupo)) %>%
    group_by(Grupo) %>%
    mutate(Proporcion = Cantidad / sum(Cantidad)) %>%
    ungroup()
  
  cruce <- abundancia_especie %>%
    inner_join(dram_largo, by = "Especie_final", relationship = "many-to-many")
  
  cobertura <- cruce %>%
    distinct(Grupo, Especie_final, Proporcion) %>%
    group_by(Grupo) %>%
    summarise(Cobertura_DRAM = round(100 * sum(Proporcion), 1))
  
  cruce <- cruce %>%
    group_by(Grupo, Funcion) %>%
    mutate(Peso = Proporcion / sum(Proporcion)) %>%
    ungroup()
  
  resultado <- cruce %>%
    group_by(Grupo, Funcion) %>%
    summarise(Nivel_ponderado_pct = round(100 * sum(Peso * Valor, na.rm = TRUE), 1), .groups = "drop")
  
  list(resultado = resultado, cobertura = cobertura)
}

graficar_heatmap_ponderado <- function(resumen, titulo) {
  ggplot(resumen, aes(x = Grupo, y = Funcion, fill = Nivel_ponderado_pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(Nivel_ponderado_pct, "%")), size = 2, color = "black") +
    scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0, 100)) +
    labs(title = titulo, x = NULL, y = NULL, fill = "% ponderado\npor abundancia") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 6))
}

graficar_barras_top_ponderado <- function(resumen, top_n = 15, titulo) {
  rangos <- resumen %>%
    group_by(Funcion) %>%
    summarise(Rango = max(Nivel_ponderado_pct) - min(Nivel_ponderado_pct)) %>%
    filter(Rango > 0)   # descarta funciones sin ninguna variación entre grupos
  
  funciones_top <- rangos %>%
    slice_max(Rango, n = top_n, with_ties = FALSE) %>%
    pull(Funcion)
  
  datos_plot <- resumen %>%
    filter(Funcion %in% funciones_top) %>%
    left_join(rangos, by = "Funcion")
  
  ggplot(datos_plot, aes(x = reorder(Funcion, Rango), y = Nivel_ponderado_pct, fill = Grupo)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    labs(title = titulo, x = NULL, y = "% ponderado por abundancia", fill = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Envuelve todo el proceso: calcular + graficar + exportar, para un tipo
# de columnas (rutas o booleanas) y un grupo (Habitat o Provincia)
analizar_y_exportar <- function(gtdbtk_sin_rejects, dram_especie, columnas, columna_grupo,
                                etiqueta_tipo, top_n = 15) {
  
  dram_largo <- pasar_a_largo(dram_especie, columnas)
  resultado  <- calcular_nivel_ponderado(gtdbtk_sin_rejects, dram_largo, columna_grupo)
  
  cat("\n===== Cobertura DRAM:", etiqueta_tipo, "x", columna_grupo, "=====\n")
  print(resultado$cobertura)
  
  resumen <- resultado$resultado
  
  p_heatmap <- graficar_heatmap_ponderado(
    resumen, paste0(etiqueta_tipo, " ponderado por abundancia, por ", columna_grupo)
  )
  print(p_heatmap)
  
  p_barras <- graficar_barras_top_ponderado(
    resumen, top_n = top_n,
    titulo = paste0(etiqueta_tipo, ": mayor diferencia entre ", columna_grupo, "s")
  )
  print(p_barras)
  
  prefijo <- paste0(tolower(gsub(" ", "_", etiqueta_tipo)), "_", tolower(columna_grupo))
  
 #ggsave(paste0("plots/heatmap_", prefijo, ".png"), p_heatmap, width = 8, height = max(6, length(columnas) * 0.35), dpi = 300)
 #ggsave(paste0("plots/heatmap_", prefijo, ".pdf"), p_heatmap, width = 8, height = max(6, length(columnas) * 0.35))
 #ggsave(paste0("plots/barras_top_", prefijo, ".png"), p_barras, width = 9, height = 7, dpi = 300)
 #ggsave(paste0("plots/barras_top_", prefijo, ".pdf"), p_barras, width = 9, height = 7)
  
  wb <- createWorkbook()
  addWorksheet(wb, "Nivel_ponderado")
  writeData(wb, "Nivel_ponderado", resumen %>% pivot_wider(names_from = Grupo, values_from = Nivel_ponderado_pct))
  addWorksheet(wb, "Cobertura_DRAM")
  writeData(wb, "Cobertura_DRAM", resultado$cobertura)
  #saveWorkbook(wb, paste0("visualize/dram_", prefijo, ".xlsx"), overwrite = TRUE)
  
  write_tsv(resumen, paste0("visualize/dram_", prefijo, ".tsv"))
  write_tsv(resultado$cobertura, paste0("visualize/dram_cobertura_", prefijo, ".tsv"))
  
  resumen
}

##=======================================================
######################################################||#
## TIPO 1: Rutas metabólicas (numéricas 0-1)           ||#
######################################################||#
##=======================================================

# Rutas x Hábitat -------------------------------------------------------------
resumen_rutas_habitat <- analizar_y_exportar(
  gtdbtk_sin_rejects, dram_especie, columnas_rutas, "Habitat",
  etiqueta_tipo = "Rutas metabolicas"
)

# Rutas x Provincia -------------------------------------------------------------
resumen_rutas_provincia <- analizar_y_exportar(
  gtdbtk_sin_rejects, dram_especie, columnas_rutas, "Provincia",
  etiqueta_tipo = "Rutas metabolicas"
)

##=======================================================
######################################################||#
## TIPO 2: CAZy / genes marcadores (True/False)        ||#
######################################################||#
##=======================================================

# CAZy/genes x Hábitat -------------------------------------------------------------
resumen_booleanas_habitat <- analizar_y_exportar(
  gtdbtk_sin_rejects, dram_especie, columnas_booleanas, "Habitat",
  etiqueta_tipo = "CAZy y genes", top_n = 20
)

# CAZy/genes x Provincia -------------------------------------------------------------
resumen_booleanas_provincia <- analizar_y_exportar(
  gtdbtk_sin_rejects, dram_especie, columnas_booleanas, "Provincia",
  etiqueta_tipo = "CAZy y genes", top_n = 20
)
