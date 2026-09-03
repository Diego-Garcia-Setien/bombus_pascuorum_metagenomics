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
library(ragg)

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

# Reemplazar la función existente para aumentar el tamaño de texto:
graficar_heatmap_taxa <- function(datos_largos, columna_taxa, titulo) {
  ggplot(datos_largos, aes(x = .data[[columna_taxa]], y = Funcion, fill = Valor)) +
    geom_tile(color = "gray90", linewidth = 0.3) +
    scale_fill_gradient(
      low = "whitesmoke", 
      high = "darkgreen", 
      limits = c(0, 1), 
      labels = scales::percent
    ) +
    labs(
      title = titulo,
      x = NULL,
      y = NULL,
      fill = "Completitud /\nPresencia"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15, face = "italic", color = "black"), # Géneros/Especies
      axis.text.y = element_text(size = 15, color = "black"), # Rutas y Genes
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right"
    )
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



##=======================================================
######################################################||#
## Perfil metabólico por ESPECIE y por GÉNERO (Separados)
######################################################||#
##=======================================================

# 1. Asegurar la extracción correcta de Género y Especie --------------------
# Nota: Si en tu dataset gtdbtk_sin_rejects existe la columna 'Genero',
# puedes usarla directamente. Si no, extraemos la primera palabra de Especie_final.
dram_taxonomia <- dram_especie %>%
  mutate(
    Especie = Especie_final,
    Genero = if("Genero" %in% colnames(gtdbtk_sin_rejects)) {
      gtdbtk_sin_rejects$Genero[match(genome, gtdbtk_sin_rejects$Name)]
    } else {
      str_extract(Especie_final, "^[A-Za-z0-9_-]+") # Extrae la primera palabra/género
    }
  )

# 2. Función genérica para graficar Heatmaps por Nivel Taxonómico --------------
graficar_heatmap_taxa <- function(datos_largos, columna_taxa, titulo) {
  ggplot(datos_largos, aes(x = .data[[columna_taxa]], y = Funcion, fill = Valor)) +
    geom_tile(color = "gray90", linewidth = 0.2) +
    scale_fill_gradient(
      low = "whitesmoke", 
      high = "darkgreen", 
      limits = c(0, 1), 
      labels = scales::percent
    ) +
    labs(
      title = titulo,
      x = NULL,
      y = NULL,
      fill = "Completitud"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "italic", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right"
    )
}

##=======================================================
######################################################||#
## Preparación de datos y extracción de Taxonomía     ||#
######################################################||#
##=======================================================
  
# 1. Crear dram_taxonomia extrayendo Género y Especie
dram_taxonomia <- dram_especie %>%
  mutate(
    Especie = Especie_final,
    Genero = if ("Genero" %in% colnames(gtdbtk_sin_rejects)) {
      gtdbtk_sin_rejects$Genero[match(genome, gtdbtk_sin_rejects$Name)]
    } else {
      str_extract(Especie_final, "^[A-Za-z0-9_-]+")
    }
  )

# Corregir el nombre con doble 's' en la lista ordenada
orden_generos <- c(
  "Snodgrassella", "Gilliamella", "Bombiscardovia", "Xylocopilactobacillus",
  "Lactobacillus", "Fructobacillus", "Microbacterium", "Pantoea",
  "Entomomonas", "Rahnella", "Saccharibacter", "Sediminibacterium"
)

# 2. Asegurar la extracción limpia del género eliminando 'g__' o prefijos de GTDB
dram_taxonomia <- dram_especie %>%
  mutate(
    Especie = Especie_final,
    Genero = if ("Genero" %in% colnames(gtdbtk_sin_rejects)) {
      gtdbtk_sin_rejects$Genero[match(genome, gtdbtk_sin_rejects$Name)]
    } else {
      str_extract(Especie_final, "^[A-Za-z0-9_-]+")
    },
    # Eliminar prefijos de GTDB (ej. 'g__Snodgrassella' -> 'Snodgrassella')
    Genero = str_remove(Genero, "^g__"),
    # Eliminar variantes/subíndices de GTDB (ej. 'Snodgrassella_A' -> 'Snodgrassella')
    Genero = str_remove(Genero, "_[A-Z]$") 
  )


## -----------------------------------------------------------------------------
## A) ANÁLISIS A NIVEL DE ESPECIE
## -----------------------------------------------------------------------------

dram_especie_largo <- dram_taxonomia %>%
  pivot_longer(
    cols = c(all_of(columnas_rutas), all_of(columnas_booleanas)),
    names_to = "Funcion",
    values_to = "Valor"
  )

p_rutas_especie <- dram_especie_largo %>%
  filter(Funcion %in% columnas_rutas) %>%
  graficar_heatmap_taxa(columna_taxa = "Especie", titulo = "Rutas Metabólicas por Especie")

p_genes_especie <- dram_especie_largo %>%
  filter(Funcion %in% columnas_booleanas) %>%
  graficar_heatmap_taxa(columna_taxa = "Especie", titulo = "Genes Marcadores y CAZy por Especie")

#print(p_rutas_especie)
#print(p_genes_especie)

##=======================================================
######################################################||#
## B) ANÁLISIS A NIVEL DE GÉNERO (Limpio, sin filas en 0)
######################################################||#
##=======================================================

dram_genero_resumen <- dram_taxonomia %>%
  filter(!is.na(Genero), Genero != "NA") %>%
  group_by(Genero) %>%
  summarise(
    across(all_of(columnas_rutas), ~ round(mean(.x, na.rm = TRUE), 2)),
    across(all_of(columnas_booleanas), ~ if_else(any(.x == 1, na.rm = TRUE), 1, 0)),
    .groups = "drop"
  ) %>%
  mutate(Genero = factor(Genero, levels = orden_generos)) %>%
  filter(!is.na(Genero)) %>%
  arrange(Genero)

# Pasar a formato largo Y FILTRAR funciones vacías
dram_genero_largo <- dram_genero_resumen %>%
  pivot_longer(
    cols = -Genero,
    names_to = "Funcion",
    values_to = "Valor"
  ) %>%
  group_by(Funcion) %>%
  filter(sum(Valor, na.rm = TRUE) > 0) %>% # Retira las filas completamente vacías/en cero
  ungroup()

# 4. Heatmap: Rutas por Género
p_rutas_genero <- dram_genero_largo %>%
  filter(Funcion %in% columnas_rutas) %>%
  graficar_heatmap_taxa(columna_taxa = "Genero", titulo = "Rutas Metabólicas por Género")

# 5. Heatmap: Genes Marcadores por Género
# 5. Heatmap: Genes Marcadores por Género (con leyenda en 2 bloques discretos)
p_genes_genero <- dram_genero_largo %>%
  filter(Funcion %in% columnas_booleanas) %>%
  # Convertimos el Valor en factor para que ggplot dibuje bloques discretos
  mutate(Valor = factor(Valor, levels = c(1, 0), labels = c("Presente", "Ausente"))) %>%
  ggplot(aes(x = Genero, y = Funcion, fill = Valor)) +
  geom_tile(color = "gray90", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Ausente" = "whitesmoke", "Presente" = "darkgreen"),
    drop = FALSE
  ) +
  labs(
    title = "Genes Marcadores por Género",
    x = NULL,
    y = NULL,
    fill = "Estado"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "italic", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )
print(p_rutas_genero)
print(p_genes_genero)

## -----------------------------------------------------------------------------
## C) EXPORTACIÓN DE GRÁFICOS Y EXCEL
## -----------------------------------------------------------------------------

#ggsave("plots/heatmap_rutas_especie.png", p_rutas_especie, width = 10, height = 8, dpi = 300)
#ggsave("plots/heatmap_genes_especie.png", p_genes_especie, width = 10, height = 10, dpi = 300)
ggsave("plots/heatmap_rutas_genero.png", p_rutas_genero, width = 12, height = 10, dpi = 300, device = ragg::agg_png)
ggsave("plots/heatmap_genes_genero.png", p_genes_genero, width = 12, height = 10, dpi = 300, device = ragg::agg_png)

wb_taxa <- createWorkbook()
addWorksheet(wb_taxa, "Perfil_Especie")
writeData(wb_taxa, "Perfil_Especie", dram_taxonomia %>% select(-Genero))
addWorksheet(wb_taxa, "Perfil_Genero")
writeData(wb_taxa, "Perfil_Genero", dram_genero_resumen)
#saveWorkbook(wb_taxa, "visualize/dram_perfil_taxonomico.xlsx", overwrite = TRUE)

## -----------------------------------------------------------------------------
## D) EXPORTACIÓN FINAL EN 4 ARCHIVOS TSV
## -----------------------------------------------------------------------------

# 1. Rutas metabólicas por ESPECIE
rutas_especie <- dram_taxonomia %>% select(Especie, all_of(columnas_rutas))
#write_tsv(rutas_especie, "visualize/dram_rutas_por_especie.tsv")

# 2. Genes marcadores por ESPECIE
genes_especie <- dram_taxonomia %>% select(Especie, all_of(columnas_booleanas))
#write_tsv(genes_especie, "visualize/dram_genes_por_especie.tsv")

# 3. Rutas metabólicas por GÉNERO
rutas_genero <- dram_genero_resumen %>% select(Genero, all_of(columnas_rutas))
#write_tsv(rutas_genero, "visualize/dram_rutas_por_genero.tsv")

# 4. Genes marcadores por GÉNERO
genes_genero <- dram_genero_resumen %>% select(Genero, all_of(columnas_booleanas))
#write_tsv(genes_genero, "visualize/dram_genes_por_genero.tsv")
