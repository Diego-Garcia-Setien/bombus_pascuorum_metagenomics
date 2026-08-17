# Este es un script para realizar estimaciones corporales a partir
# de la distancia intertegular mediante el comando bodysize () de el 
# paquete pollimetry, en las muestras del año 2024 y 2025.
#
#
# Este script sirve para realizar distintos análisis
#
# Para ambos años por separado: 
#
#   - Un modelo lineal mediante el comando lm(Est.Weight ~ Habitat),
#     para testear si el peso estimando varía según el hábitat (Natural/Agrícola/Urbano)
#   - Un ANOVA para obtener la significancia del modelo.
#   - Se realizan gráficos de los efectos estimados por hábitat, con media y
#     desviación estandar.
#
# Para la combinación de ambos años:
#
#   - Se unen ambos datasets de 2024 y 2025 mediante (bind_rows).
#   - Se realiza un modelo de interacción mediante lm(Est.Weight ~ Habitat * Año) para
#     para conocer si el efecto del hábitat sobre el peso depende del año (o viceversa).
#     Se realizan gráficos de los efectos estimados por hábitat y año, con media y
#     desviación estandar.
#
# Realiza análisis por provincia
#
#   - Para ambos años por separado, mediante lm(Est.Weight ~ Provincias * Habitat)
#     se testea si el peso varía según la interacción entre provincia y hábitat.Se realizan
#     los gráficos de efectos con hábitat en el eje X y provincia como color, uno por año,
#     con la misma escala en el eje Y (limits=c(40,150)) para poder compararlos visualmente, 
#     combinados al final en una sola figura (plot24 + plot25).
#
#
#
library(pollimetry)
library(pollimetrydata)
library(dplyr)
library(stringr)
library(readxl)
library(writexl)
library(effects)
library(ggplot2)
library(tidyverse)
library(car)
library(ggeffects)
library(patchwork)
library(ggstatsplot)

#data25
B_pascuorum_samples25 <- read_excel("B.pascuorum_samples_metadata.xlsx")
#View(B_pascuorum_samples25)

#Cambiamos el nombre de intertegular distance a ITD, omitimos las final sin valores (NA),
# cambiamos a region Europa

B_pascuorum_ITD25 <-B_pascuorum_samples25 %>% 
    rename(ITD = `Intertegular_distance (mm)`,
           Sex = sex) %>%
    mutate(Species = str_replace_all(Species, " " , "_")) %>%
    mutate(Region="Europe") %>%
    select(Code,Locality, Provincias,Habitat, Species, Sex, Region, ITD)

B_pascuorum_ITD25 <- na.omit(B_pascuorum_ITD25)
B_pascuorum_ITD25$ITD<- as.numeric(B_pascuorum_ITD25$ITD)

set.seed(42)  # <-- para que bodysize() dé siempre el mismo resultado
body_size.B_pascuorum25 <- bodysize(
  x = B_pascuorum_ITD25,
  taxa = "bee",
  type = "ITD"
)
body_size.B_pascuorum25$Provincias[body_size.B_pascuorum25$Provincias == "Navarra"]<-"Gipuzkoa"

body_size.B_pascuorum25$Provincias[body_size.B_pascuorum25$Provincias == "Guipuzkoa"]<-"Gipuzkoa"
write_xlsx(body_size.B_pascuorum25,"body_size_B.pascuorum_25.xlsx")



datos_bp_bs25 <- read_excel("body_size_B.pascuorum_25.xlsx")

datos_bp_bs25$Año<-"2025"

#Análisis25

mod25=lm(Est.Weight~Habitat, data = datos_bp_bs25)
summary(mod25)
plot(mod25)
Anova(mod25)

plot(allEffects(mod25))

#Sacar la media
datos_bp_bs25 %>% 
  group_by(Habitat) %>%
  summarise(mediaBody=mean(Est.Weight), 
            error=sd(Est.Weight))

#datos24

B_pascuorum_samples24 <- read_excel("B.pascuorum_samples_metadata.xlsx", 
                                    sheet = "data_2024")
#View(B_pascuorum_samples24)

#Cambiamos el nombre de intertegular distance a ITD, omitimos las final sin valores (NA),
# cambiamos a region Europa

B_pascuorum_ITD24 <-B_pascuorum_samples24 %>% 
  rename(ITD = `Intertegular_distance (mm)`,
         Sex = sex) %>%
  mutate(Species = str_replace_all(Species, " " , "_")) %>%
  mutate(Region="Europe") %>%
  select(Code,Locality, Provincias,Habitat, Species, Sex, Region, ITD)

B_pascuorum_ITD24 <- na.omit(B_pascuorum_ITD24)
B_pascuorum_ITD24$ITD<- as.numeric(B_pascuorum_ITD24$ITD)

set.seed(42)  # <-- para que bodysize() dé siempre el mismo resultado
body_size.B_pascuorum24 <- bodysize(
  x = B_pascuorum_ITD24,
  taxa = "bee",
  type = "ITD"
)
body_size.B_pascuorum24$Provincias[body_size.B_pascuorum24$Provincias == "Navarra"]<-"Gipuzkoa"

body_size.B_pascuorum24$Provincias[body_size.B_pascuorum24$Provincias == "Guipuzkoa"]<-"Gipuzkoa"

write_xlsx(body_size.B_pascuorum24,"body_size_B.pascuorum_24.xlsx")



datos_bp_bs24 <- read_excel("body_size_B.pascuorum_24.xlsx")

datos_bp_bs24$Año<-"2024"

#Gráfico de modelo 25
eff25 <- ggpredict(mod25, terms = "Habitat")
plot2025 <- ggplot(eff25, aes(x = x, y = predicted, color = x)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_color_manual(values = c("green3","gold3","blue")) +
  labs(x = "Hábitats", y = NULL, color = "Hábitat", title = "2025") +
  theme_classic()

#Análisis24
mod24 <- lm(Est.Weight~Habitat, data = datos_bp_bs24)
summary(mod24)
plot(mod24)
Anova(mod24)
plot(allEffects(mod24))

#Sacar la media
datos_bp_bs24 %>% 
  group_by(Habitat) %>%
  summarise(mediaBody=mean(Est.Weight), 
            error=sd(Est.Weight))

#Gráfico de modelo 24
eff24 <- ggpredict(mod24, terms = "Habitat")
plot2024 <- ggplot(eff24, aes(x = x, y = predicted, color = x)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_color_manual(values = c("green3","gold3","blue")) +
  labs(x = "Hábitats", y = "Peso estimado", color = "Hábitat", title = "2024") +
  theme_classic()

# Combinar los dos gráficos con una única leyenda compartida
plot2024 + plot2025 + plot_layout(guides = "collect")

#Grafico 24 y 25 



#Juntar 2025 + 2024

datos_unidos <- bind_rows (datos_bp_bs25, datos_bp_bs24)

mod2524=lm(Est.Weight~Habitat * Año, data = datos_unidos)
summary(mod2524)

plot(mod2524)

Anova(mod2524)

#Sacar la media de los dos años
datos_unidos %>% 
  group_by(Habitat) %>%
  summarise(mediaBody=mean(Est.Weight), 
            error=sd(Est.Weight))

datos_unidos %>% 
  group_by(Año) %>%
  summarise(mediaBody=mean(Est.Weight), 
            error=sd(Est.Weight))

#Gráfico de modelo de los dos años

eff2524<- ggpredict(mod2524, terms = c ("Habitat","Año"))
ggplot(eff2524, aes(x = x, y = predicted, color = group)) + 
  geom_point(size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                                       width = 0.2) +
  scale_color_manual(values = c("green3","gold3"))+
  labs(x = "Hábitats", y = "Peso estimado", color = "Año")+ theme_classic()

plot (allEffects(mod2524))




#Análsis entre provincias

# Año 24

modProvincias24=lm(Est.Weight~Provincias * Habitat, data = body_size.B_pascuorum24)
summary(modProvincias24)

plot(modProvincias24)
Anova(modProvincias24)


plot (allEffects(modProvincias24))

# Año 25
modProvincias25=lm(Est.Weight~Provincias * Habitat, data = body_size.B_pascuorum25)
summary(modProvincias25)

plot(modProvincias25)
Anova(modProvincias25)

plot (allEffects(modProvincias25))


# Gráfico año 2025

eff25_Provincia<- ggpredict(modProvincias25, terms = c ("Habitat","Provincias"))
plot25=ggplot(eff25_Provincia, aes(x = x, y = predicted, color = group)) + 
  geom_point(size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                                       width = 0.2) +
  scale_color_manual(values = c("green3","gold3","blue"))+ scale_y_continuous(limits=c(40,150))+
  labs(title= "2025", x = "Hábitats", y = "", color = "Provincias")+ theme_classic()


# Gráfico año 2024

eff24_Provincias<- ggpredict(modProvincias24, terms = c ("Habitat","Provincias"))
plot24=ggplot(eff24_Provincias, aes(x = x, y = predicted, color = group)) + 
  geom_point(size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                                       width = 0.2) +
  scale_color_manual(values = c("green3","gold3","blue"))+ scale_y_continuous(limits=c(40,150))+
  labs(title="2024",x = "Hábitats", y = "Peso estimado", color = "Provincias")+ theme_classic()+theme(legend.position = "none")

plot24+plot25

# Otro modelo de gráfico para el análisis de Hábitat x Provincia x Año

modProvinciaHabitatAño <- lm(Est.Weight ~ Provincias * Habitat * Año, data = datos_unidos)
summary(modProvinciaHabitatAño)
plot(modProvinciaHabitatAño)
Anova(modProvinciaHabitatAño)

plot(allEffects(modProvinciaHabitatAño))

# Medias por combinación de Provincia y Hábitat (ambos años juntos)
datos_unidos %>%
  group_by(Provincias, Habitat) %>%
  summarise(mediaBody = mean(Est.Weight),
            error = sd(Est.Weight))

# Medias por Provincia, Hábitat Y Año (para ver si el patrón cambia entre años)
datos_unidos %>%
  group_by(Provincias, Habitat, Año) %>%
  summarise(mediaBody = mean(Est.Weight),
            error = sd(Est.Weight))

# Gráfico de efectos: Hábitat en el eje X, color = Provincia, un panel por Año
eff_provincia_habitat_año <- ggpredict(modProvinciaHabitatAño,
                                       terms = c("Habitat", "Provincias", "Año"))

ggplot(eff_provincia_habitat_año, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap(~ facet) +
  scale_color_manual(values = c("green3", "gold3", "blue")) +
  labs(x = "Hábitats", y = "Peso estimado", color = "Provincias") +
  theme_classic()

# Análisis de Provincia y Hábitat con los dos años combinados como un solo dataset ---

modProvinciaHabitat_combinado <- lm(Est.Weight ~ Provincias * Habitat, data = datos_unidos)
summary(modProvinciaHabitat_combinado)
plot(modProvinciaHabitat_combinado)
Anova(modProvinciaHabitat_combinado)

plot(allEffects(modProvinciaHabitat_combinado))

# Medias por combinación de Provincia y Hábitat (2024+2025 juntos)
datos_unidos %>%
  group_by(Provincias, Habitat) %>%
  summarise(mediaBody = mean(Est.Weight),
            error = sd(Est.Weight))

# Gráfico de efectos: Hábitat en el eje X, color = Provincia (un solo gráfico, sin separar por año)
eff_provincia_habitat_combinado <- ggpredict(modProvinciaHabitat_combinado,
                                             terms = c("Habitat", "Provincias"))

ggplot(eff_provincia_habitat_combinado, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c("green3", "gold3", "blue")) +
  labs(x = "Hábitats", y = "Peso estimado", color = "Provincias") +
  theme_classic()

# Test de normalidad mediante Shapiro-Wilk del peso estimado por hábitat

# 2024
datos_bp_bs24 %>%
  group_by(Habitat) %>%
  summarise(p_valor_shapiro = shapiro.test(Est.Weight)$p.value)

# 2025
datos_bp_bs25 %>%
  group_by(Habitat) %>%
  summarise(p_valor_shapiro = shapiro.test(Est.Weight)$p.value)

# Combinado (2024 + 2025)
datos_unidos %>%
  group_by(Habitat) %>%
  summarise(p_valor_shapiro = shapiro.test(Est.Weight)$p.value)

# Test de normalidad mediante Shapirp_Wilk de los residuos del modelo de habitats

# 2024
shapiro.test(residuals(mod24))

# 2025
shapiro.test(residuals(mod25))

# Combinado (2024 + 2025)
shapiro.test(residuals(mod2524))

# Gráficos de habitats mediante ggbetweenstats, no paramétricos ya que los datos no superan
# el test de normalidad de Shapiro-Wilk

ggbetweenstats(
  data = datos_bp_bs24, x = Habitat, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Hábitat (2024)"
)

ggbetweenstats(
  data = datos_bp_bs25, x = Habitat, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Hábitat (2025)"
)

ggbetweenstats(
  data = datos_unidos, x = Habitat, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Hábitat (2024 + 2025)"
)

# Test de normalidad mediante Shapirp_Wilk de los residuos del modelo de provincias

# Modelos de un solo factor: Provincias
mod_prov24 <- lm(Est.Weight ~ Provincias, data = datos_bp_bs24)
mod_prov25 <- lm(Est.Weight ~ Provincias, data = datos_bp_bs25)
mod_prov_combinado <- lm(Est.Weight ~ Provincias, data = datos_unidos)

# Shapiro-Wilk sobre los residuos de cada uno del modelo de provincias
shapiro.test(residuals(mod_prov24))
shapiro.test(residuals(mod_prov25))
shapiro.test(residuals(mod_prov_combinado))

# Gráficos de provincias mediante ggbetweenstats, no paramétricos ya que los datos no superan
# el test de normalidad de Shapiro-Wilk

ggbetweenstats(
  data = datos_bp_bs24, x = Provincias, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Provincia (2024)"
)

ggbetweenstats(
  data = datos_bp_bs25, x = Provincias, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Provincia (2025)"
)

ggbetweenstats(
  data = datos_unidos, x = Provincias, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Provincia (2024 + 2025)"
)


# Shapiro-Wilk sobre los residuos de cada uno del modelo de provincias x habitat


shapiro.test(residuals(modProvincias24))
shapiro.test(residuals(modProvincias25))
shapiro.test(residuals(modProvinciaHabitat_combinado))


# Gráficos de provincias x hábitat mediante ggbetweenstats, no paramétricos ya que los datos no superan
# el test de normalidad de Shapiro-Wilk


datos_bp_bs24 <- datos_bp_bs24 %>%
  mutate(Provincia_Habitat = paste(Provincias, Habitat, sep = " - "))

datos_bp_bs25 <- datos_bp_bs25 %>%
  mutate(Provincia_Habitat = paste(Provincias, Habitat, sep = " - "))

datos_unidos <- datos_unidos %>%
  mutate(Provincia_Habitat = paste(Provincias, Habitat, sep = " - "))

ggbetweenstats(
  data = datos_bp_bs24, x = Provincia_Habitat, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Provincia x Hábitat (2024)"
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggbetweenstats(
  data = datos_bp_bs25, x = Provincia_Habitat, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Provincia x Hábitat (2025)"
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggbetweenstats(
  data = datos_unidos, x = Provincia_Habitat, y = Est.Weight,
  type = "nonparametric", pairwise.display = "significant",
  title = "Peso estimado por Provincia x Hábitat (2024 + 2025)"
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
