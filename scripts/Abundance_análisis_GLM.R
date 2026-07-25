library(dplyr)    # Manipulación de datos
library(tidyr)    # Limpieza de datos  
library(ggplot2)  # Visualización básica
library(MASS)     # glm.nb (binomial negativo)
library(AER)      # dispersiontest
library(lmtest)   # lrtest, waldtest
library(car)      # vif, influencePlot
library(caret)    # findCorrelation, preProcess
library(corrplot) # Visualizar correlaciones
library(GGally)   # ggpairs para matrices
library(sf)       # Simple Features (shapefiles)
library(sp)       # Spatial points/data
library(boot)     # Validación cruzada
library(rsq)      # R-squared para GLM
library(performance) # check_model, model_performance

setwd("/Users/luisalbertovalencia/Documents/Maestria/Análisis/Datos/N1/Base _final")
datosN<- read.csv("Mosquitos_AGEB.csv")
head(datosN)

# 1. Seleccionar solo variables numéricas para correlación
vars_numericas <- c("PM25_campo","PM10_campo","Temp","RH","POB_TOT","P6A14NA","SBASC","PSDSS","OVSDE","OVSEE","OVSAE","OVPT","OVSREF","OVSINT","OVSCEL","OVHAC","IM_2020","IMN_202","Altitud","CO",
                  "EII","NDVI","NO2","O3","PH","PM10","PM25","SO2","TMAX_02","TMAX_07","TMAX_PROM","TMIN_02","TMIN_07","TMIN_PROM","AREA","PERIMETER","NDWI","LST","NTL","PREC","VPD")


# Matriz de correlación
cor_matrix <- cor(datosN[, vars_numericas], use = "complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)
write.csv(cor_matrix, "matriz_correlaciones.csv")
# Encontrar pares con correlación > 0.7
high_cor_pairs <- which(abs(cor_matrix) > 0.7 & upper.tri(cor_matrix), arr.ind = TRUE)

cat("Pares de variables altamente correlacionadas (|r| > 0.7):\n")
if(nrow(high_cor_pairs) > 0) {
  for(i in 1:nrow(high_cor_pairs)) {
    var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
    cor_value <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
    cat(sprintf("%s - %s: r = %.2f\n", var1, var2, cor_value))
  }
} else {
  cat("No hay pares con |r| > 0.7\n")
}



library(lme4)
library(MuMIn)
####Selección de campo
modelo_campo <- glmer.nb(
  Abundancia_vectores ~ PM25_campo + PM10_campo + Temp + RH +
    (1|Location), data = datosN)
summary(modelo_campo)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_campo)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm

####Selección de sociales
names(datosN)
vars_modelo <- c("POB_TOT","P6A14NA","SBASC","PSDSS",
                 "OVSDE","OVSEE","OVSAE","OVPT","OVSREF","OVSINT","OVSCEL",
                 "OVHAC","IM_2020","IMN_202")

datos_modelo <- datosN[complete.cases(datosN[, vars_modelo]), ]
colSums(is.na(datosN))
modelo_social1 <- glmer.nb(
  Abundancia_vectores ~ OVPT + OVSREF + OVSINT + OVSCEL + OVHAC + IM_2020 + IMN_202 +
    (1|Location), data = datos_modelo)
summary(modelo_social1)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_social1)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm

modelo_social2 <- glmer.nb(
  Abundancia_vectores ~ POB_TOT + P6A14NA + SBASC + PSDSS + OVSDE + OVSEE + OVSAE  +
    (1|Location), data = datos_modelo)
summary(modelo_social2)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_social2)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm

####Selección de ambiental
modelo_ambiental1 <- glmer.nb(
  Abundancia_vectores ~  NDWI + LST + NTL + PREC + VPD +
    (1|Location), data = datosN)
summary(modelo_ambiental1)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_ambiental1)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm

vars <- c("Altitud","EII","NDVI","PH","TMAX_PROM","TMIN_PROM")
datosN[vars] <- scale(datosN[vars])
modelo_ambiental2 <- glmer.nb(
  Abundancia_vectores ~ Altitud + EII + NDVI + PH + TMAX_PROM +
    (1|Location), data = datosN)
summary(modelo_ambiental2)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_ambiental2)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm

####Selección de contaminantes
modelo_contaminantes <- glmer.nb(
  Abundancia_vectores ~ CO + NO2 + O3 + PM10 + PM25 + SO2 +
    (1|Location), data = datosN)
summary(modelo_contaminantes)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_contaminantes)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm
model.sel(top_modelos_lm)
sw(top_modelos_lm)



####Modelo final
datosN$OVHAC[is.na(datosN$OVHAC)] <- 0
datosN$OVSEE[is.na(datosN$OVSEE)] <- 0
modelo_final <- glmer.nb(
  Abundancia_vectores ~ CO + PM25 + NTL + PREC + Altitud + NDVI + PH + TMAX_PROM + OVHAC + OVSEE +
    (1|Location), data = datosN)
summary(modelo_final)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_final)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm
model.sel(top_modelos_lm)
sw(top_modelos_lm)




#####
# 1. Crear dataset eliminando filas con NA en OVHAC o OVSEE
datosN<- read.csv("Mosquitos_AGEB.csv")
datosN$UsoSuelo <- as.factor(datosN$UsoSuelo)
datosN$Temporada <- as.factor(datosN$Temporada)
datosN_noNA <- datosN[!is.na(datosN$OVHAC) & !is.na(datosN$OVSEE), ]
nrow(datosN)
nrow(datosN_noNA)
# Seleccionar variables continuas
vars_cont <- c("CO", "PM25", "NTL", "PREC", "Altitud", 
               "NDVI", "PH", "TMAX_PROM", "OVHAC")

# Escalar
datosN_noNA[vars_cont] <- scale(datosN_noNA[vars_cont])

modelo_final <- glmer.nb(
  Abundancia_vectores ~ CO + PM25 + NTL + PREC + Altitud + NDVI + PH +
     OVHAC + UsoSuelo + Temporada +
    (1 | Location),
  data = datosN_noNA
)
modelo_nulo <- glmer(
  Abundancia_vectores ~ 1 + (1 | Location),
  family = negative.binomial(1.4566),
  data = datosN_noNA
)
AIC(modelo_nulo)
car::vif(modelo_final)
summary(modelo_final)
options(na.action = "na.fail")
todos_modelos_lm <- dredge(modelo_final)
top_modelos_lm <- get.models(todos_modelos_lm, subset = delta < 2)
top_modelos_lm
modelo_avg_lm <- model.avg(top_modelos_lm, fit = TRUE)
summary(modelo_avg_lm)
imp <- sw(modelo_avg_lm)
imp
tabla_importancia <- data.frame(
  Variable = names(imp),
  Importancia = as.numeric(imp)
)

tabla_importancia <- tabla_importancia[order(-tabla_importancia$Importancia), ]
tabla_importancia
best_model <- get.models(todos_modelos_lm, 1)[[1]]
length(top_modelos_lm)
confint(modelo_avg_lm, full = TRUE)   # Incluye coef = 0 si no está en un modelo
confint(modelo_avg_lm, full = FALSE) 



M1 <- glmer.nb(
  Abundancia_vectores ~ CO + NTL + PREC + Altitud + NDVI +
    OVHAC + Temporada +
    (1 | Location),
  data = datosN_noNA
)
car::vif(M1)
# VALIDACIÓN DE SUPUESTOS DEL GLMM (NB) CON DHARMa

library(DHARMa)

# 1. Simulación de residuos
par(mfrow = c(3, 3))
sim_res <- simulateResiduals(
  fittedModel = M1,
  n = 1000
)
par(
  mar = c(4,4,2,1),  # márgenes más limpios
  cex = 1.3          # aumenta tamaño de todo el texto
)
plot(sim_res)
testUniformity(sim_res)
testDispersion(sim_res)
testZeroInflation(sim_res)
plotResiduals(sim_res, datosN_noNA$NTL)
plotResiduals(sim_res, datosN_noNA$NDVI)
plotResiduals(sim_res, datosN_noNA$OVHAC)
plotResiduals(sim_res, datosN_noNA$PREC)
(1 | Location/Temporada)
poly(datosN_noNA$NTL, 2)
M1 <- glmer.nb(
  Abundancia_vectores ~ CO + poly(NTL, 2) + PREC + Altitud + NDVI +
    OVHAC + Temporada +
    (1 | Location),
  data = datosN_noNA
)
# ================================
# LIBRERÍAS
# ================================
library(ggplot2)
library(ggeffects)
library(patchwork)

# ================================
# 1. Dataset sin NAs y variables originales
# ================================
datosN_noNA_original <- datosN[!is.na(datosN$OVHAC) & !is.na(datosN$OVSEE), ]

# ================================
# 2. Ajustar el mejor modelo en variables originales
# ================================
library(lme4)
best_model_original <- glmer.nb(
  Abundancia_vectores ~ NDVI + NTL + OVHAC + PREC + 
     + Temporada +
    (1 | Location),
  data = datosN_noNA_original
)

# ================================
# 3. Librerías para predicción y gráficas
# ================================
library(ggeffects)
library(ggplot2)
library(patchwork)

# ================================
# 4. Predicciones de cada variable
# ================================
pred_co    <- ggpredict(best_model_original, terms = "NDVI", bias_correction = TRUE)
pred_ntl   <- ggpredict(best_model_original, terms = "NTL", bias_correction = TRUE)
pred_ovhac <- ggpredict(best_model_original, terms = "OVHAC", bias_correction = TRUE)
pred_prec  <- ggpredict(best_model_original, terms = "PREC", bias_correction = TRUE)

# ================================
# 5. Gráficas
# ================================
# CO
graf_co <- ggplot() +
  geom_point(data = datosN_noNA_original,
             aes(x = NDVI, y = Abundancia_vectores),
             size = 3, alpha = 0.6, color = "firebrick") +
  geom_line(data = pred_co,
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  geom_ribbon(data = pred_co,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey30") +
  labs(x = "NDVI", y = "Vector abundance") +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, NA)) +  
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 0.8)
  )

# NTL
graf_ntl <- ggplot() +
  geom_point(data = datosN_noNA_original,
             aes(x = NTL, y = Abundancia_vectores),
             size = 3, alpha = 0.6, color = "darkorange") +
  geom_line(data = pred_ntl,
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  geom_ribbon(data = pred_ntl,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey30") +
  labs(x = "Nighttime Light (nW/cm²/sr)", y = "Vector abundance") +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 0.8)
  )

# OVHAC
graf_ovhac <- ggplot() +
  geom_point(data = datosN_noNA_original,
             aes(x = OVHAC, y = Abundancia_vectores),
             size = 3, alpha = 0.6, color = "forestgreen") +
  geom_line(data = pred_ovhac,
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  geom_ribbon(data = pred_ovhac,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey30") +
  labs(x = "OVHAC", y = "Vector abundance") +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 0.8)
  )

# PREC
graf_prec <- ggplot() +
  geom_point(data = datosN_noNA_original,
             aes(x = PREC, y = Abundancia_vectores),
             size = 3, alpha = 0.6, color = "dodgerblue") +
  geom_line(data = pred_prec,
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  geom_ribbon(data = pred_prec,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey30") +
  labs(x = "Precipitation (mm)", y = "Vector abundance") +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  theme(
    text = element_text(family = "Arial", size = 20),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 0.8)
  )

# ================================
# 6. Combinar gráficos y mostrar
# ================================
library(patchwork)
fig_final <- (graf_co + graf_ntl) / (graf_ovhac + graf_prec)
fig_final

# ================================
# 7. Guardar figura
# ================================
ggsave("Vector_abundance_best_model_original.tif",
       plot = fig_final,
       width = 10, height = 10, dpi = 600)
