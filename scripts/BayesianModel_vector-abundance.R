# =============================================================================
# MODELO BAYESIANO ESPACIAL (INLA-SPDE) - ABUNDANCIA DE MOSQUITOS VECTORES CDMX
# Variables seleccionadas por GLMM con model averaging (ΔAICc ≤ 2)
# Variables con IC 95% significativo: NTL*, OVHAC*, NDVI*
# Variables de soporte moderado incluidas: PREC, CO (modelo completo)
# =============================================================================

library(sf)
library(dplyr)
library(INLA)
library(ggplot2)
library(viridis)

# =============================================================================
# 1. CARGA DE DATOS
# =============================================================================
setwd("/Users/luisalbertovalencia/Documents/Maestria/Análisis/Datos/N1/Base _final/Bayesian_vectores/")
# Datos de sitios de muestreo con abundancia observada
datos_mosquitos <- read.csv(
  "Mosquitos_AGEB.csv",
  stringsAsFactors = FALSE
)
str(datos_mosquitos)
# Shapefile de AGEBs de CDMX (debe contener las covariables: NTL, OVHAC, NDVI, PREC, CO)
agebs_cdmx <- st_read("AGEB_modelo_vectores.shp")

# Verificar que las columnas necesarias existen en el shapefile de AGEBs
variables_requeridas <- c("NTL", "OVHAC", "NDVI", "PREC")
cat("Columnas en AGEBs:", paste(names(agebs_cdmx), collapse = ", "), "\n")
stopifnot(all(variables_requeridas %in% names(agebs_cdmx)))

# =============================================================================
# 2. PREPARACIÓN ESPACIAL
# =============================================================================

# Convertir puntos de muestreo a objeto espacial (WGS84)
datos_sf <- st_as_sf(
  datos_mosquitos,
  coords = c("Longitud", "Latitud"),
  crs = 4326
)

# Proyectar a UTM zona 14N (EPSG:32614) — adecuado para CDMX, unidades en metros
datos_sf_utm  <- st_transform(datos_sf,  32614)
agebs_utm     <- st_transform(agebs_cdmx, 32614)

# =============================================================================
# 3. EXTRACCIÓN DE COVARIABLES EN LOS SITIOS DE MUESTREO
#    (unir valores de las AGEBs a los puntos de campo)
# =============================================================================

datos_con_vars <- st_join(datos_sf_utm, agebs_utm[, variables_requeridas])

# Variable respuesta: abundancia de vectores
# Verificar el nombre exacto de la columna en tu CSV
# Aquí se asume "abun_vectores" — ajusta si es diferente
variable_respuesta <- "Abundancia_vectores"

variables_requeridas <- c("NTL.y", "OVHAC.y", "NDVI.y", "PREC.y")

variables_modelo <- c(variable_respuesta, variables_requeridas)

datos_modelo <- datos_con_vars %>%
  dplyr::select(all_of(variables_modelo), geometry) %>%
  filter(complete.cases(across(all_of(variables_modelo))))

cat("Observaciones completas para el modelo:", nrow(datos_modelo), "\n")
glimpse(datos_modelo)

# =============================================================================
# 4. ESTANDARIZACIÓN DE COVARIABLES CONTINUAS
#    Importante: INLA es sensible a escalas muy distintas entre variables
#    La estandarización también hace las prioris más interpretables
# =============================================================================

estandarizar <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

datos_modelo <- datos_modelo %>%
  mutate(
    NTL_z   = estandarizar(NTL.y),
    OVHAC_z = estandarizar(OVHAC.y),
    NDVI_z  = estandarizar(NDVI.y),
    PREC_z  = estandarizar(PREC.y)
  )

# Guardar medias y SDs para después re-escalar las predicciones si es necesario
escala_params <- datos_modelo %>%
  st_drop_geometry() %>%
  summarise(
    across(all_of(variables_requeridas),
           list(media = ~mean(.x, na.rm=TRUE), sd = ~sd(.x, na.rm=TRUE)))
  )
mean_y <- mean(datos_modelo$Abundancia_vectores)
var_y  <- var(datos_modelo$Abundancia_vectores)

mean_y
var_y
# =============================================================================
# 5. CONSTRUCCIÓN DE LA MALLA INLA-SPDE
#    La malla debe cubrir toda la extensión de las AGEBs (dominio de CDMX)
#    max.edge: distancias en metros. CDMX ≈ 50 km de diámetro
#    - Interior (max.edge[1]): resolución principal del modelo
#    - Exterior (max.edge[2]): zona buffer para evitar efectos de borde
# =============================================================================

# Boundary del área de estudio (polígono exterior de todas las AGEBs)
limite_cdmx <- st_union(agebs_utm) %>%
  st_convex_hull() %>%
  st_coordinates() %>%
  as.matrix()

# QUEDARSE SOLO CON X Y Y
limite_cdmx <- limite_cdmx[, 1:2]

boundary_inla <- inla.mesh.segment(limite_cdmx)

# Crear malla
# Ajusta max.edge según la densidad de tus puntos:
#   - Si tienes puntos muy separados: aumenta max.edge[1] (ej. 3000-5000)
#   - Si tienes puntos muy densos: puedes bajar a 1000-1500
malla_cdmx <- inla.mesh.2d(
  loc      = st_coordinates(datos_modelo),   # puntos guían la resolución interior
  boundary = boundary_inla,
  max.edge = c(2000, 8000),   # metros: interior=2km, buffer exterior=8km
  cutoff   = 500,             # distancia mínima entre nodos (evita sobre-densificación)
  offset   = c(5000, 10000)   # extensión del buffer interno y externo
)

# Convertir correctamente boundary
boundary_inla <- inla.sp2segment(limite_cdmx)

# Crear malla
malla_cdmx <- inla.mesh.2d(
  loc      = st_coordinates(datos_modelo),
  boundary = boundary_inla,
  max.edge = c(2000, 8000),
  cutoff   = 500,
  offset   = c(5000, 10000),   # ajuste importante
  min.angle = 30
)
malla_cdmx <- inla.mesh.2d(
  loc = st_coordinates(datos_modelo),              # tus puntos
  boundary = boundary_inla,       # polígono
  max.edge = c(1000, 5000),        # km
  cutoff = 0.2,              
  offset = c(5000, 10000),         
  min.angle = 30
)
cat("Número de nodos en la malla:", malla_cdmx$n, "\n")
par(mfrow = c(1, 1))
# Visualizar malla para verificar que cubre bien el área
plot(malla_cdmx, asp = 1, main = "Malla INLA-SPDE para CDMX")
points(st_coordinates(datos_modelo), col = "red", pch = 19, cex = 0.7)
plot(st_geometry(agebs_utm), add = TRUE, border = "blue", lwd = 0.3)

# =============================================================================
# 6. DEFINICIÓN DEL MODELO ESPACIAL (SPDE - Matérn)
#    alpha = 2 → suavidad de Matérn ν = 1 (diferenciable una vez)
#    Es el estándar para fenómenos ecológicos continuos
# =============================================================================

spde <- inla.spde2.pcmatern(
  mesh  = malla_cdmx,
  alpha = 2,
  prior.range = c(5000, 0.05),
  prior.sigma = c(1, 0.01)
)

# Índice para el efecto espacial aleatorio (campo latente gaussiano)
indice_espacial <- inla.spde.make.index(
  "espacial",
  n.spde = spde$n.spde
)

# =============================================================================
# 7. MATRICES DE PROYECCIÓN A
#    Proyectan el campo definido en la malla hacia las ubicaciones observadas
# =============================================================================

A_obs <- inla.spde.make.A(
  malla_cdmx,
  loc = st_coordinates(datos_modelo)
)

# =============================================================================
# 8. STACK DE DATOS PARA OBSERVACIONES
# =============================================================================

stack_obs <- inla.stack(
  data = list(y = datos_modelo[[variable_respuesta]]),
  A    = list(A_obs, 1),        # A_obs para el campo espacial, 1 para efectos fijos
  effects = list(
    espacial = indice_espacial,
    data.frame(
      intercepto = 1,           # intercepto explícito (recomendado con SPDE)
      NTL_z      = datos_modelo$NTL_z,
      OVHAC_z    = datos_modelo$OVHAC_z,
      NDVI_z     = datos_modelo$NDVI_z,
      PREC_z     = datos_modelo$PREC_z
    )
  ),
  tag = "observaciones"
)

# =============================================================================
# 9. FÓRMULA DEL MODELO
#    Variables principales (IC 95% significativo): NTL_z, OVHAC_z, NDVI_z
#    Variables de soporte moderado: PREC_z, CO_z
#    -1 en fórmula porque definimos intercepto explícitamente en el stack
# =============================================================================

formula_modelo <- y ~
  -1 +                           # sin intercepto automático (está en el stack)
  intercepto +                   # intercepto explícito
  NTL_z   +                      # Night-Time Lights (+, IC* significativo)
  OVHAC_z +                      # Ocupación/hábitat (+, IC* significativo)
  NDVI_z  +                      # Vegetación (-, IC* significativo)
  PREC_z  +                      # Precipitación (-, IC marginal)
  f(espacial, model = spde)      # Campo espacial latente (Matérn)

# =============================================================================
# 10. AJUSTE DEL MODELO INLA
#     Familia: Poisson (como en tu GLMM original)
#     Si hay sobredispersión puedes cambiar a "nbinomial"
# =============================================================================

modelo_inla <- inla(
  formula_modelo,
  family = "poisson",           # consistente con tu GLMM de selección de variables
  data   = inla.stack.data(stack_obs),
  control.predictor = list(
    A       = inla.stack.A(stack_obs),
    compute = TRUE,
    link    = 1                 # link canónico (log para Poisson)
  ),
  control.compute = list(
    dic    = TRUE,              # DIC: criterio de información deviance
    waic   = TRUE,              # WAIC: alternativa robusta al DIC
    cpo    = TRUE,              # CPO: leave-one-out para validación
    config = TRUE               # guarda configuración para muestras posteriores
  ),
  control.fixed = list(
    # Prioris débilmente informativas para efectos fijos (escala log)
    # Con variables estandarizadas, efectos típicos están entre -2 y 2
    mean = list(default = 0),
    prec = list(default = 0.1)  # sd ≈ 3.16 en escala del coeficiente
  ),
  verbose = FALSE               # cambia a TRUE si necesitas diagnóstico
)

# =============================================================================
# 11. REVISIÓN DE RESULTADOS
# =============================================================================

cat("\n============================================================\n")
cat("RESUMEN DEL MODELO INLA\n")
cat("============================================================\n")
summary(modelo_inla)

cat("\nEFECTOS FIJOS (escala log):\n")
print(round(modelo_inla$summary.fixed, 4))

cat("\nHIPERPARÁMETROS DEL CAMPO ESPACIAL:\n")
print(round(modelo_inla$summary.hyperpar, 4))

cat("\nCRITERIOS DE AJUSTE:\n")
cat("  DIC :", round(modelo_inla$dic$dic,  2), "\n")
cat("  WAIC:", round(modelo_inla$waic$waic, 2), "\n")

# Bondad de ajuste con CPO (log-score)
# Valores más cercanos a 0 indican mejor ajuste predictivo
cpo_vals <- modelo_inla$cpo$cpo
cat("  Log-score CPO (media):", round(mean(log(cpo_vals), na.rm=TRUE), 4), "\n")

# =============================================================================
# 12. GRÁFICO DE EFECTOS FIJOS
# =============================================================================

efectos_df <- as.data.frame(modelo_inla$summary.fixed)
efectos_df$variable <- rownames(efectos_df)
efectos_df <- efectos_df[efectos_df$variable != "intercepto", ]

# Etiquetas legibles
etiquetas <- c(
  NTL_z   = "NTL (luz nocturna)",
  OVHAC_z = "OVHAC (hábitat)",
  NDVI_z  = "NDVI (vegetación)",
  PREC_z  = "Precipitación"
)
efectos_df$etiqueta <- etiquetas[efectos_df$variable]

p_efectos <- ggplot(efectos_df, aes(x = reorder(etiqueta, mean), y = mean)) +
  geom_point(size = 3, color = "#2c7bb6") +
  geom_errorbar(
    aes(ymin = `0.025quant`, ymax = `0.975quant`),
    width = 0.25, linewidth = 0.8, color = "#2c7bb6"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title    = "Efectos fijos del modelo INLA-SPDE",
    subtitle = "Variables estandarizadas — escala log (Poisson)",
    x        = NULL,
    y        = "Coeficiente posterior (media ± IC 95%)"
  )
print(p_efectos)


# residuos tipo Pearson aproximados
residuos <- (datos_modelo$Abundancia_vectores - 
               modelo_inla$summary.fitted.values$mean[1:nrow(datos_modelo)]) /
  sqrt(modelo_inla$summary.fitted.values$mean[1:nrow(datos_modelo)])

var(residuos)
# =============================================================================
# 13. PREDICCIÓN ESPACIAL SOBRE TODAS LAS AGEBs
# =============================================================================

# Centroides de las AGEBs como ubicaciones de predicción
centroides_agebs <- st_centroid(agebs_utm)

# Estandarizar covariables del grid usando los parámetros de los datos observados
# (muy importante: usar la misma escala que en el entrenamiento)
estandarizar_con_params <- function(x, variable) {
  media <- escala_params[[paste0(variable, "_media")]]
  desv  <- escala_params[[paste0(variable, "_sd")]]
  (x - media) / desv
}



########################################################
#### Pruebas de modelos
#######################################################
datos_modelo$ID_sitio <- as.factor(datos_con_vars$ID_sitio[match(
  st_coordinates(datos_modelo)[,1],
  st_coordinates(datos_con_vars)[,1]
)])
datos_modelo$ID_sitio <- as.factor(datos_modelo$ID_sitio)
stack_obs <- inla.stack(
  data = list(y = datos_modelo[[variable_respuesta]]),
  A    = list(A_obs, 1),
  effects = list(
    espacial = indice_espacial,
    data.frame(
      intercepto = 1,
      NTL_z      = datos_modelo$NTL_z,
      OVHAC_z    = datos_modelo$OVHAC_z,
      NDVI_z     = datos_modelo$NDVI_z,
      PREC_z     = datos_modelo$PREC_z,
      ID_sitio   = datos_modelo$ID_sitio   # 🔥 AQUÍ
    )
  ),
  tag = "observaciones"
)
formula_interaccion <- y ~
  -1 +
  intercepto +
  NTL_z +
  OVHAC_z +
  NDVI_z +
  PREC_z +
  NTL_z:NDVI_z +        # 🔥 interacción
  f(espacial, model = spde)
modelo_inla_int <- inla(
  formula_interaccion,
  family = "poisson",
  data   = inla.stack.data(stack_obs),
  control.predictor = list(
    A = inla.stack.A(stack_obs),
    compute = TRUE
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.fixed = list(
    mean = list(default = 0),
    prec = list(default = 0.1)
  )
)
modelo_nb <- inla(
  formula_modelo,
  family = "nbinomial",
  data   = inla.stack.data(stack_obs),
  control.predictor = list(
    A = inla.stack.A(stack_obs),
    compute = TRUE
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
data.frame(
  Modelo = c("Poisson_SPDE", "NB_SPDE"),
  DIC  = c(modelo_inla$dic$dic, modelo_nb$dic$dic),
  WAIC = c(modelo_inla$waic$waic, modelo_nb$waic$waic)
)
formula_nb_sitio <- y ~
  -1 +
  intercepto +
  NTL_z +
  OVHAC_z +
  NDVI_z +
  PREC_z +
  f(espacial, model = spde) +
  f(ID_sitio, model = "iid")   # 🔥 efecto aleatorio de sitio
modelo_nb_sitio <- inla(
  formula_nb_sitio,
  family = "nbinomial",
  data   = inla.stack.data(stack_obs),
  control.predictor = list(
    A = inla.stack.A(stack_obs),
    compute = TRUE
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  )
)
comparacion_final <- data.frame(
  Modelo = c("NB_SPDE", "NB_SPDE + sitio"),
  DIC = c(
    modelo_nb$dic$dic,
    modelo_nb_sitio$dic$dic
  ),
  WAIC = c(
    modelo_nb$waic$waic,
    modelo_nb_sitio$waic$waic
  )
)

print(comparacion_final)
#### Modelo sin campo espacial
formula_sin_espacial <- Abundancia_vectores ~
  NTL_z +
  OVHAC_z +
  NDVI_z +
  PREC_z
modelo_inla_no_spatial <- inla(
  formula_sin_espacial,
  family = "poisson",
  data   = datos_modelo,   # 🔥 aquí NO uses stack
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.fixed = list(
    mean = list(default = 0),
    prec = list(default = 0.1)
  )
)
comparacion <- data.frame(
  Modelo = c("Con espacial", "Con espacial + interacción", "Sin espacial"),
  DIC = c(
    modelo_inla$dic$dic,
    modelo_inla_int$dic$dic,
    modelo_inla_no_spatial$dic$dic
  ),
  WAIC = c(
    modelo_inla$waic$waic,
    modelo_inla_int$waic$waic,
    modelo_inla_no_spatial$waic$waic
  )
)

print(comparacion)

####################################################
#Pruebas de priors
####################################################
# ============================================================
# SENSIBILIDAD DE PRIORS (NB + SPDE, SIN SITIO)
# ============================================================

library(INLA)
library(dplyr)

# Valores de rango a probar (en metros)
rangos <- c(2000, 5000, 10000)

# Lista para guardar modelos
resultados_priors <- list()

# Loop
for (r in rangos) {
  
  cat("\n==============================\n")
  cat("Probando prior.range =", r, "m\n")
  cat("==============================\n")
  
  # 1. Definir SPDE con nuevo prior
  spde_tmp <- inla.spde2.pcmatern(
    mesh  = malla_cdmx,
    alpha = 2,
    prior.range = c(r, 0.05),
    prior.sigma = c(1, 0.01)
  )
  
  # 2. Índice espacial
  indice_tmp <- inla.spde.make.index(
    "espacial",
    n.spde = spde_tmp$n.spde
  )
  
  # 3. Stack actualizado
  stack_tmp <- inla.stack(
    data = list(y = datos_modelo[[variable_respuesta]]),
    A    = list(A_obs, 1),
    effects = list(
      espacial = indice_tmp,
      data.frame(
        intercepto = 1,
        NTL_z      = datos_modelo$NTL_z,
        OVHAC_z    = datos_modelo$OVHAC_z,
        NDVI_z     = datos_modelo$NDVI_z,
        PREC_z     = datos_modelo$PREC_z
      )
    ),
    tag = paste0("prior_", r)
  )
  
  # 4. Modelo NB + SPDE
  modelo_tmp <- inla(
    formula_modelo,
    family = "nbinomial",
    data   = inla.stack.data(stack_tmp),
    control.predictor = list(
      A = inla.stack.A(stack_tmp),
      compute = TRUE
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE
    ),
    verbose = FALSE
  )
  
  # Guardar modelo
  resultados_priors[[as.character(r)]] <- modelo_tmp
}


############################################
##### Pruebas de INLA AJUSTES DE LAPLACE
#############################################
# ============================================================
# COMPARACIÓN DE ESTRATEGIAS INLA
# ============================================================
estrategias <- c("gaussian", "simplified.laplace", "laplace")
resultados_estrategias <- list()
for (s in estrategias) {
  
  cat("\n==============================\n")
  cat("Probando estrategia:", s, "\n")
  cat("==============================\n")
  
  modelo_tmp <- inla(
    formula_modelo,
    family = "nbinomial",
    data   = inla.stack.data(stack_obs),
    control.predictor = list(
      A = inla.stack.A(stack_obs),
      compute = TRUE
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE
    ),
    control.inla = list(
      strategy = s
    ),
    verbose = FALSE
  )
  
  resultados_estrategias[[s]] <- modelo_tmp
}
comparacion_estrategias <- data.frame(
  Estrategia = estrategias,
  DIC  = sapply(resultados_estrategias, function(m) m$dic$dic),
  WAIC = sapply(resultados_estrategias, function(m) m$waic$waic)
)
print(comparacion_estrategias)


# ============================================================
# COMPARACIÓN DE MODELOS
# ============================================================

comparacion_priors <- data.frame(
  Range_m = rangos,
  DIC  = sapply(resultados_priors, function(m) m$dic$dic),
  WAIC = sapply(resultados_priors, function(m) m$waic$waic)
)

print(comparacion_priors)

# ============================================================
# EXTRA: revisar parámetros espaciales
# ============================================================

for (r in names(resultados_priors)) {
  cat("\n--- Prior range:", r, "m ---\n")
  print(resultados_priors[[r]]$summary.hyperpar)
}





###########################################
####Modelo con dependencia espacial es el mejor, es el que se ocupa a continuación
####################################################
#Modelo de predicción espacial
agebs_pred <- agebs_utm %>%
  mutate(
    NTL_z   = estandarizar_con_params(NTL,   "NTL.y"),
    OVHAC_z = estandarizar_con_params(OVHAC, "OVHAC.y"),
    NDVI_z  = estandarizar_con_params(NDVI,  "NDVI.y"),
    PREC_z  = estandarizar_con_params(PREC,  "PREC.y")
  ) %>%
  filter(complete.cases(across(all_of(
    c("NTL_z", "OVHAC_z", "NDVI_z", "PREC_z")
  ))))

# Coordenadas de predicción (centroides de AGEBs)
coords_pred <- st_coordinates(st_centroid(agebs_pred))

# Matriz de proyección para predicción
A_pred <- inla.spde.make.A(malla_cdmx, loc = coords_pred)

# Stack de predicción
stack_pred <- inla.stack(
  data = list(y = NA),           # NA indica que queremos predecir aquí
  A    = list(A_pred, 1),
  effects = list(
    espacial = indice_espacial,
    data.frame(
      intercepto = 1,
      NTL_z      = agebs_pred$NTL_z,
      OVHAC_z    = agebs_pred$OVHAC_z,
      NDVI_z     = agebs_pred$NDVI_z,
      PREC_z     = agebs_pred$PREC_z
    )
  ),
  tag = "prediccion"
)

# Combinar stacks (observaciones + predicción)
stack_completo <- inla.stack(stack_obs, stack_pred)

# Re-ajustar modelo con stack completo
# control.mode reutiliza el modo posterior del modelo previo (más eficiente)
modelo_pred <- inla(
  formula_modelo,
  family = "nbinomial",   # 🔥 CAMBIO 1
  data   = inla.stack.data(stack_completo),
  control.predictor = list(
    A       = inla.stack.A(stack_completo),
    compute = TRUE,
    link    = 1
  ),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE),
  control.mode    = list(
    theta   = modelo_nb$mode$theta,  # 🔥 CAMBIO 2
    restart = FALSE
  ),
  control.inla = list(strategy = "laplace"),
  verbose = FALSE
)

# =============================================================================
# 14. EXTRACCIÓN Y VISUALIZACIÓN DE PREDICCIONES
# =============================================================================

# Índices de las filas del stack que corresponden a "prediccion"
idx_pred <- inla.stack.index(stack_completo, tag = "prediccion")$data

# Media posterior de la abundancia predicha (en escala de la respuesta, no log)
agebs_pred$abundancia_media  <- modelo_pred$summary.fitted.values$mean[idx_pred]
agebs_pred$abundancia_sd     <- modelo_pred$summary.fitted.values$sd[idx_pred]
agebs_pred$abundancia_q025   <- modelo_pred$summary.fitted.values$`0.025quant`[idx_pred]
agebs_pred$abundancia_q975   <- modelo_pred$summary.fitted.values$`0.975quant`[idx_pred]

# Amplitud del intervalo de credibilidad (incertidumbre)
agebs_pred$ic_amplitud <- agebs_pred$abundancia_q975 - agebs_pred$abundancia_q025

# --- Mapa de abundancia predicha ---
p_mapa_media <- ggplot() +
  geom_sf(data = agebs_pred, aes(fill = abundancia_media), color = NA) +
  geom_sf(
    data   = datos_modelo,
    aes(size = .data[[variable_respuesta]]),
    color  = "white", shape = 21, fill = "red", alpha = 0.8
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Abundancia\npredicha\n(media posterior)"
  ) +
  scale_size_continuous(name = "Abundancia\nobservada", range = c(1, 5)) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Predicción de abundancia de mosquitos vectores — CDMX",
    subtitle = "Modelo Bayesiano INLA-SPDE · NTL + OVHAC + NDVI + PREC",
    caption  = "Puntos: sitios de muestreo  |  Polígonos: AGEBs"
  )
print(p_mapa_media)

# --- Mapa de incertidumbre (amplitud IC 95%) ---
p_mapa_ic <- ggplot() +
  geom_sf(data = agebs_pred, aes(fill = ic_amplitud), color = NA) +
  scale_fill_viridis_c(
    option = "mako", direction = -1,
    name   = "Amplitud\nIC 95%\n(incertidumbre)"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Incertidumbre de la predicción — CDMX",
    subtitle = "Amplitud del intervalo de credibilidad al 95%"
  )
print(p_mapa_ic)

# =============================================================================
# 15. EXPORTAR RESULTADOS
# =============================================================================

# Guardar shapefile con predicciones por AGEB
st_write(
  agebs_pred,
  "predicciones_mosquitos_AGEB_CDMX.gpkg",
  driver = "GPKG",
  delete_dsn = TRUE
)

# Guardar tabla de efectos fijos
write.csv(
  modelo_inla$summary.fixed,
  "efectos_fijos_inla_mosquitos.csv",
  row.names = TRUE
)
write.csv(
  data.frame(
    DIC  = modelo_inla$dic$dic,
    WAIC = modelo_inla$waic$waic,
    CPO  = mean(log(modelo_inla$cpo$cpo), na.rm = TRUE)
  ),
  "metricas_modelo.csv",
  row.names = FALSE
)
write.csv(
  modelo_inla$summary.hyperpar,
  "hiperparametros_spde.csv",
  row.names = TRUE
)

cat("\nArchivos exportados correctamente.\n")
cat("Proceso finalizado.\n")
