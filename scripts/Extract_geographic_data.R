# =============================================================================
# PASO 0: DIAGNÓSTICO DE SISTEMAS DE COORDENADAS
# Verificar CRS de todos los archivos antes de la extracción

library(exactextractr)
library(sf)
library(terra)
setwd("/Users/luisalbertovalencia/Documents/Maestria/Análisis/Datos/N1/Vectores_sel/")
# --- Ajusta estas rutas a donde tengas tus archivos ---
ruta_agebs  <- "poligono_ageb_urbanas_cdmx.shp"
ruta_prec   <- "precip_acum_2024_2025.tif"
ruta_ntl    <- "NTL.tif"
ruta_ndvi   <- "NDVI.tif"

# =============================================================================
# SHAPEFILES
# =============================================================================

agebs <- st_read(ruta_agebs, quiet = TRUE)

cat("============================================================\n")
cat("SHAPEFILE: poligono_ageb_urbanas_cdmx.shp\n")
cat("  CRS      :", st_crs(agebs)$input, "\n")
cat("  EPSG     :", st_crs(agebs)$epsg, "\n")
cat("  Geometría:", unique(st_geometry_type(agebs)), "\n")
cat("  N filas  :", nrow(agebs), "\n")
cat("  Columnas :", paste(names(agebs), collapse = ", "), "\n")
cat("  Bbox     :"); print(st_bbox(agebs))

# =============================================================================
# RASTERS
# =============================================================================

for (info in list(
  list(ruta = ruta_prec, nombre = "precip_acum_2024_2025.tif"),
  list(ruta = ruta_ntl,  nombre = "NTL.tif"),
  list(ruta = ruta_ndvi, nombre = "NDVI.tif")
)) {
  r <- rast(info$ruta)
  cat("\n============================================================\n")
  cat("RASTER:", info$nombre, "\n")
  cat("  CRS        :", crs(r, describe = TRUE)$name, "\n")
  cat("  EPSG       :", crs(r, describe = TRUE)$code, "\n")
  cat("  Resolución :", paste(round(res(r), 4), collapse = " x "), "unidades del CRS\n")
  cat("  Dimensiones:", nrow(r), "filas x", ncol(r), "columnas\n")
  cat("  N capas    :", nlyr(r), "\n")
  cat("  Rango vals :", round(minmax(r)[1], 3), "a", round(minmax(r)[2], 3), "\n")
  cat("  Bbox       :"); print(ext(r))
}

cat("\n============================================================\n")
cat("Diagnóstico completo. Comparte los resultados para continuar.\n")
cat("============================================================\n")
library(terra)

# Paso 1: definir rutas (strings)
ruta_prec <- "precip_acum_2024_2025.tif"
ruta_ntl  <- "NTL.tif"
ruta_ndvi <- "NDVI.tif"

# Paso 2: cargar como objetos raster
r_prec <- rast(ruta_prec)
r_ntl  <- rast(ruta_ntl)
r_ndvi <- rast(ruta_ndvi)

# Paso 3: ahora sí calcular min/max
r_prec <- setMinMax(r_prec)
r_ntl  <- setMinMax(r_ntl)
r_ndvi <- setMinMax(r_ndvi)

# Verificar
cat("PREC:", global(r_prec, "range", na.rm = TRUE)[[1]], "\n")
cat("NTL :", global(r_ntl,  "range", na.rm = TRUE)[[1]], "\n")
cat("NDVI:", global(r_ndvi, "range", na.rm = TRUE)[[1]], "\n")

# =============================================================================
# 2. VERIFICAR QUE LOS EXTENTS SE SOLAPAN
# =============================================================================

bbox_agebs <- st_bbox(agebs)
for (nombre_r in c("PREC", "NTL", "NDVI")) {
  r_tmp <- get(paste0("r_", tolower(nombre_r)))
  ext_r <- ext(r_tmp)
  solapa <- (ext_r$xmax > bbox_agebs["xmin"]) &
    (ext_r$xmin < bbox_agebs["xmax"]) &
    (ext_r$ymax > bbox_agebs["ymin"]) &
    (ext_r$ymin < bbox_agebs["ymax"])
  cat(nombre_r, "- solapa con AGEBs:", solapa, "\n")
}

# =============================================================================
# 3. EXTRACCIÓN DE PREC — valor por centroide (raster muy grueso ~5km)
#    Cada AGEB toma el valor del píxel en su centroide
# =============================================================================

cat("\nExtrayendo PREC por centroide...\n")

centroides <- st_centroid(agebs)
coords_centroides <- st_coordinates(centroides)

prec_vals <- terra::extract(
  r_prec,
  vect(centroides)   # convierte sf a SpatVector para terra
)[, 2]               # columna 2 = valores (columna 1 es ID)

agebs$PREC <- prec_vals
cat("  NAs en PREC:", sum(is.na(agebs$PREC)), "de", nrow(agebs), "AGEBs\n")

# =============================================================================
# 4. EXTRACCIÓN DE NTL — promedio zonal por polígono (~500m resolución)
# =============================================================================

cat("\nExtrayendo NTL por promedio zonal...\n")

ntl_vals <- exact_extract(
  r_ntl,
  agebs,
  fun = "mean",         # promedio ponderado por fracción de píxel cubierta
  progress = FALSE
)

agebs$NTL <- ntl_vals
cat("  NAs en NTL:", sum(is.na(agebs$NTL)), "de", nrow(agebs), "AGEBs\n")

# =============================================================================
# 5. EXTRACCIÓN DE NDVI — promedio zonal por polígono (~10m resolución)
# =============================================================================

cat("\nExtrayendo NDVI por promedio zonal...\n")

# NDVI es el raster más pesado (~5500x5500), puede tardar 1-2 minutos
ndvi_vals <- exact_extract(
  r_ndvi,
  agebs,
  fun = "mean",
  progress = TRUE       # barra de progreso activada
)

agebs$NDVI <- ndvi_vals
cat("  NAs en NDVI:", sum(is.na(agebs$NDVI)), "de", nrow(agebs), "AGEBs\n")



# =============================================================================
# 7. RESUMEN FINAL Y DIAGNÓSTICO
# =============================================================================

cat("\n============================================================\n")
cat("RESUMEN DE COVARIABLES EXTRAÍDAS:\n")
cat("============================================================\n")

agebs_vars <- agebs %>%
  st_drop_geometry() %>%
  dplyr::select(CVEGEO, PREC, NTL, NDVI)

print(summary(agebs_vars))

# Contar AGEBs con todos los datos completos
completas <- agebs_vars %>% filter(complete.cases(.))
cat("\nAGEBs con todas las variables completas:", nrow(completas), "de", nrow(agebs), "\n")

# =============================================================================
# 8. EXPORTAR SHAPEFILE ENRIQUECIDO
# =============================================================================

# Remover columna auxiliar de área antes de guardar
agebs_final <- agebs %>% dplyr::select(-area_total, -area_ovhac)

st_write(
  agebs_final,
  "agebs_cdmx_con_covariables.gpkg",
  driver     = "GPKG",
  delete_dsn = TRUE
)


datos_mosquitos <- read.csv("AGEB_modelo_CDMX.csv")
str(datos_mosquitos)


library(dplyr)

datos_mosquitos <- datos_mosquitos %>%
  rename(CVEGEO = CVEGEO_1)
str(agebs$CVEGEO)
str(datos_mosquitos$CVEGEO)
agebs_final <- agebs %>%
  left_join(datos_mosquitos, by = "CVEGEO")
nrow(agebs)        # antes
nrow(agebs_final)  # después (deben ser iguales)
colSums(is.na(agebs_final[, c("POB_TOTAL", "OVSEE", "OVHAC", "IM_2020")]))
sum(agebs$CVEGEO %in% datos_mosquitos$CVEGEO)
sum(!agebs$CVEGEO %in% datos_mosquitos$CVEGEO)
head(agebs_final %>% st_drop_geometry())
sum(duplicated(datos_mosquitos$CVEGEO))
library(sf)

st_write(agebs_final, "AGEB_modelo_vectores.shp", delete_layer = TRUE)
