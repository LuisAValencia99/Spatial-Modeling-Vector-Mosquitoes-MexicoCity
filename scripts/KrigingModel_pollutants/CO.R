library(sf)
library(sp)
library(MASS)
library(scatterplot3d)
library(geoR)
library(fields)
library(maps)
library(grid)
library(spatstat)
library(data.table)
library(chron)
library(zoo)
library(reshape2)
library(DataCombine)
library(dplyr)
library(openair)
library(ggplot2)
library(mgcv)
library(tibble)
library(anytime)
library(lubridate)
library(tidyr)
library(tidyverse)

setwd("/Users/luisalbertovalencia/Documents/Maestria/Mapas/Contaminantes/CO/")
COOR<-read.csv("Coordenadas.csv")
CO1<-read.csv("CO_TOTAL.csv")

#### Unificar fechas
unique(CO1[["FECHA"]])

# Corregir las fechas al formato estándar YYYY-MM-DD
CO1$FECHA <- ifelse(
  grepl("/", CO1$FECHA), 
  as.character(dmy(CO1$FECHA)), 
  as.character(ymd(CO1$FECHA))  
)
unique(CO1$FECHA)



CDMX<-st_read("09ent.shp")
coords <- st_coordinates(CDMX_geo)
st_crs(CDMX)
CDMX_geo <- st_transform(CDMX, crs = 4326)
st_crs(CDMX_geo)  
coords <- st_coordinates(CDMX_geo)
head(st_coordinates(CDMX_geo))


#Checamos las coordenadas geográficas
print(CDMX_geo)
plot(CDMX_geo$geometry, col = "transparent", border = "black", lwd = 1)
puntos_sf <- st_as_sf(COOR, coords = c("Longitud", "Latitud"), crs = 4326)
st_crs(puntos_sf)  
ggplot() +
  geom_sf(data = CDMX_geo, fill = "lightgray", color = "black") +  
  geom_sf(data = puntos_sf, color = "red", size = 3) +  
  ggtitle("Estaciones de RAMA") +
  theme_minimal()


##### Promedios diarios
DailyPollutantAverage <- function(data, date_column = "FECHA", hour_column = "HORA") {
  data[[date_column]] <- as.Date(data[[date_column]])
  data <- data %>%
    mutate(across(where(is.numeric), ~ replace(., . == -99, NA)))
  daily_avg <- data %>%
    group_by_at(date_column) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
  
  return(daily_avg)
}
CO_Daverage<-DailyPollutantAverage(CO1)         

#### Filtro de años
CO_Daverage$FECHA <- as.Date(CO_Daverage$FECHA)
CO_fechas <- CO_Daverage %>%
  filter(format(FECHA, "%Y") %in% c("2017", "2018", "2019", "2022", "2023", "2024"))
#####Promedios totales
OverallStats <- function(data, date_column = "FECHA", hour_column = "HORA") {
  data <- data %>%
    dplyr::select(-dplyr::all_of(c(date_column, hour_column)))
  
  data <- data %>%
    dplyr::mutate(across(where(is.numeric), ~ replace(., . == -99, NA)))
  
  stats <- data %>%
    dplyr::summarise(across(where(is.numeric), 
                            list(mean = ~ mean(., na.rm = TRUE), 
                                 sd = ~ sd(., na.rm = TRUE)), 
                            .names = "{.col}_{.fn}"))
  
  stats <- tidyr::pivot_longer(stats, everything(), 
                               names_to = c("Variable", ".value"), 
                               names_sep = "_")
  return(stats)
}
CO_stats<-OverallStats(CO_fechas)

####Verificar el porcentaje de datos faltantes con el filtrado
Missing1 <- function(x, missing_value = NaN) { 
  x <- x %>%
    mutate(across(where(is.character), as.numeric))  
  missing_counts <- sapply(x, function(col) sum(is.nan(col), na.rm = TRUE))  
  total_rows <- nrow(x)
  missing_percentages <- (missing_counts / total_rows) * 100
  table_final <- data.frame(
    Variable = names(missing_counts),
    Missing_Values = missing_counts,
    Missing_Percentage = round(missing_percentages, 2)  
  )
  table_final <- table_final[!table_final$Variable %in% c("FECHA", "HORA"), ]
  
  return(table_final)
}
CO_F<-Missing1(CO_fechas)

#####
CO_F1 <- CO_fechas[, !colnames(CO_fechas) %in% c("HORA",
                                                 "AJU", "ARA", "AZC", "CES", "COY", "GAM", "IMP", "LAG",
                                                 "PLA", "SUR", "TAC", "TAX", "TEC", "TPN", "VAL", "SJA",
                                                 "CUT", "LLA", "FAR", "SFE", "SAC", "CAM", "HGM", "TAH",
                                                 "INN", "XAL", "ACO", "TLI"
)]



# Reorganizar la tabla y limpiar los datos
data_reorganized <- CO_F1 %>% 
  pivot_longer(
    cols = -c(FECHA),  
    names_to = "Sitio",     
    values_to = "Valor"     
  ) %>%
  mutate(Valor = abs(Valor)) %>%
  filter(!is.na(Valor) & !is.nan(Valor) & Valor != -99 & Valor != 0 & Valor <= 10)
##### GEORREFERENCIAR
CO_coor<- data_reorganized %>% 
  left_join(COOR, by = "Sitio")

#Chcamos de nuevo las coordenadas geográficas
plot(CO_coor$Longitud,CO_coor$Latitud)
plot(CDMX$geometry, col = "transparent", border = "black", lwd = 1)


####Análisis

#Grafica de dispersión para todas las variables.
pairs(CO_coor[, c("Latitud", "Longitud", "Valor")])

#Otro tipo de gráfica con curva suavizada sobrepuesta.  Nos ayudan a ver tendencia y 
# a checar homogeneidad de varianza
CO_coor_clean <- na.omit(CO_coor)
plot(CO_coor_clean$Longitud, CO_coor_clean$Valor)
lines(lowess(CO_coor_clean$Longitud, CO_coor_clean$Valor), col=2)


#Histogramas.  Nos permiten checar el número de modas
######TOTAL
hist(CO_coor$Valor,nclass=30)
hist(CO_coor_clean$Valor,nclass=30)

CO_coor$Mes <- format(CO_coor$FECHA, "%m")
CO_coor_clean$Mes <- format(CO_coor_clean$FECHA, "%m")
mes <- c("Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", 
         "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre")
par(mfrow = c(3, 2))
for (i in 4:8) {  
  hist(
    CO_coor$Valor[CO_coor$Mes == sprintf("%02d", i)],  
    nclass = 30,  
    xlab = "Concentración CO (ppm)",
    ylab = "Frecuencia",  
    main = mes[i],  
    col = "lightblue",  
    border = "black"  
  )
}


#Boxplots
par(mfrow=c(1,1))

plot(CO_coor$Mes,CO_coor$Valor)

boxplot(split(CO_coor_clean$Valor,as.factor(CO_coor_clean$Sitio)))
boxplot(split(CO_coor_clean$Valor,as.factor(CO_coor_clean$Mes)))
####Análisis espacial
# coordenadas y el otro elemento son los valores de las variables de interes
CO.geo<-as.geodata(cbind(CO_coor$Longitud,CO_coor$Latitud,CO_coor$Valor))
plot(CO.geo)
coordinates(CO.geo) <- ~Longitud + Latitud
proj4string(CO.geo) <- CRS("+proj=longlat +datum=WGS84")
CO.geo_utm <- spTransform(CO.geo, CRS("+proj=utm +zone=14 +datum=WGS84"))
print(proj4string(CO.geo_utm))
#Tendencia lineal
CO.lm<-lm(Valor~Longitud+Latitud,data=CO_coor)
summary(CO.lm)
CO.lm<-lm(Valor~Longitud+Latitud+I(Longitud^2)+I(Latitud^2)+I(Longitud*Latitud),data=CO_coor)
summary(CO.lm)

anova(CO.lm)
#removiendo la tendencia con polinomio de orden 1
variog.CO1<- variog(CO.geo,max.dist=10,trend="1st")
plot(variog.CO1)

#removiendo la tendencia con polinomio de orden 2
variog.CO2<- variog(CO.geo,max.dist=10,trend="2nd")
plot(variog.CO2)

# Si hay problema de coordenadas duplicadas

dup.coords(CO.geo)

CO.geo<-jitterDupCoords(CO.geo, max=.01)
variog.CO1<- variog(CO.geo,max.dist=10,trend="1st")
plot(variog.CO1)

# Método por mínimos cuadrados no lineales

vario.CO.mod1<-variofit(variog.CO1,ini.cov.pars=c(2.14,1.04),cov.model="exponential")
lines(vario.CO.mod1)

# Metodo a ojímetro
vario.CO.mod2<-eyefit(variog.CO1)

#### an introduction of statistical learning with applications in R
# Evaluamos el variograma ajustado mediante validación cruzada
valcruz.CO<-xvalid(CO.geo,model=vario.CO.mod1)
plot(valcruz.CO)

################################
#Prediccion Espacial
# Primero obtenemos el mapa de puntos de muestreo

plot(CO.geo$coords,xlab="Longitud",ylab="Latitud")
summary(CO.geo$coords)

#creamos una red de 10 mil puntos donde vamos a predecir
plot(CO.geo$coords)
mired<-expand.grid(long=seq(-99.39,-98.88,length=100),lat=seq(19.00,19.70,length=110))
points(mired,col=2,pch=19,cex=.4)

############# Usar el variograma que se haya elegido

coords <- st_coordinates(CDMX_geo)
plot(coords[,1], coords[,2], pch = ".")

####
CO.krig1<-krige.conv(CO.geo,loc=mired,krige=krige.control(obj.model=vario.CO.mod1))

image(CO.krig1,add=T)
plot(CO.geo$coords,xlab="Longitud",ylab="Latitud")
image(CO.krig1, main = "Kriging CO", xlab = "Longitud", ylab = "Latitud", col = heat.colors(100))
contour(CO.krig1, add = TRUE, col = "blue", lwd = 1)
plot(st_geometry(CDMX_geo), add = TRUE, border = "black", lwd = 1.5)


library(raster)

krig_df <- data.frame(
  x = mired[,1],
  y = mired[,2],
  z = CO.krig1$predict
)
krig_raster <- rasterFromXYZ(krig_df)
crs(krig_raster) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
crs(krig_raster) <- st_crs(CDMX_geo)$wkt

writeRaster(krig_raster, filename = "kriging_CO.tif", format = "GTiff", overwrite = TRUE)
crs(krig_raster)



