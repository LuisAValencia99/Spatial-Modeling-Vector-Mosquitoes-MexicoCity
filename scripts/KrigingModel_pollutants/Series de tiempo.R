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
setwd("/Users/luisalbertovalencia/Documents/Maestria/Calidad de aire/Jalisco/Final/")
Jaliso<-read_csv("newdatos2.csv")
sum(Jaliso_ancho$CO <= 0)
Jalisco_filtrado <- Jaliso_ancho[!is.na(Jaliso_ancho$Fecha),]

CO_NL_filtrado <- Jalisco_filtrado %>%
  select(Fecha, Hora, CO, Hoja) %>%
  filter(!is.na(PM10) & !PM10 %in% c("IF", "IO", "IR", "ND", "VE", "SE", "NaN"))
print(CO_NL_filtrado)

sum(CO_NL_filtrado$PM10 == "IF")
sum(CO_NL_filtrado$PM10 <0)
CO_NL_filtrado <- CO_NL_filtrado %>%
  filter(PM10 > 0)
print(CO_NL_filtrado)

# Encontrar la estación con más datos
estacion_maxima <- CO_NL_filtrado %>%
  group_by(Hoja) %>%
  summarise(horas_con_datos = n()) %>%
  arrange(desc(horas_con_datos)) %>%
  slice(1)  
max_datos <- estacion_maxima$horas_con_datos
porcentaje_comparativo <- CO_NL_filtrado %>%
  group_by(Hoja) %>%
  summarise(horas_con_datos = n()) %>%
  mutate(porcentaje = (horas_con_datos / max_datos) * 100) %>%
  arrange(desc(porcentaje))  
print(porcentaje_comparativo)

#Eliminar estaciones con menos del 70% de datos
CONL<-CO_NL_filtrado[!CO_NL_filtrado$Hoja %in% c("PIN", "ATM"),]

#Revisar que los valores no superen los limites del equipo 


str(CONL$PM10)

CONL$PM10 <- as.numeric(CONL$PM10)
CONL <- CONL %>% filter(!is.na(PM10))
head(CONL)
CONL <- CONL %>%
  filter(PM10 < 0.5)
sum(CONL$PM10 > 100)

# Promedios diarios
CONL$Fecha <- as.Date(CONL$Fecha, format = "%d/%m/%y")
promedio_diario <- CONL %>%
  group_by(Fecha) %>%
  summarise(promedio_valor = mean(PM25, na.rm = TRUE))
print(promedio_diario)
plot(promedio_diario)
write.csv(promedio_diario, file = "PM25_1.csv", row.names = FALSE)


#Para PM
CONL$Fecha <- as.Date(CONL$Fecha, format = "%d/%m/%y")
promedio_diario <- CONL %>%
  group_by(Fecha) %>%
  summarise(
    promedio_valor = mean(PM10, na.rm = TRUE),
    percentil_90 = quantile(PM10, 0.9, na.rm = TRUE) 
  )
print(promedio_diario)
plot(promedio_diario$Fecha, promedio_diario$promedio_valor, type = "l", col = "blue", 
     xlab = "Fecha", ylab = "Promedio Valor", main = "Promedio Diario de Valor")
plot(promedio_diario$Fecha, promedio_diario$percentil_90, type = "l", col = "blue", 
     xlab = "Fecha", ylab = "Percentil; Valor", main = "Percentil 90 Diario de Valor")
write.csv(promedio_diario, file = "PM10_1.csv", row.names = FALSE)
