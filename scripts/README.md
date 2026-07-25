# Bayesian spatial modeling of mosquito vector abundance in Mexico City

This repository contains the scripts used to analyze the environmental and sociodemographic factors associated with the abundance of medically important mosquito vectors in Mexico City through generalized linear models (GLMs) and Bayesian spatial modeling.

## Repository structure

### 📂 KrigingModels_pollutants

Contains the scripts used to generate kriging interpolation models for the six criteria air pollutants monitored across Mexico City.

Each script corresponds to one pollutant and includes:

- Data preparation
- Variogram fitting
- Ordinary kriging interpolation
- Raster generation for subsequent analyses

These interpolated surfaces were later used as environmental predictors.

---

### 📄 extract_data.R

This script integrates environmental and sociodemographic information at the AGEB (Área Geoestadística Básica) level.

The workflow includes:

- Extraction of raster values (e.g., climatic and environmental variables)
- Extraction of information from shapefiles
- Integration of census-derived sociodemographic variables
- Construction of the final dataset used for statistical analyses

---

### 📄 Abundance_analysis.R

This script performs the statistical analysis to identify variables associated with mosquito vector abundance.

Main steps include:

- Data exploration
- Variable selection
- Generalized Linear Models (GLMs)
- Model selection using Akaike Information Criterion (AIC)

The environmental and sociodemographic indicators selected for the final analysis were:

- NDVI (Normalized Difference Vegetation Index)
- NTL (Night-Time Lights)
- OVHAC (Household overcrowding)
- Precipitation

Given the overdispersed nature of mosquito abundance data, analyses were conducted using a **Negative Binomial** distribution.

---

### 📄 Bayesian_model-vector_abundance.R

This script implements the Bayesian spatial model for mosquito vector abundance.

The workflow includes:

- Spatial mesh construction for Mexico City
- Bayesian spatial modeling
- Prediction at the AGEB level
- Generation of predicted abundance maps
- Classification of predicted abundance into four categories

Bayesian models were evaluated using:

- Deviance Information Criterion (DIC)
- Watanabe-Akaike Information Criterion (WAIC)

The response variable was modeled using a **Negative Binomial** distribution to account for overdispersion.

---

## Workflow

```
Air pollution data
        │
        ▼
Kriging interpolation
        │
        ▼
Environmental rasters
        │
        ▼
extract_data.R
        │
        ▼
Integrated AGEB dataset
        │
        ▼
Abundance_analysis.R
(GLMs + Variable selection)
        │
        ▼
Selected indicators
(NDVI, NTL, OVHAC, Precipitation)
        │
        ▼
Bayesian_model-vector_abundance.R
        │
        ▼
Spatial prediction and abundance maps
```

---

## Statistical methods

- Ordinary Kriging
- Generalized Linear Models (GLMs)
- Negative Binomial regression
- Bayesian spatial modeling (INLA)
- Model selection using AIC
- Model comparison using DIC and WAIC

---

## Study area

Mexico City, Mexico.

Predictions were generated at the AGEB level using environmental, atmospheric, and sociodemographic indicators.

---

## Software

Analyses were performed in **R** using packages for:

- Spatial analysis
- Geostatistics
- Bayesian inference (INLA)
- Raster processing
- Vector data manipulation
- Statistical modeling

---

## Notes

The scripts are organized following the analytical workflow used in the study, from environmental data preprocessing to the final Bayesian spatial prediction of mosquito vector abundance.
