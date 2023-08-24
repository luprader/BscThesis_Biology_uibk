# this file creates an ensemble of sdm models for each year

# libraries uesd
library(terra)
library(dplyr)

pa_mod = readRDS("R/data/modelling/pa_mod_vars.rds")

data = select(pa_mod, matches("[[:digit:]]"))
data_sc = data.frame(scale(data)) # scale all variables (-mean, /stdev)
scaling = rbind("mean" = colMeans(data), "sd" = apply(data, 2, sd))
data_sc$Presence = factor(pa_mod$Presence)

# build glm with native data
mod_data = data_sc[which(pa_mod$Area == 'as'), ]
mod_glm = glm(Presence ~ ., data = mod_data, family = "binomial")

# create prediction raster from pca and clim
data_r = rast("R/data/modelling/pca_rasters/pca_2002_eu.grd")
clim_r = rast("R/data/cropped_rasters/CHELSA_bio_merged_1981-2010_eu.grd")
data_r = c(subset(clim_r, grep("bio", colnames(data), val = TRUE)), data_r)
data_r_sc = scale(data_r, center = scaling["mean",], scale = scaling["sd", ])

pred_rast = predict(data_r_sc, mod_glm)