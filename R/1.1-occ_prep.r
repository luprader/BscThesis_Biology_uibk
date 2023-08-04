# This file subsets the cleaned occurence points per year and generates
# pseudo-absence points in a x radius around the presences.
# Afterwards, the corresponding bioclim and land cover values
# are extracted to create one merged dataframe for modelling

library(terra)
library(sf)
