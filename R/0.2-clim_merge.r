#!!attention!!
# this file cannot be run from the initial data files, one would have to supply
# the separate chelsa bio layers, since they are very large in size

# this file contains the code used to prepare the climate data from CHELSA
# it merges all traditional bioclim layers for each temporal extent

library(raster)

## prepare CHELSA data

# merge all traditional bioclim layers for period 1981-2010 and period 2011-2040 (bio+ has slight extent deviation)
filenames <- list.files(path = "/mnt/f/BscBio_LP23_data/chelsa_bio_layers/1981-2010", full.names = TRUE)
merged <- raster()
for (f in seq(1:length(filenames))) {
  layer <- raster(filenames[f])
  layername <- sub(".*CHELSA_", "", filenames[f])
  layername <- sub("_1981.*", "", layername)
  merged <- stack(l = layer, merged)
  names(merged)[f] <- layername # extract layer name
}
writeRaster(merged, "/mnt/f/BscBio_LP23_data/chelsa_bio_merged/chelsa_bio_1981-2010_global.grd", overwrite = TRUE)
print("merged bio for 1981-2010")

filenames <- list.files(path = "/mnt/f/BscBio_LP23_data/chelsa_bio_layers/2011-2040", full.names = TRUE)
merged <- raster()
for (f in seq(1:length(filenames))) {
  layer <- raster(filenames[f])
  layername <- sub(".*CHELSA_", "", filenames[f])
  layername <- sub("_2011.*", "", layername)
  merged <- stack(l = layer, merged)
  names(merged)[f] <- layername # extract layer name
}
writeRaster(merged, "/mnt/f/BscBio_LP23_data/chelsa_bio_merged/chelsa_bio_2011-2040_global.grd", overwrite = TRUE)
print("merged bio for 2011-2040")
