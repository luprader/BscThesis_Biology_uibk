#!!attention!!
# this file cannot be run from the initial data files, one would have to supply
# the separate chelsa bio layers

# this file contains the code used to prepare the climate data from CHELSA
# it merges all traditional bioclim layers for each temporal extent,
# crops the raster to the spatial extent of europe and the native range and
# reprojects the raster to match the land cover layers

library(terra)
# file paths should be mnt to external (raw files)

t_start = Sys.time()
## merge all traditional bioclim layers for period 1981-2010
filenames <- list.files(path = "R/data/CHELSA_bio_layers/1981-2010", full.names = TRUE)

merged <- rast() # initialise raster
for (f in filenames) {
    layer <- rast(f) # load layer in question
    name <- strsplit(f, split = "_") # extract variable name
    names(layer) <- name[[1]][4] # rename layer
    merged <- c(merged, layer)
}
## crop and reproject with cropped land cover as reference
# load references
lc_eu <- rast("R/data/cropped_rasters/Cop_LC_1992_eu.grd")
lc_as <- rast("R/data/cropped_rasters/Cop_LC_1992_as.grd")

# crop and reproject raster
clim_eu <- crop(merged, ext(lc_eu), snap = "out")
clim_eu <- resample(clim_eu, lc_eu, method = "bilinear")
clim_as <- crop(merged, ext(lc_as), snap = "out")
clim_as <- resample(clim_as, lc_as, method = "bilinear")

# save croped rasters
fname <- "R/data/cropped_rasters/CHELSA_bio_merged_1981-2010_eu.grd"
writeRaster(clim_eu, filename = fname, overwrite = TRUE)
unlink(paste(fname, ".aux.xml", sep = "")) # remove unnecessary .xml file
fname <- "R/data/cropped_rasters/CHELSA_bio_merged_1981-2010_as.grd"
writeRaster(clim_as, filename = fname, overwrite = TRUE)
unlink(paste(fname, ".aux.xml", sep = "")) # remove unnecessary .xml file

t_euend = Sys.time()
cat('clim 1981-2010: ')
print(t_euend - t_start)

## merge all traditional bioclim layers for period 2011-2040
filenames <- list.files(path = "R/data/CHELSA_bio_layers/2011-2040", full.names = TRUE)

merged <- rast() # initialise raster
for (f in filenames) {
    layer <- rast(f) # load layer in question
    name <- strsplit(f, split = "_") # extract variable name
    names(layer) <- name[[1]][4] # rename layer
    merged <- c(merged, layer)
}
## crop and reproject with cropped land cover as reference

# crop and reproject raster
clim_eu <- crop(merged, ext(lc_eu), snap = "out")
clim_eu <- resample(clim_eu, lc_eu, method = "bilinear")
clim_as <- crop(merged, ext(lc_as), snap = "out")
clim_as <- resample(clim_as, lc_as, method = "bilinear")

# save croped rasters
fname <- "R/data/cropped_rasters/CHELSA_bio_merged_2011-2040_eu.grd"
writeRaster(clim_eu, filename = fname, overwrite = TRUE)
unlink(paste(fname, ".aux.xml", sep = "")) # remove unnecessary .xml file
fname <- "R/data/cropped_rasters/CHELSA_bio_merged_2011-2040_as.grd"
writeRaster(clim_as, filename = fname, overwrite = TRUE)
unlink(paste(fname, ".aux.xml", sep = "")) # remove unnecessary .xml file

t_asend = Sys.time()
cat('clim 2011-2040: ')
print(t_asend - t_euend)