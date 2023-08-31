#!!attention!!
# this file cannot be run from the initial data files, one would have to supply
# the merged chelsa bioclim layers and the land cover .nc files

# this file resamples all land cover layers to a coarser resolution
# and crops it to the geographical extents of europe and the native range
library(terra)
library(sf)
library(ncdf4)

## CHELSA has resolution of ~ 1km, Copernicus LCC 300m
## => aggregate 3x3 raster cells of LCC

# create boundary boxes for europe and asia (native range)
sf_use_s2(FALSE) # to get rectangular intersect boxes instead of a projection
# extent of europe (chosen)
eu_box <- c(xmin = -25, ymin = 38, xmax = 65, ymax = 72)
class(eu_box) <- "bbox"
eu_box <- st_as_sfc(eu_box)
eu_box <- st_as_sf(eu_box, crs = 4326)

# extent of native range (Orlova-Bienkowskaja, Ukrainsky & Brown, 2015)
as_box <- c(xmin = 70, ymin = 20, xmax = 150, ymax = 65)
class(as_box) <- "bbox"
as_box <- st_as_sfc(as_box)
as_box <- st_as_sf(as_box, crs = 4326)


# load gobal LCC .tif files and crop them to native and europe extent
filenames <- list.files(path = "R/data/Cop_LC_raw", full.names = TRUE)
stime_tot <- Sys.time()
for (f in filenames) {
    stime_f <- Sys.time()
    layer <- rast(f, "lccs_class") # load file

    # create new filenames
    fn <- strsplit(f, "[-]")
    y <- fn[[1]][8] # extract year of land cover layer
    dest <- "R/data/cropped_rasters" # destination directory
    n_eu <- paste("Cop_LC_", y, "_eu.tif", sep = "") # eu filename
    n_as <- paste("Cop_LC_", y, "_as.tif", sep = "") # as filename

    # crop to eu and as extent and resample LC layer to coarser resolution
    layer_cr <- crop(layer, ext(eu_box), snap = "out")
    layer_agg <- aggregate(layer_cr, fact = 3, fun = "modal", na.rm = TRUE)
    fn <- file.path(dest, n_eu)
    writeRaster(layer_agg, filename = fn, overwrite = TRUE)

    layer_cr <- crop(layer, ext(as_box), snap = "out")
    layer_agg <- aggregate(layer_cr, fact = 3, fun = "modal", na.rm = TRUE)
    fn <- file.path(dest, n_as)
    writeRaster(layer_agg, filename = fn, overwrite = TRUE)

    etime_f <- Sys.time()
    cat(y, ": ")
    print(etime_f - stime_f)
}
etime_tot <- Sys.time()
cat("Total LC: ")
print(etime_tot - stime_tot)
