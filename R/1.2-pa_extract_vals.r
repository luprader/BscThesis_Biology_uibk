# This file extracts the correct land cover and bioclim values for each year,
# using all presence and generated absence points.

# used libraries
library(terra)

# load generated absences
pa <- readRDS("R/data/occurence_data/axyridis_pa.rds")

# initialize extracted pa dataframe
pa_ext <- data.frame()

# for each geographical area
for (area in unique(pa$Area)) {
    # subset to area
    pa_a <- subset(pa, Area == area)

    ## climate and lc for < 2002
    # subset to year
    pa_y <- subset(pa_a, Year < 2002)

    # raster paths
    dir <- "R/data/cropped_rasters/"
    clim_p <- paste(dir, "CHELSA_bio_merged_1981-2010_", area, ".grd", sep = "")
    lc_p <- paste(dir, "Cop_LC_2002_", area, ".grd", sep = "")

    # create pa SpatVector for extracting
    pa_y_v <- vect(pa_y, geom = c("Lon", "Lat"), crs = "epsg:4326")

    # extract clim and lc values
    clim_l <- rast(clim_p)
    lc_l <- rast(lc_p)
    pa_y <- cbind(pa_y, extract(clim_l, pa_y_v, ID = FALSE))
    pa_y <- cbind(pa_y, extract(lc_l, pa_y_v, ID = FALSE))

    # add to total dataframe
    pa_ext <- rbind(pa_ext, pa_y)

    ## climate and lc for 2002 - 2010
    for (y in 2002:2010) {
        # subset to year
        pa_y <- subset(pa_a, Year == y)

        # raster paths
        dir <- "R/data/cropped_rasters/"
        clim_p <- paste(dir, "CHELSA_bio_merged_1981-2010_", area, ".grd", sep = "")
        lc_p <- paste(dir, "Cop_LC_", y, "_", area, ".grd", sep = "")

        # create pa SpatVector for extracting
        pa_y_v <- vect(pa_y, geom = c("Lon", "Lat"), crs = "epsg:4326")

        # extract clim and lc values
        clim_l <- rast(clim_p)
        lc_l <- rast(lc_p)
        pa_y <- cbind(pa_y, extract(clim_l, pa_y_v, ID = FALSE))
        pa_y <- cbind(pa_y, extract(lc_l, pa_y_v, ID = FALSE))

        # add to total dataframe
        pa_ext <- rbind(pa_ext, pa_y)
    }

    ## climate and lc for 2011 - 2020
    for (y in 2002:2010) {
        # subset to year
        pa_y <- subset(pa_a, Year == y)

        # raster paths
        dir <- "R/data/cropped_rasters/"
        clim_p <- paste(dir, "CHELSA_bio_merged_2011-2040_", area, ".grd", sep = "")
        lc_p <- paste(dir, "Cop_LC_", y, "_", area, ".grd", sep = "")

        # create pa SpatVector for extracting
        pa_y_v <- vect(pa_y, geom = c("Lon", "Lat"), crs = "epsg:4326")

        # extract clim and lc values
        clim_l <- rast(clim_p)
        lc_l <- rast(lc_p)
        pa_y <- cbind(pa_y, extract(clim_l, pa_y_v, ID = FALSE))
        pa_y <- cbind(pa_y, extract(lc_l, pa_y_v, ID = FALSE))

        # add to total dataframe
        pa_ext <- rbind(pa_ext, pa_y)
    }

    ## climate and lc for > 2020
    # subset to year
    pa_y <- subset(pa_a, Year > 2020)

    # raster paths
    dir <- "R/data/cropped_rasters/"
    clim_p <- paste(dir, "CHELSA_bio_merged_2011-2040_", area, ".grd", sep = "")
    lc_p <- paste(dir, "Cop_LC_2020_", area, ".grd", sep = "")

    # create pa SpatVector for extracting
    pa_y_v <- vect(pa_y, geom = c("Lon", "Lat"), crs = "epsg:4326")

    # extract clim and lc values
    clim_l <- rast(clim_p)
    lc_l <- rast(lc_p)
    pa_y <- cbind(pa_y, extract(clim_l, pa_y_v, ID = FALSE))
    pa_y <- cbind(pa_y, extract(lc_l, pa_y_v, ID = FALSE))

    # add to total dataframe
    pa_ext <- rbind(pa_ext, pa_y)
}

# save extracted dataframe
saveRDS(pa_ext, file = "R/data/occurence_data/axyridis_pa_vals_extracted.rds")
