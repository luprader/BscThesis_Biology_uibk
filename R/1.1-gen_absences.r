# This file generates pseudo-absence points in a defined radius and a minimum
# distance away from the presence points.

# used libraries
library(terra)
library(dplyr)
source("R/0.0-functions.r", encoding = "UTF-8") # self written functions used

tot_time <- Sys.time()
# load cleaned occurences
occs <- readRDS("R/data/occurence_data/axyridis_clean.rds")
occs <- subset(occs, Year == 2019)
# load reference lc layers
lc_eu <- rast("R/data/cropped_rasters/Cop_LC_2002_eu.grd")
lc_as <- rast("R/data/cropped_rasters/Cop_LC_2002_as.grd")
lc_ref <- merge(lc_eu, lc_as)
# ext(lc_ref) returns:
# SpatExtent : -25, 150, 19.9916666666667, 72 (xmin, xmax, ymin, ymax)
# lc_ref extent as vector for subdiv function
ref_ext_v <- c(-25, 150, 19.9916666666667, 72)

# create SpatVector object of presence points
occs_v <- vect(occs, geom = c("Lon", "Lat"), crs = crs(lc_ref))

# generate subdividing extents to improve computation time
subexts <- lp_subdiv_pts(occs_v, end_ptcount = 2000, ref_ext_v)
# generate absences for each sub extent and merge
pa <- data.frame() # initialize pa df
count = 0
for (e in seq_len(nrow(subexts))) {
    cat("\r", "|", count, "|")
    count = count + 1
    ext_e <- vect(ext(subexts[e, ]), crs = crs(lc_ref))
    # subset occurences inside the extent in question
    occs_c <- crop(occs_v, ext(ext_e))
    # generate absences inside the extent
    pa_e <- lp_gen_abs(occs_c, n_abs = 5, min_d = 1000, max_d = 10000, lc_ref)
    # merge with already computed points
    pa <- rbind(pa, pa_e)
}

# save complete pa data separately
saveRDS(pa, file = "R/data/occurence_data/axyridis_pa.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("\n", "absence generation completed ", td, "secs", "\n")
