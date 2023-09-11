# This file generates pseudo-absence points in a defined radius and a minimum
# distance away from the presence points.

# used libraries
library(terra)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
source("R/0.0-functions.r", encoding = "UTF-8") # self written functions used

tot_time <- Sys.time()
set.seed(4326) # consistent randomness
# load cleaned occurrences
occs <- readRDS("R/data/occurrence_data/axyridis_clean.rds")
occs <- subset(occs, Year >= 2002) # remove insignificant historic presences

ao <- data.frame() # initialize ao df
# values for absence generation
years <- c(2002:2022)
n_abs <- 5
min_d <- 1000
max_d <- 18000
# path for land cover rasters
lc_p <- "R/data/cropped_rasters/Cop_LC_"

# create SpatVector object of presence points
ref <- rast("R/data/cropped_rasters/Cop_LC_2002_as.tif")
occs_v <- vect(occs, geom = c("Lon", "Lat"), crs = crs(ref))
rm(ref)

# generate for asia
cat("as: \n")
pres_v <- subset(occs_v, occs_v$Area == "as")
for (y in years) {
    # choose correct lc reference
    if (y > 2020) {
        lc_ref <- rast(paste0(lc_p, 2020, "_as.tif"))
    } else {
        lc_ref <- rast(paste0(lc_p, y, "_as.tif"))
    }
    # generate absences
    ao_y <- lp_gen_abs(pres_v, y, n_abs, min_d, max_d, lc_ref)
    ao <- rbind(ao, ao_y)
    rm(lc_ref)
}

# generate for europe
cat("eu, sub extents in parallel: \n")
pres_v <- subset(occs_v, occs_v$Area == "eu")
# save to access in foreach
saveRDS(pres_v, file = "R/data/occurrence_data/pres_v.rds")
rm(occs_v)
# subset europe
# ext() for eu returns:
# SpatExtent : -25, 65, 34.9916666666667, 72 (xmin, xmax, ymin, ymax)
# t_ref extent as vector for subdiv function
t_ext <- c(-25, 65, 34.9916666666667, 72)
subexts <- lp_subdiv_pts(pres_v, 20000, t_ext)

# prepare for parallelization
e_s <- seq_len(nrow(subexts)) # for iteration of foreach
cl <- makeCluster(detectCores())
# load libraries in cl
clusterEvalQ(cl, lapply(c("terra", "dplyr"), library, character.only = TRUE))
registerDoParallel(cl)
# parallelized for loop
ao_eu <- foreach(e = e_s, .combine = rbind, .inorder = FALSE) %dopar% {
    ext_e <- ext(subexts[e, ]) # get extent
    pres_v_c <- crop(readRDS("R/data/occurrence_data/pres_v.rds"), ext_e)
    ao_e <- data.frame() # initialize ao df
    for (y in years) {
        # choose correct lc reference
        if (y > 2020) {
            lc_ref_c <- crop(rast(paste0(lc_p, 2020, "_eu.tif")), ext_e)
        } else {
            lc_ref_c <- crop(rast(paste0(lc_p, y, "_eu.tif")), ext_e)
        }
        # generate absences
        ao_y <- lp_gen_abs(pres_v_c, y, n_abs, min_d, max_d, lc_ref_c)
        ao_e <- rbind(ao_e, ao_y)
        rm(lc_ref_c)
    }
    return(ao_e)
}
stopCluster(cl)
unlink("R/data/occurrence_data/pres_v.rds") # remove saved pres_v again

ao <- rbind(ao, ao_eu) # merge all ao

# create pa dataframe
po <- occs
po$Presence <- "present"
pa <- rbind(po, ao)

# save complete pa data
saveRDS(pa, file = "R/data/occurrence_data/axyridis_pa.rds")
saveRDS(subexts, file = "R/data/plotting/axyridis_abs_gen_subexts.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("\n", "absence generation completed", td, "secs", "\n")
