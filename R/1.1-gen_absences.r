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


# use half the eu data
eu_sub <- subset(occs, Area == "eu")
eu_sub <- eu_sub[sample(nrow(eu_sub), as.integer(nrow(eu_sub) / 2)), ]
occs <- rbind(subset(occs, Area == "as"), eu_sub)

ao <- data.frame() # initialize ao df
# values for absence generation
years <- c(2002:2022)
n_abs <- 3
# path for land cover rasters
lc_p <- "R/data/cropped_rasters/Cop_LC_"

# create SpatVector object of presence points
ref <- rast("R/data/cropped_rasters/Cop_LC_2002_as.tif")
occs_v <- vect(occs, geom = c("Lon", "Lat"), crs = crs(ref))
rm(ref)

# generate for native range
cat("as: \n")
pres_v <- subset(occs_v, occs_v$Area == "as")

# generate basic absence part with no subdiv
n_bg <- n_abs * 1 / 3 # relative amount of random backgrount wanted
cat("base absences with no subdiv: \n")
for (y in years) {
    # choose correct lc reference
    if (y > 2020) {
        lc_ref <- rast(paste0(lc_p, 2020, "_as.tif"))
    } else {
        lc_ref <- rast(paste0(lc_p, y, "_as.tif"))
    }

    # generate absences
    ao_y <- lp_gen_abs(pres_v, y, n_bg, lc_ref)
    ao <- rbind(ao, ao_y)
    rm(lc_ref)
    gc()
}

# subset native range
cat("subdiv for bias correction: \n")
# t_ref extent for subdiv function
t_ref <- ext(rast("R/data/cropped_rasters/Cop_LC_2002_as.tif"))
t_ext <- as.vector(c(t_ref$xmin, t_ref$xmax, t_ref$ymin, t_ref$ymax))
subexts <- lp_subdiv_pts(pres_v, round(0.3 * nrow(pres_v)), t_ext)

# generate additional absences in subextents
# prepare for parallelization
# save to access in foreach
saveRDS(pres_v, file = "R/data/occurrence_data/pres_v.rds")
e_s <- seq_len(nrow(subexts)) # for iteration of foreach
cl <- makeCluster(detectCores() - 2)
# load libraries in cl
clusterEvalQ(cl, lapply(c("terra", "dplyr"), library, character.only = TRUE))
registerDoParallel(cl)
# parallelized for loop
ao_as_sd <- foreach(e = e_s, .combine = rbind, .inorder = FALSE) %dopar% {
    ext_e <- ext(subexts[e, ]) # get extent
    pres_v_c <- crop(readRDS("R/data/occurrence_data/pres_v.rds"), ext_e)
    ao_e <- data.frame() # initialize ao df
    for (y in years) {
        # choose correct lc reference
        if (y > 2020) {
            lc_ref_c <- crop(rast(paste0(lc_p, 2020, "_as.tif")), ext_e, snap = "near")
        } else {
            lc_ref_c <- crop(rast(paste0(lc_p, y, "_as.tif")), ext_e, snap = "near")
        }

        # generate absences
        ao_y <- lp_gen_abs(pres_v_c, y, n_abs - n_bg, lc_ref_c)
        ao_e <- rbind(ao_e, ao_y)
        rm(lc_ref_c)
        gc()
    }
    return(ao_e)
}

# generate for europe
cat("eu, sub extents in parallel: \n")
pres_v <- subset(occs_v, occs_v$Area == "eu")
# save to access in foreach
saveRDS(pres_v, file = "R/data/occurrence_data/pres_v.rds")
rm(occs_v)

# generate basic absence part with no subdiv
n_bg <- n_abs * 1 / 3 # relative amount of random backgrount wanted
cat("base absences with no subdiv: \n")
ao_eu <- data.frame()
for (y in years) {
    # choose correct lc reference
    if (y > 2020) {
        lc_ref <- rast(paste0(lc_p, 2020, "_eu.tif"))
    } else {
        lc_ref <- rast(paste0(lc_p, y, "_eu.tif"))
    }

    # generate absences
    ao_y <- lp_gen_abs(pres_v, y, n_bg, lc_ref)
    ao_eu <- rbind(ao_eu, ao_y)
    rm(lc_ref)
    gc()
}

# subset europe
cat("subdiv for bias correction: \n")
# t_ref extent for subdiv function
t_ref <- ext(rast("R/data/cropped_rasters/Cop_LC_2002_eu.tif"))
t_ext <- as.vector(c(t_ref$xmin, t_ref$xmax, t_ref$ymin, t_ref$ymax))
subexts <- lp_subdiv_pts(pres_v, round(0.3 * nrow(pres_v)), t_ext)

# prepare for parallelization
e_s <- seq_len(nrow(subexts)) # for iteration of foreach
cl <- makeCluster(detectCores() - 2)
# load libraries in cl
clusterEvalQ(cl, lapply(c("terra", "dplyr"), library, character.only = TRUE))
registerDoParallel(cl)
# parallelized for loop
ao_eu_sd <- foreach(e = e_s, .combine = rbind, .inorder = FALSE) %dopar% {
    ext_e <- ext(subexts[e, ]) # get extent
    pres_v_c <- crop(readRDS("R/data/occurrence_data/pres_v.rds"), ext_e)
    ao_e <- data.frame() # initialize ao df
    for (y in years) {
        # choose correct lc reference
        if (y > 2020) {
            lc_ref_c <- crop(rast(paste0(lc_p, 2020, "_eu.tif")), ext_e, snap = "near")
        } else {
            lc_ref_c <- crop(rast(paste0(lc_p, y, "_eu.tif")), ext_e, snap = "near")
        }

        # generate absences
        ao_y <- lp_gen_abs(pres_v_c, y, n_abs - n_bg, lc_ref_c)
        ao_e <- rbind(ao_e, ao_y)
        rm(lc_ref_c)
        gc()
    }
    return(ao_e)
}
stopCluster(cl)
unlink("R/data/occurrence_data/pres_v.rds") # remove saved pres_v again

ao <- rbind(ao, ao_as_sd, ao_eu, ao_eu_sd) # merge all ao
# create pa dataframe
po <- occs
po$Presence <- "present"
pa <- rbind(po, ao)

# save complete pa data
saveRDS(pa, file = "R/data/occurrence_data/axyridis_pa.rds")
saveRDS(subexts, file = "R/data/plotting/axyridis_abs_gen_subexts.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("\n", "absence generation completed", td, "secs", "\n")
