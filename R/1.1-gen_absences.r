# This file generates pseudo-absence points in a defined radius and a minimum
# distance away from the presence points.

# used libraries
library(terra)
library(dplyr)
source("R/0.0-functions.r", encoding = "UTF-8") # self written functions used

tot_time <- Sys.time()
set.seed(4326) # consistent randomness
# load cleaned occurrences
occs <- readRDS("R/data/occurrence_data/axyridis_clean.rds")
occs <- subset(occs, Year >= 2002) # remove insignificant historic presences

ao <- data.frame() # initialize ao df
# values for absence generation
n_abs <- 5
min_d <- 1000
max_d <- 18000
# path for land cover rasters
lc_p <- "R/data/cropped_rasters/Cop_LC_"

# create SpatVector object of presence points
ref <- rast("R/data/cropped_rasters/Cop_LC_2002_as.grd")
occs_v <- vect(occs, geom = c("Lon", "Lat"), crs = crs(ref))
rm(ref)

# generate for asia
cat("as: \n")
pres_v <- subset(occs_v, occs_v$Area == "as")

for (i in 2002:2020) {
    lc_ref <- rast(paste(lc_p, i, "_as.grd", sep = ""))
    ao_y <- lp_gen_abs(pres_v, i, n_abs, min_d, max_d, lc_ref)
    ao <- rbind(ao, ao_y)
    rm(lc_ref)
}
# 2021 and 2022 using 2020 lc
lc_ref <- rast(paste(lc_p, 2020, "_as.grd", sep = ""))
ao_y <- lp_gen_abs(pres_v, 2021, n_abs, min_d, max_d, lc_ref)
ao <- rbind(ao, ao_y)
ao_y <- lp_gen_abs(pres_v, 2022, n_abs, min_d, max_d, lc_ref)
ao <- rbind(ao, ao_y)
rm(lc_ref)

# generate for europe
cat("eu: \n")
pres_v <- subset(occs_v, occs_v$Area == "eu")
rm(occs_v)
# subset europe
# ext() for eu returns:
# SpatExtent : -25, 65, 34.9916666666667, 72 (xmin, xmax, ymin, ymax)
# t_ref extent as vector for subdiv function
t_ext <- c(-25, 65, 34.9916666666667, 72)
subexts <- lp_subdiv_pts(pres_v, 20000, t_ext)


prog <- 1
for (e in seq_len(nrow(subexts))) {
    cat("||", prog, "||", "\n")
    prog <- prog + 1
    ext_e <- ext(subexts[e, ])

    for (i in 2002:2020) {
        lc_ref_c <- crop(rast(paste(lc_p, i, "_eu.grd", sep = "")), ext_e)
        pres_v_c <- crop(pres_v, ext(subexts[e, ]))
        ao_y <- lp_gen_abs(pres_v_c, i, n_abs, min_d, max_d, lc_ref_c)
        ao <- rbind(ao, ao_y)
        rm(list = c("lc_ref_c", "pres_v_c"))
    }
    # 2021 and 2022 using 2020 lc
    lc_ref <- rast(paste(lc_p, 2020, "_eu.grd", sep = ""))
    ao_y <- lp_gen_abs(pres_v, 2021, n_abs, min_d, max_d, lc_ref)
    ao <- rbind(ao, ao_y)
    ao_y <- lp_gen_abs(pres_v, 2022, n_abs, min_d, max_d, lc_ref)
    ao <- rbind(ao, ao_y)
    rm(lc_ref)
    gc()
}

# create pa dataframe
po <- occs
po$Presence <- "present"
pa <- rbind(po, ao)

# save complete pa data
saveRDS(pa, file = "R/data/occurrence_data/axyridis_pa.rds")
saveRDS(subexts, file = "R/data/plotting/axyridis_abs_gen_subexts.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("\n", "absence generation completed", td, "secs", "\n")
