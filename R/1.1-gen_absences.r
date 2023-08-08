# This file generates pseudo-absence points in a defined radius and a minimum
# distance away from the presence points.

# used libraries
library(terra)
library(dplyr)

n_tot <- 5 # amount of background points per occurrence point
ao_range <- 10000 # range for absence generation in m
pdist <- 500 # minimum distance of absences to a presence point in m

# load cleaned occurences
pres <- readRDS("R/data/occurence_data/axyridis_clean.rds")
# load reference lc layers
lc_eu <- rast("R/data/cropped_rasters/Cop_LC_2002_eu.grd")
lc_as <- rast("R/data/cropped_rasters/Cop_LC_2002_as.grd")
lc_ref <- merge(lc_eu, lc_as)

# initialize df for absences
names <- c("Lat", "Lon", "Year", "CoordUncert", "Area", "Presence")
abs <- data.frame(matrix(nrow = 0, ncol = length(names)))
colnames(abs) <- names

# create SpatVector object of presence points
pres_v <- vect(pres, geom = c("Lon", "Lat"), crs = crs(lc_ref))
# create circles around each point
circs <- buffer(pres_v, ao_range)

stime <- Sys.time()
pb <- txtProgressBar(min = 0, max = length(circs), initial = 0, style = 3, width = 50)
prog <- 0

set.seed(4326) # have consistent randomness
for (i in seq_along(circs)) {
    c <- circs[i] # select one single circle each time
    pts <- spatSample(c, n_tot) # generate random points inside
    # extract lc values
    pts <- cbind(pts, extract(lc_ref, pts, method = "simple", ID = FALSE))
    # test for lc = water or NA (out of cropped area)
    pts <- pts["lccs_class" != 210 & !is.na("lccs_class"), ]
    # have to be more than 0.5km away from other presences
    pts <- pts[apply(distance(pts, pres_v), 1, function(x) any(x > pdist)), ]
    pts$lccs_class <- NULL # remove lc column

    # generate replacements if needed
    whilecount <- 0
    while (nrow(pts) < n_tot) {
        whilecount <- whilecount + 1
        n <- n_tot - nrow(pts)
        pts_n <- spatSample(circs, n)
        pts_n <- cbind(pts_n, extract(lc_ref$lccs_class, pts_n))
        pts_n <- pts_n["lccs_class" != 210 & !is.na("lccs_class"), ]
        pts_n <- pts_n[apply(distance(pts_n, pres_v), 1, function(x) any(x > pdist)), ]
        pts$lccs_class <- NULL # remove lc column
        pts <- rbind(pts, pts_n)
    }
    pts_df <- as.data.frame(pts, geom = "XY") # turn SpatVector back to df
    pts_df <- rename(pts_df, c("Lon" = "x", "Lat" = "y"))
    pts_df$Presence <- "absent" # add presence column
    # add generated points to total dataframe
    abs <- rbind(abs, pts_df)
    # map progress
    prog <- prog + 1
    setTxtProgressBar(pb, prog)
}
close(pb)

cat("presences:", nrow(pres), "absences:", nrow(abs), "\n")
cat("tot whilecount:", whilecount, " ")
print(Sys.time() - stime) 

# save po, ao and pa data separately
saveRDS(pres, file = "R/data/occurence_data/axyridis_po.rds")
saveRDS(abs, file = "R/data/occurence_data/axyridis_ao.rds")
pa = rbind(pres, abs)
saveRDS(pa, file = "R/data/occurence_data/axyridis_pa.rds")