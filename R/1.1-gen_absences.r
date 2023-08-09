# This file generates pseudo-absence points in a defined radius and a minimum
# distance away from the presence points.

# used libraries
library(terra)
library(dplyr)

n_tot <- 5 # amount of background points per occurrence point
ao_range <- 10000 # range for absence generation in m
pdist <- 1000 # minimum distance of absences to a presence point in m

# load cleaned occurences
pres <- readRDS("R/data/occurence_data/axyridis_clean.rds")
# pres = subset(pres, Year == 2002)
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
circs_r <- buffer(pres_v, ao_range)
circs_d <- buffer(pres_v, pdist)
circs_d <- circs_d[, -(1:ncol(circs_d))] # remove value cols for merging
ctime <- Sys.time()
# merge distance requirements to not have unnecessary iterations
circs_d_c <- combineGeoms(circs_d[1:2], circs_d[-2]) 
circs_d_c <- combineGeoms(circs_d_c[1], circs_d_c[2])
circs_rd <- erase(circs_r, circs_d_c) # actually allowed space for absences
cat("buffer gen:")
print(Sys.time() - ctime)

set.seed(4326) # have consistent randomness
# progress bar in terminal
pb <- txtProgressBar(
    min = 0,
    max = nrow(pres),
    initial = 0,
    style = 3,
    width = 50
)
prog <- 0
gtime <- Sys.time()
## generate n_tot absence points per circle
whilecount <- 0 # how often replacements had to be generated
for (i in seq_along(circs_rd)) {
    c <- circs_rd[i]
    pts <- spatSample(c, n_tot) # generate random points inside
    # extract lc values
    pts <- cbind(pts, extract(lc_ref, pts, method = "simple", ID = FALSE))
    # test for lc = water or NA (out of cropped area)
    pts <- pts["lccs_class" != 210 & !is.na("lccs_class"), ]
    pts$lccs_class <- NULL # remove lc column

    # generate replacements if needed
    while (nrow(pts) < n_tot) {
        whilecount <- whilecount + 1
        n <- n_tot - nrow(pts)
        pts_n <- spatSample(c, n)
        pts_n <- cbind(pts_n, extract(lc_ref$lccs_class, pts_n))
        pts_n <- pts_n["lccs_class" != 210 & !is.na("lccs_class"), ]
        pts$lccs_class <- NULL # remove lc column
        pts <- rbind(pts, pts_n)
    }
    pts_df <- as.data.frame(pts, geom = "XY") # turn SpatVector back to df
    pts_df <- rename(pts_df, c("Lon" = "x", "Lat" = "y"))
    # add generated points to total dataframe
    abs <- rbind(abs, pts_df)
    # map progress
    prog <- prog + 1
    setTxtProgressBar(pb, prog)
}
close(pb)

cat("presences:", nrow(pres), "absences:", nrow(abs), "\n")
cat("tot whilecount:", whilecount, "|absence gen:")
print(Sys.time() - gtime)

pres$Presence <- "present"
abs$Presence <- "absent"
pa <- rbind(pres, abs)
# save po, ao and pa data separately
saveRDS(pres, file = "R/data/occurence_data/axyridis_po.rds")
saveRDS(abs, file = "R/data/occurence_data/axyridis_ao.rds")
saveRDS(pa, file = "R/data/occurence_data/axyridis_pa.rds")
