# This file subsets the cleaned occurence points per year and generates
# pseudo-absence points in a x radius around the presences.
# Afterwards, the corresponding bioclim and land cover values
# are extracted to create one merged dataframe for modelling

# used libraries
library(terra)

# load cleaned occurences
occs <- readRDS("R/data/occurence_data/axyridis_clean.rds")
# load reference lc layers
lc_eu <- rast("R/data/cropped_rasters/Cop_LC_2002_eu.grd")
lc_as <- rast("R/data/cropped_rasters/Cop_LC_2002_as.grd")

years <- c("2002") # unique(occs$Year)

for (y in years) {
    # extract points of certain year
    occs_y <- subset(occs, Year == y)
    # create SpatVector object
    occs_v <- vect(occs_y, geom = c("Lat", "Lon"), crs = "epsg:4326")
    # create 10 km circles around each point
    circs <- buffer(occs_v, 10000)

    absence_pts <- data.frame(
        Year = y,
        Presence = "absent",
        Area = occs_y$Area
    )


    for (c in circs) { # for each circle separately
        n_tot <- 100 # amount of background points per occurrence point
        pts <- spatSample(c, n_tot) # generate random points inside
        pts <- cbind(pts, extract(lc_ref$lccs_class, pts)) # extract lc values
        pts <- pts[lccs_class != 210, ] # test for lc = water
        # have to be more than 0.5km away from other presences
        pts <- pts[apply(distance(pts, occs_v), 1, function(x) any(x > 500)), ]

        # generate replacements if needed
        while (nrow(pts) < n_tot) {
            n <- n_tot - nrow(pts)
            pts_n <- spatSample(circs, n)
            pts_n <- cbind(pts_n, extract(lc_ref$lccs_class, pts_n))
            pts_n <- pts_n[lccs_class != 210]
            pts_n <- pts_n[apply(distance(pts_n, occs_v), 1, function(x) any(x > 500)), ]
            pts <- rbind(pts, pts_n)
        }
        pts$lccs_class <- NULL # remove lc column
    }
}
points()
