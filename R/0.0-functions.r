# this file includes all self written functions used in this project

# libraries used
library(terra)
# (lp_subdiv_pts, lp_gen_abs, lp_ext_vals, lp_clean_lc, lp_pca_proj_lc)
library(dplyr)
# (lp_gen_abs)
library(FactoMineR)
# (lp_pca_proj, lp_pca_proj_lc)

################################################################################
# function subdividing point extents until all extents have no more than the
# desired pointcount

# end_ptcount -> the desired maximal point count per extent
# has to be a number
# points -> points to subset to the desired point count
# has to be a SpatVector object (terra) with extent <= init_ext
# init_ext -> the initial extent to start off of
# has to be a vector in the form of c(xmin, xmax, ymin, ymax)

# returns a vector containing the needed extents in the form of init_ext

lp_subdiv_pts <- function(points, end_ptcount, init_ext) {
    # check if points can be subdivided (goal smaller than given)
    ext_v_init <- vect(ext(init_ext), crs = "epsg:4326")
    ptcount_init <- nrow(crop(points, ext(ext_v_init)))
    if (ptcount_init <= end_ptcount) {
        cat("Error(lp_subdiv_pts): points \u2265 end_ptcount", "\n")
        return(init_ext) # subsetting will not help, init_ext is "the best"
    }

    print("starting subdivision")
    starting_time <- Sys.time()
    # initialize result lists
    good_exts <- c()
    bad_exts <- c()

    # set init_ext as current extent to investigate
    n_ext <- init_ext
    # get x and y span of extent
    xdiff <- n_ext[2] - n_ext[1]
    ydiff <- n_ext[4] - n_ext[3]

    # decide wether to split x or y in half
    if (xdiff > ydiff) {
        xstep <- xdiff / 2
        new_exts <- rbind(
            c(n_ext[1], n_ext[1] + xstep, n_ext[3], n_ext[4]), # left of split
            c(n_ext[1] + xstep, n_ext[2], n_ext[3], n_ext[4]) # right of split
        )
    } else {
        ystep <- ydiff / 2
        new_exts <- rbind(
            c(n_ext[1], n_ext[2], n_ext[3], n_ext[3] + ystep), # left of split
            c(n_ext[1], n_ext[2], n_ext[3] + ystep, n_ext[4]) # right of split
        )
    }
    # add new_exts to all extents still needing to be split
    bad_exts <- rbind(bad_exts, new_exts)

    repeat {
        # start testing until all bad exts are gone
        # select new ext and check for last element
        if (is.null(nrow(bad_exts))) {
            n_ext <- bad_exts
        } else {
            n_ext <- bad_exts[1, ]
        }

        # compute pointcount in this extent
        n_ext_v <- vect(ext(n_ext), crs = "epsg:4326")
        ptcount_n_ext <- nrow(crop(points, ext(n_ext_v)))
        # check if desired pointcount is reached
        if (ptcount_n_ext <= end_ptcount) {
            if (ptcount_n_ext != 0) {
                # if not empty, move extent from bad to good
                good_exts <- rbind(good_exts, n_ext)
            }
        } else {
            # too many points, split extent again
            # get x and y span of extent
            xdiff <- n_ext[2] - n_ext[1]
            ydiff <- n_ext[4] - n_ext[3]
            # decide wether to split x or y in half
            if (xdiff > ydiff) {
                xstep <- xdiff / 2
                new_exts <- rbind(
                    c(n_ext[1], n_ext[1] + xstep, n_ext[3], n_ext[4]), # left
                    c(n_ext[1] + xstep, n_ext[2], n_ext[3], n_ext[4]) # right
                )
            } else {
                ystep <- ydiff / 2
                new_exts <- rbind(
                    c(n_ext[1], n_ext[2], n_ext[3], n_ext[3] + ystep), # left
                    c(n_ext[1], n_ext[2], n_ext[3] + ystep, n_ext[4]) # right
                )
            }
            # add new_exts to all extents still needing to be split
            bad_exts <- rbind(bad_exts, new_exts)
        }
        # remove old extent and check if last element is good
        if (is.null(nrow(bad_exts))) {
            bad_exts <- 0
            break # all extents are good, function is finished
        } else {
            bad_exts <- bad_exts[-1, ]
        }
        # current amount of good extents
        cat("\r", "Good extents: ", nrow(good_exts), " ")
    }
    # function end statements
    td <- difftime(Sys.time(), starting_time, units = "secs")[[1]]
    cat("\r", "Final extents: ", nrow(good_exts), "|", td, "secs", "\n")

    return(good_exts)
}

################################################################################
# function generating pseudoabsence/background points with a minimum and maximum
# distance away from other presences.
# The absences are generated for a specific year but distance to all years is
# taken into account.

# pres -> presences for which to compute absence points
# has to be a SpatVector object (terra)
# year -> year for which to generate absences
# has to be a value present in pres$Year
# n_abs -> number of absences to compute per presence in m
# has to be an integer
# min_d -> minimum distance of absences to presences in m
# has to be a number
# max_d -> maximum distance of absences to presences
# has to be a number
# lc_ref -> land cover raster to test for water/NA
# has to be a SpatRaster object (terra) with numerical values, water = 210

# returns a dataframe with the following columns:
# (Lon, Lat, Year, CoordUncert, Area, Presence)
# contains the absences generated for year

lp_gen_abs <- function(pres, year, n_abs, min_d, max_d, lc_ref) {
    starting_time <- Sys.time()

    # check if pres is empty
    if (length(pres) == 0) {
        stop(paste("no presences given \n"))
    }

    # initialize dataframe for output
    names <- c("Lat", "Lon", "Year", "CoordUncert", "Area", "Presence")
    ao <- data.frame(matrix(nrow = 0, ncol = length(names)))
    colnames(ao) <- names

    # subset presences to year in question
    pres_y <- subset(pres, pres$Year == year)
    # check if pres_y is empty
    if (length(pres_y) == 0) {
        cat("no presences for", year, "\n")
        return(ao)
    }

    # create circles around each point
    # maximum distance to presences of the year to generate for
    circs_r <- buffer(pres_y, max_d)
    # minimum distance to presences of all years
    circs_d <- buffer(pres, min_d)
    # remove min_d circles from all max_d circles
    circs_rd <- erase(circs_r, circs_d)

    ## generate n_abs absence points per circle
    wc <- 0 # how often replacements had to be generated

    #for (i in seq_along(circs_rd)) {
        cat("\r", "|", year, "|")#, i, "|") # print gen progress
        c <- vect(ext(pres_y), crs = crs(pres_y)) # normal random extent sampling
        c$Year = year
        c$CoordUncert = 0
        c$Area = pres_y$Area[1]
        pts <- spatSample(c, n_abs * nrow(pres_y)) # generate random points inside
        # extract lc values
        pts <- cbind(pts, extract(lc_ref, pts, ID = FALSE))
        # test for lc = water or NA (out of cropped area)
        pts <- pts[pts$lccs_class != 210 & !is.na(pts$lccs_class), ]
        pts$lccs_class <- NULL # remove lc column

        # generate replacements if needed
        while (nrow(pts) < (n_abs * nrow(pres_y))) {
            wc <- wc + 1
            n <- n_abs * nrow(pres_y) - nrow(pts)
            pts_n <- spatSample(c, n)
            pts_n <- cbind(pts_n, extract(lc_ref, pts_n, ID = FALSE))
            pts_n <- pts_n[pts_n$lccs_class != 210 & !is.na(pts_n$lccs_class), ]
            pts_n$lccs_class <- NULL # remove lc column
            pts <- rbind(pts, pts_n)
        }
        pts_df <- as.data.frame(pts, geom = "XY") # turn SpatVector to df
        pts_df <- rename(pts_df, c("Lon" = "x", "Lat" = "y"))
        # add generated points to total dataframe
        ao <- rbind(ao, pts_df)
    #}

    ao$Presence <- "absent"
    n_pres <- length(pres_y)
    n_abs <- nrow(ao)
    # function end statements
    td <- difftime(Sys.time(), starting_time, units = "secs")[[1]]
    cat(
        "presences:", n_pres, "absences:", n_abs,
        "wc:", wc, "|", td, "secs", "\n"
    )

    # return generated presence-absence dataframe
    return(ao)
    # free some memory
    rm(list = c("circs_r", "circs_d", "circs_rd", "c", "pts", "pts_n"))
}

################################################################################
# function extracting all point values of the desired CHELSA bioclim or
# Copernicus LCCS rasters.

# pts -> points for which to extract the values
# has to be a dataframe of the following form:
# (Lon, Lat, Year, CoordUncert, Area, Presence)
# y_clim -> time frame for which CHELSA data should be extracted
# has to be one of the following strings: "1981-2010" or "2011-2040"
# y_lc -> year for which Copernicus LCCS data should be extracted
# has to be a year between 2002 and 2020
# area -> cropped area for which to extract the data
# has to be one of the following strings: "eu" or "as"

# returns a dataframe with the extracted values as 20 new columns

lp_ext_vals <- function(pts, y_clim, y_lc, area) {
    # create points SpatVector for extracting
    pts_v <- vect(pts, geom = c("Lon", "Lat"), crs = "epsg:4326")

    # raster paths
    clim_p <- "R/data/cropped_rasters/CHELSA_bio_merged_"
    clim_p_ya <- paste(clim_p, y_clim, "_", area, ".tif", sep = "")
    lc_p <- "R/data/cropped_rasters/Cop_LC_"
    lc_p_ya <- paste(lc_p, y_lc, "_", area, ".tif", sep = "")

    # extract clim and lc values
    clim_l <- rast(clim_p_ya)
    lc_l <- rast(lc_p_ya)
    pts_ext <- cbind(pts, extract(clim_l, pts_v, ID = FALSE))
    pts_ext <- cbind(pts_ext, extract(lc_l, pts_v, ID = FALSE))

    return(pts_ext)
}
################################################################################
# function to remove water and NA when cleaning occurrences

# points -> points for which to extract the values
# has to be a dataframe of the following form:
# (Lon, Lat, Year, CoordUncert, Area)
# y_lc -> year for which Copernicus LCCS data should be extracted
# has to be a year between 2002 and 2020
# area -> cropped area for which to extract the data
# has to be one of the following strings: "eu" or "as"

# returns a dataframe with all points in water (210 lccs_class) or NA removed

lp_clean_lc <- function(points, y_lc, area) {
    # load landcover for year and area
    lc_p <- "R/data/cropped_rasters/Cop_LC_"
    lc_p_ya <- paste(lc_p, y_lc, "_", area, ".tif", sep = "")
    lc_l <- rast(lc_p_ya)

    # extract lc values and remove water or NA
    points_v <- vect(points, geom = c("Lon", "Lat"), crs = crs(lc_l))
    points <- cbind(points, extract(lc_l, points_v, ID = FALSE))
    points_r <- subset(points, lccs_class != 210 & !is.na(points$lccs_class))
    points_r$lccs_class <- NULL

    return(points_r)
    # free some memory
    rm(c("lc_l", "points_v"))
}
################################################################################
# function computing pca projected coordinates for a lccs_class vector

# lc -> vector/matrix of lccs_class values to project
# pca_res -> result from a PCA

# returns a matrix with the separate values of lc on all dimensions of pca_res

lp_pca_proj <- function(lc, pca_res) {
    # extract the lccs_classes used by the pca
    lccs_factors <- as.factor(sub("lc_", "", rownames(pca_res$var$contrib)))

    for (i in seq_along(lccs_factors)) {
        v <- lccs_factors[i]
        lc <- cbind(lc, c_name = as.numeric(lc == v)) # make binary column
        colnames(lc)[ncol(lc)] <- paste0("lc_", v) # rename column to variable
    }
    lc_bin <- lc[, -1] # remove original lccs_class column

    # project lc onto pca axes
    lc_proj <- predict.PCA(pca_res, lc_bin)$coord
    colnames(lc_proj) <- paste0("lc", seq_len(ncol(lc_proj))) # rename

    return(lc_proj)
}
################################################################################
# function projecting Cop LC layers to pca dims, writing to file

# pca_res -> result of a PCA to pass onto lp_pca_proj
# year -> year for which the Cop LC layer should be projected

# returns the year for which projection was computed
# writes the projected raster as a new file with filename = fn

lp_pca_proj_lc <- function(pca_res, year) {
    # load Cop LC layer of year
    lc_r <- rast(paste0("R/data/cropped_rasters/Cop_lc_", year, "_eu.tif"))

    # create destination file name
    fn <- paste0("R/data/modelling/pca_rasters/pca_", year, "_eu.tif")

    # project raster onto pca axes
    app(lc_r, lp_pca_proj, pca_res = pca_res, filename = fn, overwrite = TRUE)
    cat("projected", year, "lc raster onto pca axes \n")

    return(year)
}
################################################################################
