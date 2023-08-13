# this file includes all self written functions used in this project

# function_name <- function(parameters){
#  function body
# }

# libraries used
library(terra)
library(dplyr)

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
    if(ptcount_init <= end_ptcount) {
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
    cat("\r","Final extents: ",nrow(good_exts),"|",td,"secs","\n")

    return(good_exts)
}

################################################################################
# function generating pseudoabsence/background points with a minimum and maximum
# distance away from other presences.

# pres -> presences for which to compute absence points
# has to be a SpatVector object (terra)
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

lp_gen_abs <- function(pres, n_abs, min_d, max_d, lc_ref) {
    set.seed(4326) # have consistent randomness
    # check if pres is empty
    if (length(pres) == 0) {
        print('no presences given')
        break
    }
    Sys.sleep(1)
    starting_time <- Sys.time()

    # initialize dataframe for output
    names <- c("Lat", "Lon", "Year", "CoordUncert", "Area", "Presence")
    pa <- data.frame(matrix(nrow = 0, ncol = length(names)))
    colnames(pa) <- names

    # create circles around each point
    circs_r <- buffer(pres, max_d) # maximum distance to presences
    circs_d <- buffer(pres, min_d) # minimum distance to presences
    # remove min_d circles from all max_d circles
    circs_rd <- erase(circs_r, circs_d)

    ## generate n_abs absence points per circle
    wc <- 0 # how often replacements had to be generated
    cat("|gen|") # print gen progress
    for (i in seq_along(circs_rd)) {
        cat("\r", "|", i, "|") # print gen progress
        c <- circs_rd[i]
        pts <- spatSample(c, n_abs) # generate random points inside
        # extract lc values
        pts <- cbind(pts, extract(lc_ref, pts, ID = FALSE))
        # test for lc = water or NA (out of cropped area)
        pts <- pts["lccs_class" != 210 & !is.na("lccs_class"), ]
        pts$lccs_class <- NULL # remove lc column

        # generate replacements if needed
        while (nrow(pts) < n_abs) {
            wc <- wc  + 1
            n <- n_abs - nrow(pts)
            pts_n <- spatSample(c, n)
            pts_n <- cbind(pts_n, extract(lc_ref, pts_n, ID = FALSE))
            pts_n <- pts_n["lccs_class" != 210 & !is.na("lccs_class"), ]
            pts_n$lccs_class <- NULL # remove lc column
            pts <- rbind(pts, pts_n)
        }
        pts_df <- as.data.frame(pts, geom = "XY") # turn SpatVector to df
        pts_df <- rename(pts_df, c("Lon" = "x", "Lat" = "y"))
        # add generated points to total dataframe
        pa <- rbind(pa, pts_df)
    }

    # turn presence SpatVector into dataframe
    pres_df <- as.data.frame(pres, geom = "XY")
    pres_df <- rename(pres_df, c("Lon" = "x", "Lat" = "y"))
    # add presence labels and add presences to total dataframe
    pres_df$Presence <- "present"
    pa$Presence <- "absent"
    n_pres <- nrow(pres_df)
    n_abs <- nrow(pa)
    pa <- rbind(pres_df, pa)
    # function end statements
    td <- difftime(Sys.time(), starting_time, units = "secs")[[1]]
    cat("presences:", n_pres, "absences:", n_abs,
        "wc:", wc, "|", td, "secs")

    # return generated presence-absence dataframe
    return(pa)
}

################################################################################
# function computing the occupied niche of given points with environmental data
# using (Broennimann et al. 2011)
# use ecospat?
################################################################################
