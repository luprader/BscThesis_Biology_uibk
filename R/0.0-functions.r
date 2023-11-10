# this file includes all self written functions used in this project

# libraries used
library(terra)
# (lp_subdiv_pts, lp_gen_abs, lp_ext_vals, lp_clean_lc, lp_pca_proj_lc)
library(dplyr)
# (lp_gen_abs, lp_eval_mods)
library(FactoMineR)
# (lp_pca_proj, lp_pca_proj_lc)
library(rnaturalearth)
# (lp_pca_proj_lc)
library(PresenceAbsence)
# (lp_eval_mods)
library(maxnet)
# (lp_eval_mods)

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
# function evenly splitting an extent into a grid of extents.
# n -> amount of subdivisions per side (n = 1 => one split in x & y => 4 cells)
# s_ext -> the initial extent to start off of
# has to be a vector in the form of c(xmin, xmax, ymin, ymax)

# returns a vector containing the needed extents in the form of init_ext

lp_subdiv_grd <- function(n, s_ext) {
    x_step <- (s_ext[2] - s_ext[1]) / n
    y_step <- (s_ext[4] - s_ext[3]) / n
    grid_exts <- c()
    for (i in 1:n) {
        for (j in 1:n) {
            n_ext <- c(s_ext[1] + (j - 1) * x_step, s_ext[1] + j * x_step, s_ext[3] + (i - 1) * y_step, s_ext[3] + i * y_step)
            grid_exts <- rbind(grid_exts, n_ext)
        }
    }
    return(grid_exts)
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

    ## generate n_abs absence points per circle
    wc <- 0 # how often replacements had to be generated

    # for (i in seq_along(circs_rd)) {
    cat("\r", "|", year, "|") # , i, "|") # print gen progress
    c <- vect(ext(lc_ref), crs = crs(lc_ref)) # normal random extent sampling
    c$Year <- year
    c$CoordUncert <- 0
    c$Area <- pres_y$Area[1]
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
    # }

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
    # check if pts is empty
    if (nrow(pts) == 0) {
        cat("no points given")
        return(c())
    }
    # create points SpatVector for extracting
    pts_v <- vect(pts, geom = c("Lon", "Lat"), crs = "epsg:4326")

    # raster paths
    clim_p <- "R/data/cropped_rasters/CHELSA_bio_merged_"
    clim_p_ya <- paste(clim_p, y_clim, "_", area, ".tif", sep = "")
    lc_p <- "R/data/cropped_rasters/Cop_LC_"
    lc_p_ya <- paste(lc_p, y_lc, "_", area, ".tif", sep = "")

    # extract clim and lc values
    clim_l <- rast(clim_p_ya)
    pts_ext <- cbind(pts, extract(clim_l, pts_v, ID = FALSE))
    rm(clim_l)

    # extract lc as relative area covered in a buffer around each presence
    lc_l <- rast(lc_p_ya)
    lc_cs <- cellSize(lc_l) # cell size for each raster cell
    lc_b <- buffer(pts_v, 36000) # 18 km radius buffer for extraction
    # all lccs_class values used in a raster
    lc_classes <- read.csv("R/plots/Cop_LCCS_legend.csv", header = TRUE)[, 1]

    # compute relative area of each lccs_class per circle
    lc_rel <- data.frame()
    for (i in seq_along(lc_b)) {
        circ <- lc_b[i]
        # extract lccs_class and cell size for single circle
        c_vals <- cbind(extract(lc_l, circ), extract(lc_cs, circ, ID = FALSE))
        for (c in lc_classes) {
            colnm <- paste0("lc_", c) # column name for class
            c_vals_c <- subset(c_vals, lccs_class == c)
            lc_rel[i, colnm] <- sum(c_vals_c$area) / sum(c_vals$area)
        }
    }
    rm(lc_l)
    pts_ext <- cbind(pts_ext, lc_rel)
    return(pts_ext)
    gc()
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
    # project lc onto pca axes
    lc_proj <- predict.PCA(pca_res, lc)$coord
    colnames(lc_proj) <- paste0("lc", seq_len(ncol(lc_proj))) # rename

    return(lc_proj)
}
################################################################################
# function projecting Cop LC layers to pca dims, writing to file

# pca_res -> result of a PCA to pass onto lp_pca_proj
# year -> year for which the Cop LC layer should be projected
# cont_eu -> SpatVector of European country extents to remove oceans

# returns the year for which projection was computed
# writes the projected raster as a new file with filename = fn

lp_pca_proj_lc <- function(pca_res, year, cont_eu) {
    # load Cop LC layer of year
    lc_r <- rast(paste0("R/data/cropped_rasters/Cop_lc_", year, "_eu.tif"))

    # remove oceans from LC layer
    lc_c <- crop(lc_r, cont_eu, mask = TRUE)

    # create destination file name
    fn <- paste0("R/data/modelling/pca_rasters/pca_", year, "_eu_water.tif")

    # project raster onto pca axes
    app(lc_c, lp_pca_proj, pca_res = pca_res, filename = fn, overwrite = TRUE)
    cat("projected", year, "lc raster onto pca axes \n")

    rm(list = c("lc_r", "lc_c"))
    return(year)
}
################################################################################
# function generating response curves and evaluating all given models
# all models should have been trained on the same dataset

# glm -> object of class "glm"
# gam -> object of class "gam" (library "gam")
# brt -> object of class "gbm" (library "gbm")
# maxent -> object of class "maxnet" (library "maxnet")
# data -> data with which to evaluate the models
# ys -> years of which to use EU data for evaluation (each year used separately)
# must be a vector containing numeric values between 2002 and 2022
# sc -> scaling used to scale the training data prior to building models
# png_name -> filename of the response curve png

# returns a list containing the presence.absence.accuracy() results for each
# model and each value of ys as well as an ensemble prediction weighted with TSS
# generates a png with response curves for all variables in the trained range

lp_eval_mods <- function(m_glm, m_gam, m_brt, m_max, data, ys, sc, png_name) {
    ## plot data distribution compared to response curves
    m_data <- m_glm$data # data used to train the models
    # create a dummy data frame for predicting response curves
    pr_data <- data.frame(matrix(0, nrow = 100, ncol = ncol(m_data) - 1))
    colnames(pr_data) <- colnames(m_data)[colnames(m_data) != "pres"]

    png(width = 3000, height = 600 * ncol(pr_data), filename = png_name)
    par(mfrow = c(ncol(pr_data), 5), cex = 1.5)
    for (v in colnames(pr_data)) {
        # plot data distribution histograms
        x <- m_data[[v]] * sc["sd", v] + sc["mean", v] # rescaled for plot
        b <- seq(min(x), max(x), length.out = 20) # breaks
        hist(subset(x, m_data$pres == 0),
            breaks = b, main = "training data distribution", xlab = v,
            xlim = range(b), col = "grey"
        )
        hist(subset(x, m_data$pres == 1),
            breaks = b, main = "training data distribution", xlab = v,
            xlim = range(b), col = "grey35", add = TRUE
        )
        legend("topright", c("absent", "present"), fill = c("grey", "grey35"))

        # plot model response curves
        pr_data[[v]] <- seq(min(m_data[[v]]), max(m_data[[v]]), length.out = nrow(pr_data))
        x <- pr_data[[v]] * sc["sd", v] + sc["mean", v] # rescaled for plot
        # glm response
        p <- predict(m_glm, newdata = pr_data, type = "response")
        plot(x, p,
            type = "l", ylim = c(0, 1), xlim = range(x),
            xlab = v, ylab = "suitability", main = "glm response curve"
        )
        # gam response
        p <- predict(m_gam, newdata = pr_data, type = "response")
        plot(x, p,
            type = "l", ylim = c(0, 1), xlim = range(x),
            xlab = v, ylab = "suitability", main = "gam response curve"
        )
        # brt response
        p <- predict(m_brt, newdata = pr_data, type = "response")
        plot(x, p,
            type = "l", ylim = c(0, 1), xlim = range(x),
            xlab = v, ylab = "suitability", main = "brt response curve"
        )
        # maxent response
        p <- predict(m_max, newdata = pr_data, type = "logistic")
        plot(x, p,
            type = "l", ylim = c(0, 1), xlim = range(x),
            xlab = v, ylab = "suitability", main = "maxent response curve"
        )
        pr_data[[v]] <- 0 # set variable constant again
    }
    dev.off()

    ## evaluate model accuracy for used training data
    # data for evaluation
    e_data <- select(m_data, matches("[[:digit:]]"))
    for (v in colnames(e_data)) {
        e_data[[v]] <- e_data[[v]] * sc["sd", v] + sc["mean", v]
    }
    # data frame for evaluation
    th_data <- data.frame(id = seq_len(nrow(e_data)), pres = m_data$pres)
    # model predictions for e_data
    th_data$glm <- predict(m_glm, newdata = m_data, type = "response")
    th_data$gam <- predict(m_gam, newdata = m_data, type = "response")
    th_data$brt <- predict(m_brt, newdata = m_data, type = "response")
    th_data$max <- predict(m_max, newdata = m_data, type = "logistic")
    # threshhold optimising mean of sensitivity and specificity
    ths <- optimal.thresholds(th_data, opt.methods = 3)
    # compute accuracy measurements (PCC, sens, spec, Kappa) for each model
    ma <- list()
    for (m in 1:4) {
        ma[[m]] <- presence.absence.accuracy(th_data, which.model = m, ths[[m + 1]])
    }

    # create TSS weighted ensemble
    tss <- c()
    # get tss for each model
    for (m in 1:4) {
        sens <- ma[[m]]$sensitivity
        spec <- ma[[m]]$specificity
        tss <- c(tss, sens + spec - 1)
    }
    # tss[tss<0] = 0 # if tss is negative, exclude from weighting
    # get weighted average prediction with tss
    th_data$ens <- apply(th_data[, 3:6], 1, weighted.mean, w = tss)
    # compute performance of ensemble
    th <- optimal.thresholds(th_data, which.model = 5, opt.methods = 3)
    ma[[5]] <- presence.absence.accuracy(th_data, which.model = 5, th[1, 2])
    res_t <- ma
    print(res_t)

    ## evaluate model accuracies against European data of ys
    res_y <- c() # initialize results vector
    for (y in ys) {
        # data for evaluation
        data_y <- subset(data, Area == "eu" & Year == y)
        e_data <- select(data_y, matches("[[:digit:]]"))
        for (v in colnames(e_data)) { # scale with scale from m_data
            e_data[[v]] <- (e_data[[v]] - sc["mean", v]) / sc["sd", v]
        }
        # data frame for evaluation
        th_data <- data.frame(id = seq_len(nrow(e_data)), pres = data_y$Pres)
        # model predictions for e_data
        th_data$glm <- predict(m_glm, newdata = e_data, type = "response")
        th_data$gam <- predict(m_gam, newdata = e_data, type = "response")
        th_data$brt <- predict(m_brt, newdata = e_data, type = "response")
        th_data$max <- predict(m_max, newdata = e_data, type = "logistic")
        # threshhold optimising mean of sensitivity and specificity
        ths <- optimal.thresholds(th_data, opt.methods = 3)
        # compute accuracy measurements (PCC, sens, spec, Kappa) for each model
        ma <- list()
        for (m in 1:4) {
            ma[[m]] <- presence.absence.accuracy(th_data, which.model = m, ths[[m + 1]])
        }

        # create TSS weighted ensemble
        tss <- c()
        # get tss for each model
        for (m in 1:4) {
            sens <- ma[[m]]$sensitivity
            spec <- ma[[m]]$specificity
            tss <- c(tss, sens + spec - 1)
        }
        # tss[tss<0] = 0 # if tss is negative, exclude from weighting
        # get weighted average prediction with tss
        th_data$ens <- apply(th_data[, 3:6], 1, weighted.mean, w = tss)
        # compute performance of ensemble
        th <- optimal.thresholds(th_data, which.model = 5, opt.methods = 3)
        ma[[5]] <- presence.absence.accuracy(th_data, which.model = 5, th[1, 2])

        # merge to other ys
        res_y <- rbind(res_y, ma)

        # save th_data for re examination
        if (length(ys) == 1) {
            fname <- paste0("R/data/modelling/th_data_mods/th_data_y_", y, ".rds")
        } else {
            fname <- paste0("R/data/modelling/th_data_mods/th_data_nt_", y, ".rds")
        }
        saveRDS(th_data, file = fname)
    }
    res_y <- rbind(res_y, res_t) # add training tss results
    rownames(res_y) <- c(ys, "trained")
    return(res_y)
}
