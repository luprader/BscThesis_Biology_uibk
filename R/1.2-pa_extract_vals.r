# This file extracts the correct land cover and bioclim values for each year,
# using all presence and generated absence points.

# used libraries
library(terra)
library(parallel)
library(foreach)
library(doParallel)
source("R/0.0-functions.r", encoding = "UTF-8") # self written functions used

tot_time <- Sys.time()
# load generated absences
pa <- readRDS("R/data/occurrence_data/axyridis_pa.rds")

# only use half the data
pa_thin = c()
for (y in unique(pa$Year)) {
    sub = subset(pa, Year == y)
    pa_thin = rbind(pa_thin, sub[sample(nrow(sub), as.integer(nrow(sub)/2)), ])
}
pa = pa_thin

# initialize extracted pa dataframe
pa_ext <- data.frame()

# prepare for parallelization
years <- 2002:2022 # for iteration of foreach
cl <- makeCluster(detectCores() - 2)
# load libraries in cl
clusterEvalQ(cl, library(terra))
registerDoParallel(cl)

# for each geographical area
for (area in unique(pa$Area)) {
    # subset to area
    pa_a <- subset(pa, Area == area)
    print(area)

    # parallelized for loop over years
    pa_ys <- foreach(y = years, .combine = rbind, .inorder = FALSE) %dopar% {
        # which bioclim time frame and land cover year to use
        if (y <= 2010) { # 2002-2010
            y_clim <- "1981-2010"
            y_lc <- y
        } else if (y > 2010 && y <= 2020) { # 2011-2020
            y_clim <- "2011-2040"
            y_lc <- y
        } else { # 2021-2022
            y_clim <- "2011-2040"
            y_lc <- 2020
        }

        # extract values for subset
        pa_y <- lp_ext_vals(subset(pa_a, Year == y), y_clim, y_lc, area)
        return(pa_y)
    }
    # add to total dataframe
    pa_ext <- rbind(pa_ext, pa_ys)
}
stopCluster(cl)

# save extracted dataframe
saveRDS(pa_ext, file = "R/data/occurrence_data/axyridis_pa_vals_extracted.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("value extraction finished", td, "secs", "\n")