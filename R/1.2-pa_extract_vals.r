# This file extracts the correct land cover and bioclim values for each year,
# using all presence and generated absence points.

# used libraries
library(terra)
source("R/0.0-functions.r", encoding = "UTF-8") # self written functions used

tot_time <- Sys.time()
# load generated absences
pa <- readRDS("R/data/occurrence_data/axyridis_pa.rds")

# initialize extracted pa dataframe
pa_ext <- data.frame()

# for each geographical area
for (area in unique(pa$Area)) {
    # subset to area
    pa_a <- subset(pa, Area == area)
    print(area)

    ## climate and lc for 2002 - 2010
    for (y in 2002:2010) {
        # extract values for subset
        pa_y <- lp_ext_vals(subset(pa_a, Year == y), "1981-2010", y, area)
        # add to total dataframe
        pa_ext <- rbind(pa_ext, pa_y)

        cat("\r", y) # progress
    }

    ## climate and lc for 2011 - 2020
    for (y in 2011:2020) {
        # extract values for subset
        pa_y <- lp_ext_vals(subset(pa_a, Year == y), "2011-2040", y, area)

        # add to total dataframe
        pa_ext <- rbind(pa_ext, pa_y)

        cat("\r", y) # progress
    }

    ## climate and lc for > 2020, use 2020
    # extract values for subset
    pa_y <- lp_ext_vals(subset(pa_a, Year > 2020), "2011-2040", 2020, area)

    # add to total dataframe
    pa_ext <- rbind(pa_ext, pa_y)

    cat("\n", ">2020", "\n") # progress
}
# save extracted dataframe
saveRDS(pa_ext, file = "R/data/occurrence_data/axyridis_pa_vals_extracted.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("value extraction finished", td, "secs", "\n")
