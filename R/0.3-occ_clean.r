# This file prepares the downloaded occurence data from gbif,
# checking the dataset for low quality points and removing them
library(CoordinateCleaner)
library(terra)

# load raw gbif data
axyridis_raw <- read.csv("R/data/occurence_data/Harmonia-axyridis_gbif_raw.csv",
    header = TRUE,
    sep = "\t"
)

# only take relevant columns
axyridis_raw <- axyridis_raw[c("decimalLatitude", "decimalLongitude", "year", "coordinateUncertaintyInMeters")]
names(axyridis_raw) <- c("Lat", "Lon", "Year", "CoordUncert")

# remove relevant NAs in data
# remove occurences after 2022 and with larger uncertainty than 1 km
axyridis_sub <- subset(axyridis_raw, complete.cases(axyridis_raw) &
    Year < 2023 & CoordUncert <= 1000)
cat(length(axyridis_raw$Year) - length(axyridis_sub$Year), "points removed with subset", "\n")

# subset data into eu and as spatial extent for further cleaning
# land cover layers vor reference (xmin, xmax, ymin, ymax)
eu <- ext(rast("R/data/cropped_rasters/Cop_LC_2016_eu.grd"))
as <- ext(rast("R/data/cropped_rasters/Cop_LC_2016_as.grd"))
# some issue with crop and SpatVector with terra, so ugly subset
axyridis_eu <- subset(axyridis_sub, Lon > eu[1] & Lon < eu[2] & Lat > eu[3] & Lat < eu[4])
axyridis_eu <- subset(axyridis_eu, Year >= 1991) #first invasion according to EASIN (contradiction in Roy et al. 2015 ?)
axyridis_as <- subset(axyridis_sub, Lon > as[1] & Lon < as[2] & Lat > as[3] & Lat < as[4])

croplength = length(axyridis_eu$Year) + length(axyridis_as$Year)
cat(length(axyridis_sub$Year) - croplength, "points removed with crop", "\n")
# clean coordinates using CoordinateCleaner
# create species columns because clean_coordinates asks for that for some reason
axyridis_eu$species <- "EU"
axyridis_as$species <- "AS"

tests_used <- c("capitals", "centroids", "duplicates", "equal", "institutions", "outliers", "seas", "zeros")

options(warn = -1) # remove legacy warnings
eu_clean <- clean_coordinates(axyridis_eu, tests = tests_used, lat = "Lat", lon = "Lon")
as_clean <- clean_coordinates(axyridis_as, tests = tests_used, lat = "Lat", lon = "Lon")
options(warn = 0) # allow warnings again

# remove flagged points
axyridis_eu <- axyridis_eu[eu_clean$.summary, ]
axyridis_as <- axyridis_as[as_clean$.summary, ]

# merge and rename species column to Area for easier subsetting later
axyridis_clean <- rbind(axyridis_eu, axyridis_as)
names(axyridis_clean)[5] <- "Area"

cat(croplength - length(axyridis_clean$Year), "points removed with clean", "\n")
cat(length(axyridis_clean$Year), 'cleaned points remaining in total', "\n")
cat("EU: ", length(axyridis_eu$Year), "| AS: ", length(axyridis_as$Year), "\n")
saveRDS(axyridis_clean, file = "R/data/occurence_data/axyridis_clean.rds")
