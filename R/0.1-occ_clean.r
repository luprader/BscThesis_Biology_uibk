# This file prepares the downloaded occurence data from gbif,
# checking the dataset for low quality points and removing them
library(CoordinateCleaner)

## clean occurence presence data
# load raw gbif data
axyridis_raw <- read.csv("R/data/Harmonia-axyridis_gbif_raw.csv",
  header = TRUE,
  sep = "\t"
)
# remove relevant NAs in data
axyridis <- subset(axyridis_raw, !is.na(decimalLongitude & decimalLatitude & year &
  coordinateUncertaintyInMeters) & year < 2023 &
  coordinateUncertaintyInMeters < 1000)
# clean coordinates using CoordinateCleaner
axyridis_clean <- clean_coordinates(axyridis,
  lat = "decimalLatitude",
  lon = "decimalLongitude",
  countries = "countryCode", capitals_rad = 5000
)
# save data
saveRDS(axyridis_clean, file = "R/data/axyridis_global.rds")
print("completed occurence cleaning")
