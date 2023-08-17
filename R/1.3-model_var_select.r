# this file selects the bioclim variables to be used in modelling,
# using variance inflation factors and a pca for land cover

# used libraries
library(dplyr)
library(corrplot)
library(car) # for vif

tot_time <- Sys.time()
# load pa data
pa_ext <- readRDS("R/data/occurence_data/axyridis_pa_vals_extracted.rds")
pa_ext <- subset(pa_ext, Year == 2020)
bio <- select(pa_ext, starts_with("bio")) # only take bio vars for scaling

# add squared values as variable columns
bio_sq <- bio^2
names(bio_sq) <- paste0(names(bio), "_2")
bio <- cbind(bio, bio_sq)
# scale variables (-mean, /stdev)
bio_sc <- data.frame(scale(bio)) # scale all variables (-mean, /stdev)
bio_scaling <- rbind("mean" = colMeans(bio), "sd" = apply(bio, 2, sd))

bio_sc$Presence <- factor(pa_ext$Presence) # add factorized presence column
bio_sc <- bio_sc[, order(colnames(bio_sc))] # sort columns for readability

# create full glm for all bio vars
var_mod <- glm(Presence ~ ., data = bio_sc, family = "binomial")
# compute vifs for all variables
vifs <- vif(var_mod)

# drop components until no vif is > 10
while (max(vifs) > 10) {
    # find variable with largest vif
    highest <- names(which(vifs == max(vifs)))

    # drop quadratic term before linear
    if (! grepl("_2", highest) & paste0(highest, "_2") %in% names(bio_sc)) {
        # drop quadratic variable from dataframe
        bio_sc <- bio_sc[, -which(names(bio_sc) %in% paste0(highest, "_2"))]
        cat("dropped", paste0(highest, "_2"), "\n")
    }else {
        # drop variable from dataframe
        bio_sc <- bio_sc[, -which(names(bio_sc) %in% highest)]
        cat("dropped", highest, "\n")
    }

    # calculate new model and vifs
    var_mod <- glm(Presence ~ ., data = bio_sc, family = "binomial")
    vifs <- vif(var_mod)
}
#final bio vifs 2002
#   bio10    bio15    bio16    bio18     bio2     bio3     bio8     bio9
#2.567260 3.667203 7.988719 4.992440 4.657536 3.273481 1.260575 1.164423

#final bio vifs 2020
#    bio1    bio14    bio15  bio15_2    bio18    bio19     bio2     bio3     bio8     bio9

#1.718275 8.240707 7.150737 6.111256 3.513818 7.519414 1.647547 2.139928 2.604866 3.227555


td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("reduced model with vifs:", td, "secs", "\n")
