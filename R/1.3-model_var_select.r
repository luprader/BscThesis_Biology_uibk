# this file selects the bioclim variables to be used in modelling,
# using variance inflation factors and a pca for land cover

# used libraries
library(dplyr)
library(car) # for vif
library(FactoMineR) # for pca
library(gam) # for gam
source("R/0.0-functions.r", encoding = "UTF-8")

tot_time <- Sys.time()
# load pa data
pa_ext <- readRDS("R/data/occurrence_data/axyridis_pa_vals_extracted.rds")
pa_mod <- subset(pa_ext, Year == 2002) # only take 2002 as reference

## compute pca for lccs_class

# generate independent points for pca, using 2002 as reference year
lc_ref <- rast("R/data/cropped_rasters/Cop_LC_2002_eu.tif")

## generate random points in eu extent
ptsmax <- 5000 # desired pointcount
pts <- spatSample(lc_ref, ptsmax, xy = TRUE) # generate random points inside
pts <- subset(pts, lccs_class != 210 & !is.na(pts$lccs_class))

while (nrow(pts) < ptsmax) {
    pts_n <- spatSample(lc_ref, ptsmax - nrow(pts), xy = TRUE)
    pts_n <- subset(pts_n, lccs_class != 210 & !is.na(pts_n$lccs_class))
    pts <- rbind(pts, pts_n)
}
# (Lon, Lat, Year, CoordUncert, Area, Presence)
pts$lccs_class <- 2002
pts$Coord <- 0
pts$Area <- "eu"
pts$Presence <- "present"
colnames(pts) <- c("Lon", "Lat", "Year", "CoordUncert", "Area", "Presence")

# get relative lc class abundance for generated points
lc <- select(lp_ext_vals(pts, "1981-2010", 2002, "eu"), starts_with("lc"))

# calculate pca for relative class abundance
lc_pca <- PCA(lc, ncp = 10, scale.unit = FALSE, graph = FALSE)

# extract pca dims with cumulative variance closest to 80%
cutoff <- which.min(abs(lc_pca$eig[, 3] - 80))
lc_pca_dims <- as.data.frame(lc_pca$ind$coord[, 1:cutoff])
colnames(lc_pca_dims) <- paste0("lc", seq_len(ncol(lc_pca_dims)))

# calculate pca with cutoff for saving
lc_pca <- PCA(lc, ncp = cutoff, scale.unit = FALSE, graph = FALSE)
# save lc pca results
saveRDS(lc_pca, file = "R/data/modelling/var_select_lc_pca_res.rds")

# project lc for pa from 2002
pa_mod_lc <- select(pa_mod, starts_with("lc"))
pa_lc_proj <- as.data.frame(lp_pca_proj(pa_mod_lc, lc_pca))

# all potential variables
bio <- select(pa_mod, starts_with("bio"))
vars <- cbind(pa_lc_proj, bio)
# scale variables (-mean, /stdev)
vars_sc <- data.frame(scale(vars)) # scale all variables (-mean, /stdev)
vars_sc$Pres <- factor(pa_mod$Presence) # add factorized presence column

# create full gam for all vars
t <- Sys.time()
f <- paste0("s(", names(select(vars_sc, !"Pres")), ")", collapse = " + ") # formula for gam
f <- as.formula(paste0("Pres ~ ", f))
var_mod <- gam(f, data = vars_sc, family = "binomial")
# compute vifs for all variables
vifs <- vif(var_mod)

# drop components until no vif is > 10
while (max(vifs) > 10) {
    # find variable with largest vif
    highest <- gsub(".*s\\((.+)\\)*.", "\\1", names(which(vifs == max(vifs))))

    # drop variable from dataframe
    vars_sc <- vars_sc[, -which(names(vars_sc) %in% highest)]
    cat("dropped", highest, "\n")

    # calculate new model and vifs
    f <- paste0("s(", names(select(vars_sc, !"Pres")), ")", collapse = " + ") # formula for gam
    f <- as.formula(paste0("Pres ~ ", f))
    var_mod <- gam(f, data = vars_sc, family = "binomial")
    vifs <- vif(var_mod)
}
vifs <- as.data.frame(vifs)
# save final vifs
saveRDS(vifs, file = "R/data/modelling/var_select_vifs.rds")

# create dataframe with final variables for modelling
lc <- data.matrix(select(pa_ext, starts_with("lc")))
lc_proj <- as.data.frame(lp_pca_proj(lc, lc_pca))
# get names of final variables used
fin_vars <- gsub(".*s\\((.+)\\)*.", "\\1", rownames(vifs))
# subset selected lc vars
lc_vars <- select(lc_proj, any_of(fin_vars))
# subset selected bioclim variables
bio_vars <- select(pa_ext, any_of(fin_vars))

# merge all variables to pa and save
pa_mod_vars <- cbind(
    Area = pa_ext$Area, Year = pa_ext$Year,
    Pres = as.numeric(pa_ext$Presence == "present"),
    bio_vars, lc_vars
)
saveRDS(pa_mod_vars, file = "R/data/modelling/pa_mod_vars.rds")

td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("reduced model with vifs, selected modelling variables:", td, "secs", "\n")
