# this file prepares the pa data and rasters to use them in modelling with the
# selected variables
# the correct bioclim variables are subset and quadratic terms are added
# the landcover values are projected onto the chosen pca dimensions
# new rasters are generated with the respective pca values to enable prediction

# used libraries
library(terra)
library(dplyr)
library(FactoMineR)
library(parallel)
library(foreach)
library(doParallel)
source("R/0.0-functions.r", encoding = "UTF-8")

tot_time = Sys.time()
# load var_select results
lc_pca <- readRDS("R/data/modelling/var_select_lc_pca_res.rds")
vs_vifs <- as.data.frame(readRDS("R/data/modelling/var_select_vifs.rds"))

# prepare pa with selected vars
pa_ext <- readRDS("R/data/occurrence_data/axyridis_pa_vals_extracted.rds")

lc <- data.matrix(select(pa_ext, lccs_class))
lc_proj <- as.data.frame(lp_pca_proj(lc, lc_pca))
# subset selected lc vars
lc_vars <- select(lc_proj, any_of(rownames(vs_vifs)))

# subset selected bioclim variables
bio_vars <- select(pa_ext, any_of(rownames(vs_vifs)))
# add any present squared variables
sq_names <- grep("_2", rownames(vs_vifs), value = TRUE)
bio_vars_sq <- select(bio_vars, sub("_2", "", sq_names))^2
colnames(bio_vars_sq) <- sq_names

# merge all variables to pa and save
pa_mod_vars <- cbind(pa_ext[, 1:6], bio_vars, bio_vars_sq, lc_vars)
saveRDS(pa_mod_vars, file = "R/data/modelling/pa_mod_vars.rds")

# generate rasters with lccs_class projected to pca axes
# prepare for parallelization
years <- c(2002:2020) # for iteration of foreach
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, lapply(c("terra", "FactoMineR"), library, character.only = TRUE))
registerDoParallel(cl)
# parallelized for loop
ys <- foreach(y = years, .inorder = FALSE) %dopar% {
    yp <- lp_pca_proj_lc(lc_pca, y)
}
stopCluster(cl)

# remove unnecessary .xml files
unlink(paste0("R/data/modelling/pca_rasters/", "*.aux.xml"))

td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("data prep for model prediction finished:", td, "secs", "\n")
