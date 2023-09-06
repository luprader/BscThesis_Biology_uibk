# this file selects the bioclim variables to be used in modelling,
# using variance inflation factors and a pca for land cover

# used libraries
library(dplyr)
library(car) # for vif
library(FactoMineR) # for pca
source("R/0.0-functions.r", encoding = "UTF-8")

tot_time <- Sys.time()
# load pa data
pa_ext <- readRDS("R/data/occurrence_data/axyridis_pa_vals_extracted.rds")
pa_mod <- subset(pa_ext, Year == 2002) # only take 2002 as reference

## compute pca for lccs_class
# transform each present class into separate column
lc <- select(pa_mod, starts_with("lc"))

# calculate pca for binary lc values
lc_pca <- PCA(lc, ncp = 10, scale.unit = FALSE, graph = FALSE)

# extract pca dims with cumulative variance closest to 90%
cutoff <- which.min(abs(lc_pca$eig[, 3] - 90))
lc_pca_dims <- as.data.frame(lc_pca$ind$coord[, 1:cutoff])
colnames(lc_pca_dims) <- paste0("lc", seq_len(ncol(lc_pca_dims)))

# calculate pca with cutoff for saving
lc_pca <- PCA(lc, ncp = cutoff, scale.unit = FALSE, graph = FALSE)
# save lc pca results
saveRDS(lc_pca, file = "R/data/modelling/var_select_lc_pca_res.rds")

# add squared bioclim values as variable columns
bio <- select(pa_mod, starts_with("bio"))
bio_sq <- bio^2
names(bio_sq) <- paste0(names(bio), "_2")

# all potential variables
vars <- cbind(lc_pca_dims, bio, bio_sq)
# scale variables (-mean, /stdev)
vars_sc <- data.frame(scale(vars)) # scale all variables (-mean, /stdev)
vars_sc$Presence <- factor(pa_mod$Presence) # add factorized presence column

# create full glm for all vars
var_mod <- glm(Presence ~ ., data = vars_sc, family = "binomial")
# compute vifs for all variables
vifs <- vif(var_mod)

# drop components until no vif is > 10
while (max(vifs) > 10) {
    # find variable with largest vif
    highest <- names(which(vifs == max(vifs)))

    # drop quadratic term before linear
    if (!grepl("_2", highest) && paste0(highest, "_2") %in% names(vars_sc)) {
        # drop quadratic variable from dataframe
        vars_sc <- vars_sc[, -which(names(vars_sc) %in% paste0(highest, "_2"))]
        cat("dropped", paste0(highest, "_2"), "\n")
    } else {
        # drop variable from dataframe
        vars_sc <- vars_sc[, -which(names(vars_sc) %in% highest)]
        cat("dropped", highest, "\n")
    }
    # calculate new model and vifs
    var_mod <- glm(Presence ~ ., data = vars_sc, family = "binomial")
    vifs <- vif(var_mod)
}
vifs <- as.data.frame(vifs)
# save final vifs
saveRDS(vifs, file = "R/data/modelling/var_select_vifs.rds")

# create dataframe with final variables for modelling
lc <- data.matrix(select(pa_ext, starts_with("lc")))
lc_proj <- as.data.frame(lp_pca_proj(lc, lc_pca))
# subset selected lc vars
lc_vars <- select(lc_proj, any_of(rownames(vifs)))

# subset selected bioclim variables
bio_vars <- select(pa_ext, any_of(rownames(vifs)))
# add any present squared variables
sq_names <- grep("_2", rownames(vifs), value = TRUE)
bio_vars_sq <- select(bio_vars, sub("_2", "", sq_names))^2
colnames(bio_vars_sq) <- sq_names

# merge all variables to pa and save
pa_mod_vars <- cbind(
    Area = pa_ext$Area, Year = pa_ext$Year,
    Pres = as.numeric(pa_ext$Presence == "present"),
    bio_vars, bio_vars_sq, lc_vars
)
saveRDS(pa_mod_vars, file = "R/data/modelling/pa_mod_vars.rds")

td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("reduced model with vifs, selected modelling variables:", td, "secs", "\n")
