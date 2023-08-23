# this file selects the bioclim variables to be used in modelling,
# using variance inflation factors and a pca for land cover

# used libraries
library(dplyr)
library(car) # for vif
library(FactoMineR) # for pca
library(factoextra) # for pca visualization
library(ggpubr)

tot_time <- Sys.time()
# load pa data
pa_ext <- readRDS("R/data/occurrence_data/axyridis_pa_vals_extracted.rds")
pa_ext <- subset(pa_ext, Year == 2002) # only take 2002 as reference

## compute pca for lccs_class
# transform each present class into separate column
lc <- select(pa_ext, lccs_class)
for (v in unique(lc$lccs_class)) {
    lc <- cbind(lc, c_name = as.numeric(lc$lccs_class == v)) # make binary
    lc <- rename(lc, !!paste0("lc_", v) := c_name) # rename to variable
}
lc_bin <- lc[, -1] # separate colvals (1, 0)

# calculate pca for binary lc values
lc_pca <- PCA(lc_bin, ncp = 10, scale.unit = FALSE, graph = FALSE)

# extract pca dims with cumulative variance closest to 80%
cutoff <- which.min(abs(lc_pca$eig[, 3] - 80))
lc_pca_dims <- as.data.frame(lc_pca$ind$coord[, 1:cutoff])
colnames(lc_pca_dims) <- paste0("lc", seq_len(ncol(lc_pca_dims)))

# calculate pca with cutoff for saving
lc_pca <- PCA(lc_bin, ncp = cutoff, scale.unit = FALSE, graph = FALSE)
# save lc pca results
saveRDS(lc_pca, file = "R/data/modelling/var_select_lc_pca_res.rds")
rm(lc_pca)

# add squared bioclim values as variable columns
bio <- select(pa_ext, starts_with("bio"))
bio_sq <- bio^2
names(bio_sq) <- paste0(names(bio), "_2")

# all potential variables
vars <- cbind(lc_pca_dims, bio, bio_sq)
# scale variables (-mean, /stdev)
vars_sc <- data.frame(scale(vars)) # scale all variables (-mean, /stdev)
vars_scaling <- rbind("mean" = colMeans(vars), "sd" = apply(vars, 2, sd))

vars_sc$Presence <- factor(pa_ext$Presence) # add factorized presence column

# create full glm for all vars
var_mod <- glm(Presence ~ ., data = vars_sc, family = "binomial")
# compute vifs for all variables
vifs <- vif(var_mod)
print(vifs)

# drop components until no vif is > 10
while (max(vifs) > 10) {
    # find variable with largest vif
    highest <- names(which(vifs == max(vifs)))

    # drop quadratic term before linear
    if (!grepl("_2", highest) & paste0(highest, "_2") %in% names(vars_sc)) {
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

# save final vifs and var scaling
saveRDS(vifs, file = "R/data/modelling/var_select_vifs.rds")
saveRDS(vars_scaling, file = "R/data/modelling/var_select_vars_scaling.rds")
td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("reduced model with vifs:", td, "secs", "\n")
