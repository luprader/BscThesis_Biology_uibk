# this file creates an ensemble of sdm models for each year

# libraries uesd
library(terra)
library(dplyr)

pa <- readRDS("R/data/modelling/pa_mod_vars.rds")
pa_mod = subset(pa, Area == "as")

data = select(pa_mod, starts_with("bio"))
data_sc = data.frame(scale(data)) # scale all variables (-mean, /stdev)
scaling = rbind("mean" = colMeans(data), "sd" = apply(data, 2, sd))
data_sc$Presence = as.numeric(pa_mod$Presence == 'present')

mod_glm <- glm(Presence ~ ., data = data_sc, family = "binomial")

## create a dummy data frame for predicting response curves
pr_data <- data.frame(matrix(0, nrow = 100, ncol = ncol(data_sc)))
colnames(pr_data) = colnames(data_sc)
pr_data$Presence = NULL

par(mfrow = c(3, 2))
for (v in names(pr_data)) {
    pr_data[[v]] <- seq(-1, 1, length.out = nrow(pr_data))
    theta1 <- predict(mod_glm, newdata = pr_data, type = "response")
    x <- pr_data[[v]] * scaling["sd", v] + scaling["mean", v]
    y <- data_sc[[v]] * scaling["sd", v] + scaling["mean", v]
    plot(x, theta1, type = "l", xlab = v, ylab = "suitability", 
    ylim = c(0,1), xlim = range(y))
    points(y, data_sc$Presence)
    hist(y[which(data_sc$Presence == 0)], breaks = 20, main = NULL, xlab = v,
    xlim = range(y))
    hist(y[which(data_sc$Presence == 1)], breaks = 20, main = NULL, xlab = v,
    xlim = range(y), col = "green", add = TRUE)
    pr_data[[v]] <- 0
}
dev.off()

# create prediction raster from pca and clim
 data_r <- rast("R/data/modelling/pca_rasters/pca_2002_eu.grd")
 clim_r <- rast("R/data/cropped_rasters/CHELSA_bio_merged_1981-2010_eu.grd")
 data_r <- c(subset(clim_r, grep("bio", colnames(data), val = TRUE)), data_r)
 data_r_sc <- scale(data_r, center = scaling["mean", ], scale = scaling["sd", ])

 pred_rast <- predict(data_r_sc, mod_glm, filename = "test_pred.tif", overwrite = TRUE)
 pred_rast = rast("test_pred.tif")
 plot(pred_rast)
 pred_class = classify(pred_rast, c(0, 1))
 plot(pred_class, col = "green")