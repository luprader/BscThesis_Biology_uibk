# this file creates an ensemble of sdm models for each year

# libraries uesd
library(dplyr)
library(gam) # for gam
# library(gbm) # for brt

## use 2008 and 2012 for testing

pa <- readRDS("R/data/modelling/pa_mod_vars.rds")
pa_mod <- subset(pa, Area == "as")

# select model variables from df
data <- select(pa_mod, matches("[[:digit:]]"))

# create a dummy data frame for predicting response curves
pr_data <- data.frame(matrix(0, nrow = 100, ncol = ncol(data)))
colnames(pr_data) <- colnames(data)

# plot distribution of every variable for pa
png_name <- "R/data/modelling/response_curves/hist_native.png"
png(width = 600, height = 600 * ncol(data), filename = png_name)
par(mfrow = c(ncol(data), 1), cex = 1.5)
for (v in colnames(data)) {
    b <- seq(min(pa_mod[[v]]), max(pa_mod[[v]]), length.out = 20) # bin breaks
    hist(subset(pa_mod, Presence == "absent")[[v]],
        breaks = b, main = "native value histogram", xlab = v,
        xlim = range(b), col = "grey"
    )
    hist(subset(pa_mod, Presence == "present")[[v]],
        breaks = b, main = "native value histogram", xlab = v,
        xlim = range(b), col = "grey35", add = TRUE
    )
    legend("topright", c("absent", "present"), fill = c("grey", "grey35"))
}
dev.off()

# scale data for modelling
data_sc <- data.frame(scale(data)) # - mean, / stdev
scaling <- rbind("mean" = colMeans(data), "sd" = apply(data, 2, sd))
data_sc$Presence <- as.numeric(pa_mod$Presence == "present") # p/a to 1/0

# fit glm
f <- paste0(names(data), collapse = " + ") # formula for glm
f <- as.formula(paste0("Presence ~ ", f))
mod_glm <- glm(f, data = data_sc, family = "binomial")

# plot response curves for every variable
png_name <- "R/data/modelling/response_curves/response_glm_native.png"
png(width = 600, height = 600 * ncol(data), filename = png_name)
par(mfrow = c(ncol(data), 1), cex = 1.5)
for (v in names(data)) {
    pr_data[[v]] <- seq(min(data_sc[[v]]), max(data_sc[[v]]), length.out = nrow(pr_data))
    p <- predict(mod_glm, newdata = pr_data, type = "response")
    x <- pr_data[[v]] * scaling["sd", v] + scaling["mean", v]
    plot(x, p,
        type = "l", xlab = v, ylab = "suitability",
        main = "native glm response curve",
        ylim = c(0, 1), xlim = c(min(pa_mod[[v]]), max(pa_mod[[v]]))
    )
    pr_data[[v]] <- 0
}
dev.off()

# fit gam
f <- paste0("s(", names(data), ")", collapse = " + ") # formula for gam
f <- as.formula(paste0("Presence ~ ", f))
mod_gam <- gam(f, data = data_sc, family = "binomial")

# plot response curves for every variable
png_name <- "R/data/modelling/response_curves/response_gam_native.png"
png(width = 600, height = 600 * ncol(data), filename = png_name)
par(mfrow = c(ncol(data), 1), cex = 1.5)
for (v in names(data)) {
    pr_data[[v]] <- seq(min(data_sc[[v]]), max(data_sc[[v]]), length.out = nrow(pr_data))
    p <- predict(mod_gam, newdata = pr_data, type = "response")
    x <- pr_data[[v]] * scaling["sd", v] + scaling["mean", v]
    plot(x, p,
        type = "l", xlab = v, ylab = "suitability",
        main = "native gam response curve",
        ylim = c(0, 1), xlim = c(min(pa_mod[[v]]), max(pa_mod[[v]]))
    )
    pr_data[[v]] <- 0
}
dev.off()
