# this file creates an ensemble of sdm models for each year

# libraries uesd
library(dplyr)
library(gam) # for gam
library(gbm) # for brt
library(maxnet) # for maxent
source("R/0.0-functions.r", encoding = "UTF-8")

## use 2008 and 2012 for testing

# load modelling data
pa <- readRDS("R/data/modelling/pa_mod_vars.rds")
# model with native data
pa_mod <- subset(pa, Area == "as")

# select model variables from df
data <- select(pa_mod, matches("[[:digit:]]"))
# scale data for modelling
data_sc <- data.frame(scale(data)) # - mean, / stdev
scaling <- rbind("mean" = colMeans(data), "sd" = apply(data, 2, sd))
data_sc$pres <- pa_mod$Pres # p/a to 1/0

# fit glm
f <- paste0(names(data), collapse = " + ") # formula for glm
f <- as.formula(paste0("pres ~ ", f))
m_glm <- glm(f, data = data_sc, family = "binomial")

# fit gam
f <- paste0("s(", names(data), ")", collapse = " + ") # formula for gam
f <- as.formula(paste0("pres ~ ", f))
m_gam <- gam(f, data = data_sc, family = "binomial")

# fit brt
f <- paste0(names(data), collapse = " + ") # formula for brt
f <- as.formula(paste0("pres ~ ", f))
m_brt <- gbm(f, data = data_sc, distribution = "bernoulli", verbose = FALSE,
             n.trees = 2000, interaction.depth = 2, shrinkage = 0.01)

# fit maxent
f <- maxnet.formula(data_sc$pres, select(data_sc, !pres), classes = "default")
m_max <- maxnet(data_sc$pres, select(data_sc, !pres), formula = f)

# evaluate models
png_name = "R/plots/response_curves/native_mod_resp.png"
cmxs = lp_eval_mods(m_glm, m_gam, m_brt, m_max, pa, 2022, scaling, png_name)
print(cmxs)