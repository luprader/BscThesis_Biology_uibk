# this file creates an ensemble of sdm models for each year

# libraries uesd
library(dplyr)
library(gam) # for gam
library(gbm) # for brt
library(maxnet) # for maxent
source("R/0.0-functions.r", encoding = "UTF-8")

set.seed(4326) # consistent randomness

# load modelling data
pa <- readRDS("R/data/modelling/pa_mod_vars.rds")
# model with native data
pa_mod <- subset(pa, Area == "as")

# select model variables from df
data <- select(pa_mod, matches("[[:digit:]]"))
# scale data for modelling
data_sc <- data.frame(scale(data)) # - mean, / stdev
sc <- rbind("mean" = colMeans(data), "sd" = apply(data, 2, sd))
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
m_brt <- gbm(f,
    data = data_sc, distribution = "bernoulli", verbose = FALSE,
    n.trees = 1500, interaction.depth = 2, shrinkage = 0.01
)

# fit maxent
f <- maxnet.formula(data_sc$pres, select(data_sc, !pres), classes = "default")
m_max <- maxnet(data_sc$pres, select(data_sc, !pres), formula = f)

# evaluate native models for all years
pname <- "R/plots/response_curves/native_mod_resp.png"
rnt <- lp_eval_mods(m_glm, m_gam, m_brt, m_max, pa, 2002:2022, sc, pname)
saveRDS(rnt, file = "R/data/modelling/eval_results/eval_mod_native.rds")
head(rnt)
# build and evaluate models built iteratively
# parallel for loop?
rys <- c() # initialize results
years <- c(2004, 2008) # for testing
for (y in years) {
    t_time <- Sys.time()
    # load modelling data
    pa <- readRDS("R/data/modelling/pa_mod_vars.rds")
    # subset to eu data up to y
    pa_mod <- subset(pa, Area == "eu" & Year <= y)

    # select model variables from df
    data <- select(pa_mod, matches("[[:digit:]]"))
    # scale data for modelling
    data_sc <- data.frame(scale(data)) # - mean, / stdev
    sc <- rbind("mean" = colMeans(data), "sd" = apply(data, 2, sd))
    data_sc$pres <- pa_mod$Pres # p/a to 1/0

    # fit glm
    a <- Sys.time()
    f <- paste0(names(data), collapse = " + ") # formula for glm
    f <- as.formula(paste0("pres ~ ", f))
    m_glm <- glm(f, data = data_sc, family = "binomial")
    cat("glm", difftime(Sys.time(), a, units = "secs")[[1]], "\n")

    # fit gam
    a <- Sys.time()
    f <- paste0("s(", names(data), ")", collapse = " + ") # formula for gam
    f <- as.formula(paste0("pres ~ ", f))
    m_gam <- gam(f, data = data_sc, family = "binomial")
    cat("gam", difftime(Sys.time(), a, units = "secs")[[1]], "\n")

    # fit brt
    a <- Sys.time()
    f <- paste0(names(data), collapse = " + ") # formula for brt
    f <- as.formula(paste0("pres ~ ", f))
    m_brt <- gbm(f,
        data = data_sc, distribution = "bernoulli", verbose = FALSE,
        n.trees = 1500, interaction.depth = 2, shrinkage = 0.01
    )
    cat("brt", difftime(Sys.time(), a, units = "secs")[[1]], "\n")

    # fit maxent
    a <- Sys.time()
    f <- maxnet.formula(data_sc$pres, select(data_sc, !pres), classes = "default")
    m_max <- maxnet(data_sc$pres, select(data_sc, !pres), formula = f)
    cat("max", difftime(Sys.time(), a, units = "secs")[[1]], "\n")

    # evaluate models for following year and 2022
    a <- Sys.time()
    pname <- paste0("R/plots/response_curves/", y, "_mod_resp.png")
    ry <- lp_eval_mods(m_glm, m_gam, m_brt, m_max, pa, c(y + 1, 2022), sc, pname)
    rys <- rbind(rys, ry)
    cat("eval", difftime(Sys.time(), a, units = "secs")[[1]], "\n")
    print(Sys.time() - t_time)
}
colnames(rys) <- c("glm_f", "gam_f", "brt_f", "max_f", "glm_22", "gam_22", "brt_22", "max_22")
