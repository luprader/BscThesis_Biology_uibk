# this file creates an ensemble of sdm models for each year

# libraries uesd
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
    n.trees = 1500, interaction.depth = 1, shrinkage = 0.01
)

# fit maxent
f <- maxnet.formula(data_sc$pres, select(data_sc, !pres), classes = "lht")
m_max <- maxnet(data_sc$pres, select(data_sc, !pres), formula = f)

# evaluate native models for all years
pname <- "R/plots/response_curves/native_mod_resp.png"
rnt <- lp_eval_mods(m_glm, m_gam, m_brt, m_max, pa, 2002:2022, sc, pname)
rm(list = c("m_glm", "m_gam", "m_brt", "m_max")) # clear up memory
saveRDS(rnt, file = "R/data/modelling/eval_mod_native.rds")

# build and evaluate models built iteratively
# prepare for parallelization
years <- 2002:2020 # for iteration of foreach
cl <- makeCluster(detectCores() - 1)
# load libraries in cl
clusterEvalQ(cl, lapply(c("gam", "gbm", "maxnet"), library, character.only = TRUE))
registerDoParallel(cl)
# parallelized for loop
rys <- foreach(y = years, .combine = rbind) %dopar% {
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
        n.trees = 1500, interaction.depth = 1, shrinkage = 0.01
    )

    # fit maxent
    f <- maxnet.formula(data_sc$pres, select(data_sc, !pres), classes = "lht")
    m_max <- maxnet(data_sc$pres, select(data_sc, !pres), formula = f)

    # evaluate models for following year and 2022
    pnm <- paste0("R/plots/response_curves/", y, "_mod_resp.png")
    ry <- lp_eval_mods(m_glm, m_gam, m_brt, m_max, pa, c(y + 1, 2022), sc, pnm)

    rm(list = c("m_glm", "m_gam", "m_brt", "m_max")) # clear up memory
    return(ry)
}
stopCluster(cl)

# save evaluation results
colnames(rys) <- c("glm", "gam", "brt", "max", "ens")
saveRDS(rys, file = "R/data/modelling/eval_mod_years.rds")