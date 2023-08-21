# this file computes and visualizes niche comparisons according to 
# (Broennimann et al. 2011), using the `ecospat` library

# used libraries
library(ecospat)
library(ade4)
library(dplyr)

tot_time <- Sys.time()
# load extracted pa data
pa <- readRDS("R/data/occurence_data/axyridis_pa_vals_extracted.rds")
pa_as <- subset(pa, Area == "as")
po_as <- subset(pa_as, Presence == "present")
pa_eu <- subset(pa, Area == "eu")
po_eu <- subset(pa_eu, Presence == "present")
bio_tot <- select(pa, starts_with("bio")) # only use bioclim vars

## compare total native to each invaded year
# pca calibrated on all datapoints
pca_tot <- dudi.pca(bio_tot, scannf = FALSE, nf = 2)
# plot pca axis contribution
fname = "R/plots/niche_comp/as_eu_pca.png"
png(width = 600, height = 600, filename = fname)
ecospat.plot.contrib(contrib = pca_tot$co, eigen = pca_tot$eig)
dev.off()
# pca scores for the whole area
scores_tot <- pca_tot$li
# pca scores for native occs
scores_po_as <- suprow(pca_tot, select(po_as, starts_with("bio")))$li
# pca scores for native occs with background
scores_pa_as <- suprow(pca_tot, select(pa_as, starts_with("bio")))$li

# calculate native dynamic occurrence densities grid
grid_as <- ecospat.grid.clim.dyn(
    glob = scores_tot,
    glob1 = scores_pa_as, sp = scores_po_as, R = 100, th.sp = 0
)

# calculate overlap between total native and each EU year
for (y in 2002:2022) {
    # take new subsets
    pa_eu_y <- subset(pa_eu, Year == y)
    po_eu_y <- subset(po_eu, Year == y)
    # pca scores for EU
    scores_po_eu_y <- suprow(pca_tot, select(po_eu_y, starts_with("bio")))$li
    scores_pa_eu_y <- suprow(pca_tot, select(pa_eu_y, starts_with("bio")))$li
    # EU grid
    grid_eu_y <- ecospat.grid.clim.dyn(
        glob = scores_tot,
        glob1 = scores_pa_eu_y, sp = scores_po_eu_y, R = 100, th.sp = 0
    )
    # overlap
    ol <- ecospat.niche.overlap(grid_as, grid_eu_y, cor = TRUE)$D
    # choose color for overlap
    if (ol == 0) {
        col_stab <- "green"
    } else {
        col_stab <- "blue"
    }
    # plot niche overlap
    fname = paste0("R/plots/niche_comp/as_eu_", y, "_niche.png")
    png(width = 600, height = 600, filename = fname)
    text <- paste("Niche Overlap native / EU", y, "| D:", ol)
    ecospat.plot.niche.dyn(grid_as, grid_eu_y,
        quant = 0.1,
        col.unf = "green", col.exp = "red", col.stab = col_stab,
        interest = 2, title = text, name.axis1 = "PC1",
        name.axis2 = "PC2"
    )
    dev.off()

    fname = paste0("R/plots/niche_comp/as_eu_", y, "_eq-sim.png")
    png(width = 1200, height = 600, filename = fname)
    par(mfrow = c(1,2))
    # niche equivalency test
    eq_test <- ecospat.niche.equivalency.test(grid_as, grid_eu_y,
        rep = 100, ncores = 2 # , alternative = "greater"
    )
    # plot eq_test
    ecospat.plot.overlap.test(eq_test, "D", paste("Equivalency Native / EU", y))
    rm(eq_test)
    # niche similarity test
    sim_test <- ecospat.niche.similarity.test(grid_as, grid_eu_y,
        rep = 100, rand.type = 2, ncores = 2 # , alternative
    )
    # plot sim_test
    ecospat.plot.overlap.test(sim_test, "D", paste("Similarity Native / EU", y))
    rm(sim_test)
    dev.off()
    rm(list = c("grid_eu_y1", "grid_eu_y2"))
}
cat("compared EU years to native range \n")
rm(list = c('pca_tot', 'scores_tot', 'scores_pa_as', 'scores_po_as', 'grid_as'))
gc()


## compare EU years to their following year
bio_eu = select(pa_eu, starts_with("bio"))
# pca calibrated on eu datapoints
pca_eu <- dudi.pca(bio_eu, scannf = FALSE, nf = 2)
# plot pca axis contribution
fname = "R/plots/niche_comp/eu_years_pca.png"
png(width = 600, height = 600, filename = fname)
ecospat.plot.contrib(contrib = pca_eu$co, eigen = pca_eu$eig)
dev.off()
# pca scores for the whole area
scores_eu_tot <- pca_eu$li

#values for starting year
pa_eu_y2 <- subset(pa_eu, Year == 2002)
po_eu_y2 <- subset(po_eu, Year == 2002)
scores_po_eu_y2 <- suprow(pca_eu, select(po_eu_y2, starts_with("bio")))$li
scores_pa_eu_y2 <- suprow(pca_eu, select(pa_eu_y2, starts_with("bio")))$li
grid_eu_y2 <- ecospat.grid.clim.dyn(
    glob = scores_eu_tot,
    glob1 = scores_pa_eu_y2, sp = scores_po_eu_y2, R = 100, th.sp = 0
)

# calculate overlap between total native and each EU year
for (y in 2003:2022) {
    # take subsets for both years
    pa_eu_y1 <- pa_eu_y2
    po_eu_y1 <- po_eu_y2
    pa_eu_y2 <- subset(pa_eu, Year == y)
    po_eu_y2 <- subset(po_eu, Year == y)
    # pca scores for each year
    scores_po_eu_y1 <- scores_po_eu_y2
    scores_pa_eu_y1 <- scores_pa_eu_y2
    scores_po_eu_y2 <- suprow(pca_eu, select(po_eu_y2, starts_with("bio")))$li
    scores_pa_eu_y2 <- suprow(pca_eu, select(pa_eu_y2, starts_with("bio")))$li
    # grids for both years
    grid_eu_y1 <- grid_eu_y2
    grid_eu_y2 <- ecospat.grid.clim.dyn(
        glob = scores_eu_tot,
        glob1 = scores_pa_eu_y2, sp = scores_po_eu_y2, R = 100, th.sp = 0
    )
    # overlap
    ol <- ecospat.niche.overlap(grid_eu_y1, grid_eu_y2, cor = TRUE)$D
    # choose color for overlap
    if (ol == 0) {
        col_stab <- "green"
    } else {
        col_stab <- "blue"
    }
    # plot niche overlap
    fname = paste0("R/plots/niche_comp/eu_", y - 1, y, "_niche.png")
    png(width = 600, height = 600, filename = fname)
    text <- paste("Niche Overlap EU", y - 1, "/", y, "| D:", ol)
    ecospat.plot.niche.dyn(grid_eu_y1, grid_eu_y2,
        quant = 0.1,
        col.unf = "green", col.exp = "red", col.stab = col_stab,
        interest = 2, title = text, name.axis1 = "PC1",
        name.axis2 = "PC2"
    )
    dev.off()

    fname = paste0("R/plots/niche_comp/eu_", y - 1, y, "_eq-sim.png")
    png(width = 1200, height = 600, filename = fname)
    par(mfrow = c(1,2))
    # niche equivalency test
    eq_test <- ecospat.niche.equivalency.test(grid_eu_y1, grid_eu_y2,
        rep = 100, ncores = 2 # , alternative = "greater"
    )
    # plot eq_test
    ecospat.plot.overlap.test(eq_test, "D", paste("Equivalency", y - 1, "/", y))
    rm(eq_test)
    # niche similarity test
    sim_test <- ecospat.niche.similarity.test(grid_eu_y1, grid_eu_y2,
        rep = 100, rand.type = 2, ncores = 2 # , alternative
    )
    # plot sim_test
    ecospat.plot.overlap.test(sim_test, "D", paste("Similarity", y - 1, "/", y))
    rm(sim_test)
    dev.off()
    rm(grid_eu_y1)
}
cat("compared each year to following year \n")

td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("niche comparisons finished", td, "secs", "\n")

# calculate niche dynamic index
# di = ecospat.niche.dyn.index(grid_as, grid_eu, intersection = 0.1)$dynamic.index.w
# expansion stability unfilling
# 0.1521472 0.8478528 0.2763685
