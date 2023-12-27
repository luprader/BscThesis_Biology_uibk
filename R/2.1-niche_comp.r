# this file computes and visualizes niche comparisons according to
# (Broennimann et al. 2011), using the `ecospat` library

# used libraries
library(ecospat)
library(ade4)
library(dplyr)

tot_time <- Sys.time()
# load extracted pa data
pa <- readRDS("R/data/occurrence_data/axyridis_pa_vals_extracted.rds")
pa_as <- subset(pa, Area == "as")
po_as <- subset(pa_as, Presence == "present")
pa_eu <- subset(pa, Area == "eu")
po_eu <- subset(pa_eu, Presence == "present")
bio_tot <- select(pa, starts_with("bio")) # only use bioclim vars

## compare total native to total invaded
# pca calibrated on all datapoints
pca_tot <- dudi.pca(bio_tot, scannf = FALSE, nf = 2)
# plot pca axis contribution
fname <- "R/plots/niche_comp/as_eu_pca.png"
png(width = 600, height = 600, filename = fname)
ecospat.plot.contrib(contrib = pca_tot$co, eigen = pca_tot$eig)
dev.off()
# pca scores for the whole area
scores_tot <- pca_tot$li
# pca scores for native occs
scores_po_as <- suprow(pca_tot, select(po_as, starts_with("bio")))$li
# pca scores for native occs with background
scores_pa_as <- suprow(pca_tot, select(pa_as, starts_with("bio")))$li
# pca scores for eu occs
scores_po_eu <- suprow(pca_tot, select(po_eu, starts_with("bio")))$li
# pca scores for eu occs with background
scores_pa_eu <- suprow(pca_tot, select(pa_eu, starts_with("bio")))$li

# calculate native dynamic occurrence densities grids
grid_as <- ecospat.grid.clim.dyn(
    glob = scores_tot,
    glob1 = scores_pa_as, sp = scores_po_as, R = 100, th.sp = 0
)
grid_eu <- ecospat.grid.clim.dyn(
    glob = scores_tot,
    glob1 = scores_pa_eu, sp = scores_po_eu, R = 100, th.sp = 0
)

# compute niche overlap (Schoeners overlap metric)
ol <- round(ecospat.niche.overlap(grid_as, grid_eu, cor = TRUE)$D, digits = 3)
# calculate niche dynamic indices
di <- round(ecospat.niche.dyn.index(grid_as, grid_eu, intersection = 0.1)$dynamic.index.w, digits = 3)

# plot niche overlap
fname <- "R/plots/niche_comp/as_eu_tot_niche.png"
png(width = 600, height = 600, filename = fname)
text <- paste("Niche Overlap native / EU | D:", ol, "\n expansion stability unfilling :", di[[1]], di[[2]], di[[3]])
ecospat.plot.niche.dyn(grid_as, grid_eu,
    quant = 0.1,
    col.unf = "green", col.exp = "red", col.stab = "blue",
    interest = 2, title = text, name.axis1 = "PC1",
    name.axis2 = "PC2"
)
dev.off()

fname <- "R/plots/niche_comp/as_eu_tot_eq-sim.png"
png(width = 1200, height = 600, filename = fname)
par(mfrow = c(1, 2))
# niche equivalency test
eq_test <- ecospat.niche.equivalency.test(grid_as, grid_eu,
    rep = 100, ncores = 2 # , alternative = "greater"
)
# plot eq_test
ecospat.plot.overlap.test(eq_test, "D", "Equivalency Native / EU")
rm(eq_test)
# niche similarity test
sim_test <- ecospat.niche.similarity.test(grid_as, grid_eu,
    rep = 100, rand.type = 2, ncores = 2 # , alternative
)
# plot sim_test
ecospat.plot.overlap.test(sim_test, "D", "Similarity Native / EU")
rm(sim_test)
dev.off()
rm(grid_eu)

cat("compared total EU to native range \n")
rm(list = c("scores_pa_as", "scores_po_as", "scores_pa_eu", "scores_po_eu"))
gc()



## calculate overlap between total native and each EU year
eq_tests <- c()
sim_tests <- c()
overlap <- c()
dynamic <- c()
for (y in 2002:2022) {
    # take subsets
    pa_eu_y <- subset(pa_eu, Year == y)
    po_eu_y <- subset(po_eu, Year == y)
    # pca scores for year
    scores_po_eu_y <- suprow(pca_tot, select(po_eu_y, starts_with("bio")))$li
    scores_pa_eu_y <- suprow(pca_tot, select(pa_eu_y, starts_with("bio")))$li
    # grids for both years
    grid_eu_y <- ecospat.grid.clim.dyn(
        glob = scores_tot,
        glob1 = scores_pa_eu_y, sp = scores_po_eu_y, R = 100, th.sp = 0
    )
    # overlap and niche dynamic
    ol <- round(ecospat.niche.overlap(grid_as, grid_eu_y, cor = TRUE)$D, digits = 3)
    overlap <- rbind(overlap, ol)
    di <- round(ecospat.niche.dyn.index(grid_as, grid_eu_y, intersection = 0.1)$dynamic.index.w, digits = 3)
    dynamic <- rbind(dynamic, di)


    # plot niche overlap
    if (ol == 0) { # correct color for case of no stable overlap
        col_stab <- "#00BA38"
    } else {
        col_stab <- "#619CFF"
    }

    fname <- paste0("R/plots/niche_comp/single_ys/eu_", y, "_niche.png")
    png(width = 600, height = 600, filename = fname)
    text <- paste("Niche Overlap native / EU ", y, " | D:", ol, "\n expansion stability unfilling :", di[[1]], di[[2]], di[[3]])
    ecospat.plot.niche.dyn(grid_as, grid_eu_y,
        quant = 0.1,
        col.unf = "#00ba38", col.exp = "#f8766d", col.stab = "#3f80f1", 
        interest = 2, title = text, name.axis1 = "PC1",
        name.axis2 = "PC2"
    )
    dev.off()

    # niche equivalency test
    eq_test <- ecospat.niche.equivalency.test(grid_as, grid_eu_y,
        rep = 100, ncores = 2 # , alternative = "greater"
    )
    eq_tests <- rbind(eq_tests, eq_test) # save in array

    # niche similarity test
    sim_test <- ecospat.niche.similarity.test(grid_as, grid_eu_y,
        rep = 100, rand.type = 2, ncores = 2 # , alternative
    )
    sim_tests <- rbind(sim_tests, sim_test) # save in array
}
cat("compared each year to native \n")

# save niche sim and eq results
eq_sim <- cbind(eq_tests, sim_tests)
saveRDS(overlap, file = "R/data/modelling/niche_y_overlap.rds")
saveRDS(dynamic, file = "R/data/modelling/niche_y_dynamic.rds")
saveRDS(eq_sim, file = "R/data/modelling/niche_y_eq_sim.rds")

td <- difftime(Sys.time(), tot_time, units = "secs")[[1]]
cat("niche comparisons finished", td, "secs", "\n")
