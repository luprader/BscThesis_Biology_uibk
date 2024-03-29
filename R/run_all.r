# this file executes all scripts to generate the final data used for this project
#all_time = Sys.time()

# crop land cover layers and aggregate to bioclim resolution
source("R/0.1-lc_crop-aggr.r", encoding = "UTF-8")

# merge bioclim layers and crop/rescale to land cover extent
source("R/0.2-clim_merge-crop-res.r", encoding = "UTF-8")

# clean occurrences to desired data quality
source("R/0.3-occ_clean.r", encoding = "UTF-8")

# generate pseudo-absences for modelling
source("R/1.1-gen_absences.r", encoding = "UTF-8")

# extract raster values of the pa-data for modelling
source("R/1.2-pa_extract_vals.r", encoding = "UTF-8")

# select variables for modelling
source("R/1.3-model_var_select.r", encoding = "UTF-8")

# compare niches of the species in different areas/at different times
source("R/2.1-niche_comp.r", encoding = "UTF-8")

# build models
source("R/2.2-model_building.r", encoding = "UTF-8")

# analyse modelling results (change to .r script)
source("R/2.3-analysis.rmd", encoding = "UTF-8")

#td <- difftime(Sys.time(), all_time, units = "secs")[[1]]
#cat("building and model evaluation finished:", td, "secs", "\n")