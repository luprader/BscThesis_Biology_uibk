### up to 04/08/2023
This entry compiles everything up to the first entry.  
The Gbif data for *Harmonia axyridis* was downloaded from the website, Copernicus land cover classification layers for 1992-2020 as well. 
The CHELSA bioclim layers were downloaded with wget.  
At first, all Gbif occurrences were checked for NA values in coordinates, year and coordinate uncertainty. 
Afterwards, they were cleaned using all default tests for the function `clean_coordinates` from `CoordinateCleaner`.  
The land cover layers were aggregated to match the resolution of the climate data, so a 3x3 cell area of the original layer was aggregated to one larger cell, matching the resolution of the CHELSA bioclim layers.
The new value was the mode of the 9 grid cells.  
The global layer was cropped to two spatial extents. 
The europe extent was arbitrarily chosen: 
(xmin = -25, ymin = 35, xmax = 65, ymax = 72). 
The native extent was chosen from (Orlova-Bienkowskaja, Ukrainsky & Brown, 2015), using Fig. 1 and the mentioned locations to draw an outline around the native range:
(xmin = 70, ymin = 20, xmax = 150, ymax = 65).  
The CHELSA bioclim layers were merged, cropped and resampled to match the two new extents of the land cover layers.  
A rough introduction and methods draft was written, as well as this entry.


### 04/08/2023
The occurrence cleaning code was rewritten to first subset the occurrences into the two spatial extents before conducting tests, since for example outlier tests would be heavily influenced by the global distribution. The tests used now are :
('capitals', 'centroids', 'duplicates', 'equal', 'institutions', 'outliers', 'seas', 'zeros'), with the default settings for each test.
Using the 'urban' test resulted in very large amounts of flags, which makes sense given that the distribution of *Harmonia axyridis* is in part due to human distribution as a control agent.
To keep some of that information and still correct for some bias, 'urban' was not used as a test, but the 'capitals' test with a rather large radius of 10 km.
The subsetting before cleaning targets NAs in coordinates, year, uncertainty and also removes all points after 2022, before 1991 in Europe (first invasive occurrence) and with an uncertainty larger than 1 km.

### 05/08/2023
Some analysis was conducted to visualize the amount of occurrences per year for each spatial extent. 
It was noticed, that the amount of European occurrences prior to 2000 was greatly underestimated. 
The timeframe will probably only be 2002-2022 instead of 1992-2022 because of that, but a more detailed analysis will be done tomorrow.

### 06/08/2023
The amounts of occurrences per year were computed in more detail.
This resulted in two plots (R/plots/clean_occ_yearly_*.png) and the following values:
"Raw EU < 2002: 32 Raw Native < 2002: 35"  
"Cleaned EU < 2002: 0 Cleaned Native < 2002: 5"  
Because of this, even though the first occurrence according to EASIN was 1991 in Belgium, the temporal extent will not include years prior to 2002 anymore.
There might be value in creating an initial model with the 2 and 12 cleaned points prior to 2002, so for now they are not excluded in the cleaning process.
Only the land cover layers that were not needed anymore were removed from the used data folder.
The thought has come to mind, that it might be worthwhile to not use a general duplicate test, but instead only check for duplicates in one year.
Having the exact same coordinates in several years would indicate lazy sampling though and might not be desirable because of that.

The outline intro was refined more and the methods adapted to include the new temporal extent.

The generation of pseudo-absence/background points was started, generating 5 points in a 10 km radius around a presence point and checking for water or a distance closer than 0.5km to any other presence point. 
The presence points will be subset per year and area to be able to give them the corresponding labels.

### 08/08/2023
The generation of background points was finished.
Subsetting was not necessary, since the generated points were automatically assigned by the `buffer` and `spatSample` functions from `terra`.
The script currently takes a lot of time due to the minimum distance check as it is implemented right now.
There might be better options of doing this, for example with intersected buffer polygons instead of a normal distance calculation.

### 09/08/2023
The generation of absence/background points was improved to have an estimated runtime of 3h.
(the past script was never executed to completion, took way too long)
The point distance approach was removed and replaced with removing distance buffer circles from the range circles and generating points in these new polygons.
The distance circles are also merged into one polygon to prevent unnecessary iterations when erasing.
The current amount of absences per presence is 5, with a range of 10 km and a minimum distance of 1 km to other presences.

### 10/08/2023
Merging all minimum distance circles (1 km), failed since the combined polygon was too large to compute. 
Combining seems to not improve the runtime at all after further investigation.
The new approach now is to geographically subset all points and compute the absences for each subset.
The total extent is now split into extents until the desired density per extent is reached. (function `lp_subdiv_pts`)
Afterwards, all occurrences are cropped to those extents and then separately used for absence generation. (function `lp_gen_abs`)
The functions for these steps are defined in a separate function.r script.


### 11/08/2023
A benchmark of the absence generation regarding different desired point densities (end_ptcount) then computing sub extents was attempted.
After creating a dummy land cover layer with some NA space to include absence flagging, different subset amounts relative to the tested total point count were tested: none, 1/2, 1/10 and 1/100. 
This means that for 1000 presence points, a subset to 100, 10 and 500 was tested for computation speed.
This process was repeated for 5000 and 10000 dummy presences, resulting in the plot 'abs_gen_subdiv_bench.png'.
From this, the decision was made to run absence generation with a desired density of <10000 presences, meaning a division of ~1/12.6 for ~126000 total presences.

### 12/08/2023
An error in the absence replacement generation was fixed and absence generation was run for the first time for all presences.
The runtime was 43 minutes.

I want to take the time here to explain the subsetting in more detail, also touching on drawbacks:
Starting from the merged extent of Europe and Asia, the extent is split along its longer side.
The resulting two extents are added to a list of extents to check.
When the loop selects an extent, it calculates the amount of points inside it.
If the amount is larger than the specified cutoff (end_ptcount), the extent is split again.
Otherwise, if the count is smaller or equal, the extent is moved to a list of "good" extents.
A count of 0 results in the extent being removed from all lists.

This approach is pretty fast and computation efficient, but it comes with a drawback.
When generating absences, the main goal is to not generate absences too close to presences.
A presence point could be very close to an extent border though.
This introduces the possibility of absences being generated on the other side of the border, which should normally have been excluded for being too close.  
One could improve the subsetting process to take a minimum distance to the border into account, but this would probably result in way longer computation times, possibly not improving the total time anymore.
At this point in time, this drawback is accepted as part of the absence generation, not only because the current amount of extents (29), does not introduce too many borders.  
One has to also acknowledge that the chance for this error will be severely higher in regions with more point density, which have to also be split more often to reach the max point count.
In those areas, the error will not have a severe impact though, since the presence density is already so high.

### 13/08/2023
A first version of bioclim and land cover value extraction for the pa data was finished.
It will be rewritten with a defined function and progress statements tomorrow.
There might be benefit in running the absence generation with no subsets, since the estimated runtime would only bee double the time with subsets according to the benchmark results (so <2h).
This would remove any possible errors, but there might be issues with memory similar to a merged polygon from 10/08/2023. 
A test run will be done overnight.

### 14/08/2023
The test run without sub extents failed with an out of memory error.
Maybe there are options to fix the memory error, with detriment to runtime, but this will not be pursued for now.  
The value extraction for all years was rewritten with a defined function `lp_ext_vals`.

The maximum distance for absences was adjusted to 18 km, the estimated "typical" flight distance for *Harmonia axyridis* according to (Jeffries et al. 2013).
There are many papers discussing minimum and maximum sampling distance of absences.
The conclusions are often rather specific, so it is hard to know what the best values would be right now.
With this, the current choices are rather arbitrary, with 1 km minimum being the maximum presence location uncertainty, and typical flight distance as a maximum just to have some limit on the potential area of interest for a specific presence point.

After checking the pa data, an error with the absence generation was found.
Generated absences were never checked against water or NA.
The code was fixed and new benchmarks were run.

### 15/08/2023
The benchmark results show a different optimum of subdivision now, but not a huge increase in runtime as feared at first.
Memory limitations seem to occur with subdivisions above 10000, this will be investigated. (maybe WSL issue)
For now, some minor memory management was added to the absence generation process.

### 16/08/2023
It was found that generated absences and some presences were still in water according to their lccs_class.
This was in part due to not testing for land cover in the cleaning process.
The occurrence cleaning process was then extended to test for water or NA in the corresponding Copernicus land cover layers.
The absence generation had to be modified to use the correct land cover layer for each year, while still testing for proximity to presences from all years.
Subdivision is now only used for Europe, since Asia has so little points in comparison.
A new benchmark was conducted, with a similar optimum to yesterday (1/2).
A full computation of absences will be conducted overnight.

## 17/08/2023
Absence generation was successful, with a computation time of ~92 mins for a Europe subdivision to <20000.
The absence generation is now only conducted for years 2002-2020.
Those few points prior to 2002 are negligible, and 21/22 will only be used for model validation.
With this, a new plot of the used subdivisions was created and the procedure of variable selection was started.

There is the question of what data to use for variable selection, since variable correlation can be different for different years for example.
This was also visualized comparing correlation matrices for all bioclim variables of 2002, 2011 and 2020.  
Since the iteration should help with future models of emerging invasive species, an approach using the earliest usable year seems reasonable.
This way, the iteration can model an example of really trying to model *H. axyridis* at the time of its emergence.

Right now, bioclim variable selection uses variance inflation factors of a GLM fitted on data only from 2002 to choose the best variables for modeling. Variables, also using the squared values of each variable as an option, are dropped iteratively until no VIF in the model exceeds a value of 10.
Squared variables are also dropped before their linear counterparts.
The order of dropping seems to have a large influence on the outcome though.
For comparison, this iteration was also conducted with only 2020 occurrences used.
The resulting variables for both years are commented out in "1.3-model_var-select.r" right now.

The yearly occurrence plots were also updated with the new cleaned occurrence data.

### 18/08/2023
For variable selection a PCA for the land cover classes was conducted by creating binary columns for each present class.
The PCA results in many dimensions with small explained variances (<15%).
A PCA of the scaled bioclim variables is very successful in comparison, so maybe the method is implemented wrong.
Plots were created to visualize the PCA results for land cover.

### 21/08/2023
Niche comparison was implemented using the `ecospat` library.
It follows methods described in (Broennimann et al. 2011).
Right now, two iterations of niche comparisons are computed.
The first one compares the climatic niche of each year to all data of the native range.
The second iteration compares each year of data in Europe to the year after in order to have a measure for the niche change.
The comparison of each year with the native range might not be that beneficial in the end, since the main interest there would be the difference between the native and invaded range.

Absence generation was extended to also include 2021 and 2022, since the absences are used as background during niche comparison.

### 22/08/2023
Model variable selection was extended to now incorporate the land cover PCA results into the VIF selection process.
Land cover PCA was fixed by correcting `scale.unit` to FALSE, since a 0/1 binary variable is obviously not scaled to zero mean and unit variance.
This results in a great PCA with almost 40% of variance covered by PC1.
The VIFs of each PCA component are very low, even at the beginning of dropping variables.
With each run until now, no land cover dimension has been dropped by the VIF algorithm.

It was also discovered, that a large portion of occurrences are in areas with the lccs_class 'urban'.
This could have bad influences on the modelling capability, yet it could also be a distinct feature of *H. axyridis*.
The modelling will show if those points skew the predictions too much.

A new file named openissues.md was created to collect all concerns and unresolved issues regarding the project.
It was also verified, that the buffer generation creates equal area circles.

### 23/08/2023
The land cover PCA used in variable selection was improved to use a cutoff value of 80% for the cumulative variance covered by PCA axes.
To prepare for model building and prediction, the selected variables were added as value columns to the pa data, reprojecting land cover to the used PCA dimensions (`lp_pca_proj`).

To enable prediction in geographic space, rasters were created with layers corresponding to the projected lccs classes (`lp_pca_proj_lc`).
It uses a parallel for loop for the years of land cover data, since the projection is very CPU intensive.

## 24/08/2023
After it was attempted to run 1.1-gen_absences overnight, an error with the generation in Europe was fixed (21/22 not using sub extent).
The iteration over each sub extent in Europe was also parallelized with `foreach` to improve computation time.

With this finished, model building was started.

### 25/08/2023
A full absence generation was run, having a runtime of ~130 min.
The outline was improved, updating and expanding the methods to the current state of progress.

### 29/08/2023
It was noticed that a glm fit on the native data did not produce a working model at all.
Investigating the value histograms it was noticed that absence generation created absences mimicking the distribution of the presence points.
For now, absence generation was simplified to just take random spatial samples in the specified extent for each year.

### 30/08/2023
For all generated raster files, a format switch from .grd to .tif was made, since they take up less space and are used by `terra` for internal processes anyway.

### 31/08/2023
The function `lp_pca_proj_lc` was modified to remove oceans beforehand, which speeds up computation.
The generation of prediction rasters was removed for now, since model accuracy can be computed just from the extracted pa data as well.
With this, the creation of a modelling-ready dataframe was moved from `1.4-model_data_prep.r` to `1.3-model_var_select.r`.
With randomly sampled absences it was noticed that the bio14 layer for 2011-2040 has NA cells at one edge, which is why the geographic extent for Europe was slightly modified.

A first version of the basic model building structure was made, implementing a native fit with `glm` and `gam`.

### 01/09/2023
Model evaluation was implemented, as well as fits with `gbm` and `maxnet`.
Right now, a complete creation of a native model and evaluation with 2022 data is realized in `2.2-model_building.r`.
With this, the basic pipeline for this project is more or less in place.

### 03/09/2023
Model evaluation was changed to return a list of accuracy measurements (PCC, sensitivity, specificity and Cohens Kappa).
A test run with 2004 and 2008 as cutoff years showed that the current maxent implementation has a very long computation time, the used formula might need to be simplified.