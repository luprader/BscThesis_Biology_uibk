### up to 04/08/2023
This entry compiles everything up to the first entry.  
The Gbif data for *Harmonia axyridis* was downloaded from the website, Copernicus land cover classification layers for 1992-2020 as well. 
The CHELSA bioclim layers were downloaded with wget.  
At first, all Gbif occurrences were checked for NA values in coordinates, year and coordinate uncertainty. 
Afterwards, they were cleaned using all default tests for the function 'clean_coordinates' from 'CoordinateCleaner'.  
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
"Cleaned EU < 2002: 2 Cleaned Native < 2002: 12"  
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
Subsetting was not necessary, since the generated points were automatically assigned by the 'buffer' and 'spatSample' functions from 'terra'.
The script currently takes a lot of time due to the minimum distance check as it is implemented right now.
There might be better options of doing this, for example with intersected buffer polygons instead of a normal distance calculation.

### 09/08/2023
The generation of absence/background points was improved to have a feasible runtime of <3h.
(the past script was never executed to completion, took way too long)
The point distance approach was removed and replaced with removing distance buffer circles from the range circles and generating points in these new polygons.
The distance circles are also merged into one polygon to prevent unnecessary iterations when erasing.
The current amount of absences per presence is 5, with a range of 10 km and a minimum distance of 1 km to other presences.