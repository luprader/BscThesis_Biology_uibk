### up to 04/08/2023
This entry compiles everything up to the first entry.  
The Gbif data for *Harmonia axyridis* was downloaded from the website, Copernicus land cover classification layers for 1992-2020 as well. 
The CHELSA bioclim layers were downloaded with wget.  
At first, all Gbif occurences were checked for NA values in coordinates, year and coordinate uncertainty. 
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
The occurence cleaning code was rewritten to first subset the occurences into the two spatial extents before conducting tests, since for example outlier tests would be heavily influenced by the global distribution. The tests used now are :
('capitals', 'centroids', 'duplicates', 'equal', 'institutions', 'outliers', 'seas', 'zeros'), with the default settings for each test.
Using the 'urban' test resulted in very large amounts of flags, which makes sense given that the distribution of *Harmonia axyridis* is in part due to human distribution as a control agent.
To keep some of that information and still correct for some bias, 'urban' was not used as a test, but the 'capitals' test with a rather large radius of 10 km.
The subsetting before cleaning targets NAs in coordinates, year, uncertainty and also removes all points after 2022, before 1991 in europe (first invasive occurence) and with an uncertainty larger than 1 km.

