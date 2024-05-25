# BscThesis_Biology_uibk
### Author: Lukas Prader

## Overview: 
This repository holds the code and files for my bachelor thesis written under supervision of Lauren Talluto, with the title "Studying SDM performance throughout a time series: A case study using the invasive species Harmonia axyridis".

## !!important!!
The code can not be run from the repository itself, since the data required is quite large. 
All files numerated with 0.x run with raw downloaded files from the Gbif, CHELSA or Copernicus databases and require rather large amounts of disk space (~50 GB). 
In total, all generated files together with the raw data take up ~80 GB of space.
The file 'envidats3paths.txt' holds all links to download the CHELSA bioclim layers with wget, Gbif and Copernicus data were downloaded from the respective websites.  
Some functions (i.e. clean_coordinates) need a connection to the internet to download reference maps. 

The code has mainly been written and executed using the IDE VSCode in a WSL environment. 
In general, everything should also be executable in RStudio with the included R project file.
