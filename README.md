# BscBio_LP23
### Author: Lukas Prader

## Overview: 
This repository holds the code and files for my bachelor thesis written in the summer of SS23 under supervision of Matthew Talluto.

## !!important!!
The code can not be run from the repository itself, since the data required is quite large. 
All files numerated with 0.x run with raw downloaded files from the Gbif, CHELSA or Copernicus databases and require rather large amounts of disk space (75 GB). 
In total, all generated files take up X GB of space.
The file 'envidats3paths.txt' holds all links to download the CHELSA bioclim layers with wget, Gbif and Copernicus data were downloaded from the respective websites.  
Some functions (i.e. clean_coordinates) need a connection to the internet to download reference maps. 

The code has mainly been written and executed using the IDE VSCode in a WSL environment. 
In general, everything should also be executable in RStudio with the included R project file.