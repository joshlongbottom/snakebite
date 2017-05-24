# Snakebite project
Code for the generation of a venomous snake species richness surface, identifying populations at risk of exposure to venomous snake species

## Dependancies
Packages:
1. dismo
2. raster
3. maptools
4. rgeos
5. rgdal
6. pdftools
7. foreign
8. spatstat
9. utils
10. spocc
11. seegSDM
12. ggplot2
13. reshape2

## Description
Expert opinion range (EOR) maps for medically important snake species were obtained from the WHO antivenom database (http://apps.who.int/bloodproducts/snakeantivenoms/database/). ~'get_who_snakes.R'

Each EOR map was digitized using ArcMAP, and a shapefile was created for each unique species per map. 

These digitized ranges were then converted to raster format (~convert_shp_to_raster.R), and stacked (~gen_species_richness.R) to give a count of the number of venomous snake species per 5 km x 5 km cell.

Utilizing the spocc package, species occurrence data was obtained for each species. Using this species occurrence data, multivariate environmental suitability surfaces were generated for each species, in order to validate the existing EOR map.
