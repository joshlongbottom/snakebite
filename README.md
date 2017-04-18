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

## Description
Expert opinion range (EOR) maps for medically important snake species were obtained from the WHO antivenom database (http://apps.who.int/bloodproducts/snakeantivenoms/database/). ~'get_who_snakes'

Each EOR map was digitized using ArcMAP, and a shapefile was created for each unique species per map. 

stack these maps to give a count of the
number of venomous snake species per 5 km x 5 km cell.
