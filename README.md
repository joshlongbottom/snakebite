# Snakebite project
Code for the generation of a venomous snake species richness surface, identifying populations at risk of exposure to venomous snake species

## Description
Expert opinion range (EOR) maps for medically important snake species were obtained from the WHO antivenom database (http://apps.who.int/bloodproducts/snakeantivenoms/database/). ~'get_who_snakes.R'

Example map:
![Alt text](http://apps.who.int/bloodproducts/snakeantivenoms/database/Images/SnakesDistribution/Small/map_Acanthophis_antarcticus.png "Acanthophis antarcticus EOR map")

Each EOR map was digitized using ArcMAP, and a shapefile was created for each unique species per map.

These digitized ranges were then converted to raster format (~convert_shp_to_raster.R), and stacked (~gen_species_richness.R) to give a count of the number of venomous snake species per 5 km x 5 km cell.
Example digitization and rasterization process:


Utilizing the spocc package, species occurrence data was obtained for each species. Using this species occurrence data, multivariate environmental suitability surfaces were generated for each species, in order to validate the existing EOR map.

## Dependancies
Packages:
1. dismo (https://cran.r-project.org/web/packages/dismo/index.html)
2. raster (https://cran.r-project.org/web/packages/raster/index.html)
3. maptools (https://cran.r-project.org/web/packages/maptools/index.html)
4. rgeos (https://cran.r-project.org/web/packages/rgeos/index.html)
5. rgdal (https://cran.r-project.org/web/packages/rgdal/index.html)
6. pdftools (https://cran.r-project.org/web/packages/pdftools/index.html)
7. foreign (https://cran.r-project.org/web/packages/foreign/index.html)
8. spatstat (https://cran.r-project.org/web/packages/spatstat/index.html)
9. utils (https://cran.r-project.org/web/packages/R.utils/index.html)
10. spocc (https://cran.r-project.org/web/packages/spocc/index.html)
11. seegSDM (https://github.com/SEEG-Oxford/seegSDM)
12. ggplot2 (https://cran.r-project.org/web/packages/ggplot2/index.html)
13. reshape2 (https://cran.r-project.org/web/packages/reshape2/index.html)


