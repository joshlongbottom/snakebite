# Snakebite project
Code for the generation of a venomous snake species richness surface, identifying populations at risk of exposure to venomous snake species.
Project is currently in progress. Main themes of the project:
- Digitise existing expert opinion range maps for all medically important venomous snake species
- Redefine the ranges for each of these species
- Create a venomous snake species surface
- Compare this surface with disease burden estimates and antivenom availability data

## Description
Expert opinion range (EOR) maps for medically important snake species were obtained from the WHO antivenom database (http://apps.who.int/bloodproducts/snakeantivenoms/database/). **~'get_who_snakes.R'**

Example EOR map:

![Alt text](http://apps.who.int/bloodproducts/snakeantivenoms/database/Images/SnakesDistribution/Small/map_Acanthophis_antarcticus.png "Acanthophis antarcticus EOR map")

Each EOR map was digitized using ArcMAP, and a shapefile was created for each unique species per map.

Example digitization and rasterization process:

![Alt text](https://preview.ibb.co/kZX4Lv/Mapping_venomous_snake_species_richness.png "Digitised Micrurus lemniscatus EOR map")

Utilizing the spocc package, species occurrence data was then obtained for each species (n = ~277) **~'plot_occurrence_spocc.R'**. Through a combination of digitised EOR maps, and species occurrence data obtained via spocc, Multivariate Environmental Suitability Surfaces (MESS) were generated for each species **~'gen_mess.R'**. These surfaces were used to validate any new occurrence data outside of the currently accepted EOR, and informed a process generating a new species range.

Example ocurrence data trawl:

![Alt text](https://preview.ibb.co/hzvVtF/Grab_occurrence.png "Occurrence data grabbing example")

Example MESS surface (in progress, Acanthophis antarcticus):

![Alt text](https://image.ibb.co/j3vr0v/Example_MESS.png "Acanthophis antarcticus MESS (in progress)")

The newly created ranges were then converted to raster format (**~'convert_shp_to_raster.R'**), and stacked (**~'gen_species_richness.R'**) to give a count of the number of venomous snake species per 5 km x 5 km cell.

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


