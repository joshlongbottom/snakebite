# clear workspace
rm(list = ls())

# load required packages
require(raster)
require(rgdal)
require(maptools)

# read in shapefiles to merge 
shape1 <- shapefile("Z:/users/joshua/Snakebite/WHO EOR shapefiles/Agkistrodon_contortrix_cat1.shp")
shape2 <- shapefile("Z:/users/joshua/Snakebite/WHO EOR shapefiles/Agkistrodon_contortrix_cat2.shp")

# bind the shapefiles together
shape_bound <- union(shape1,
                     shape2)

# write out the bound shapefile
writeSpatialShape(shape_bound, 
                  "Z:/users/joshua/Snakebite/WHO EOR shapefiles/Agkistrodon_contortrix.shp")
