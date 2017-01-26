# clear workspace
rm(list = ls())

# load required libraries
require(raster)

# list all extracted shapefiles
list <- as.list(list.files("Z:/users/joshua/Snakebite/WHO EOR shapefiles/", 
                           pattern = '*shp$'))

# loop through and convert shape files 
for(i in 1:length(list)){
  
  # get file name
  file_name <- unlist(list[i])
  
  # get snake species name
  spp <- gsub('.shp$', '', file_name)
  
  # gen path to shapefile
  paths <- paste("Z:/users/joshua/Snakebite/WHO EOR shapefiles/",
                 file_name,
                 sep = "")
  
  # read in a 5 km x 5 km raster to use as a template to rasterize by
  raster_mask <- raster("Z:/users/joshua/Snakebite/raster_mask/CoastGlobal_5k.tif")
  
  # inform progress
  message(paste('rasterizing file', i, 'of', length(list), sep = " "))
  
  # open shapefile
  snake_shape <- shapefile(paths)
  
  # get extents of shapefile
  snake_extent <- extent(snake_shape)
  
  # crop raster mask extent by snake_extent
  raster_mask <- crop(raster_mask, snake_extent)
  
  # convert snake shapefile to raster
  snake_raster <- rasterize(snake_shape, raster_mask)
  
  # set all raster cells = 1
  snake_raster[!is.na(snake_raster)] <- 1
  
  # generate outpath
  outpath <- paste("Z:/users/joshua/Snakebite/WHO EOR rasters/",
                   spp,
                   sep = "")
  
  # save raster
  writeRaster(snake_raster, 
              file = outpath,
              format = "GTiff",
              overwrite = TRUE)
  
}