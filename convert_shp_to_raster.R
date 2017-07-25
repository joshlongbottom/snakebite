# clear workspace
rm(list = ls())

# load required libraries
pacman::p_load(raster, foreach, doMC)

# list all files in the 'modified_ranges' folder
# this folder should contain the 'updated range' for a species; if no range modification has occurred, 
# this folder will contain the original digitised EOR shapefile
shape_list <- list.files('output/all_ranges_incl_modified',
                         pattern = ".shp$",
                         full.names = TRUE)

# read in a 5 km x 5 km raster to use as a template to rasterize by
raster_mask <- raster('data/raw/raster/land_sea/CoastGlobal_5k.tif')

# initialize the cluster
registerDoMC(50)

# loop through and convert shape files 
raster_range <- foreach(i = 1:length(shape_list)) %dopar% {
  
  # get file name
  file_name <- shape_list[i]
  
  # get snake species name
  spp <- gsub('output/all_ranges_incl_modified/', '', file_name)
  spp <- gsub('.shp$', '', spp)
  
  # inform progress
  message(paste('rasterizing file', i, 'of', length(shape_list), sep = " "))
  
  # open shapefile
  snake_shape <- shapefile(file_name)
  
  # get extents of shapefile
  snake_extent <- extent(snake_shape)
  
  # crop raster mask extent by snake_extent
  raster_crop <- crop(raster_mask, snake_extent)
  
  # create up to three different rasters based on the medical classification of each species
  # generate a med class one shapefile
  class_one_shp <- snake_shape[snake_shape$Med_Class == 1, ]
  
  # generate a med class two shapefile (which is c2 and c3 in instances where c3 exists)
  class_two_shp <- snake_shape[snake_shape$Med_Class >= 2, ]
  
  # convert snake shapefiles to raster
  # convert class one
  if(length(class_one_shp) > 0){
    
    # rasterize to extent 
    class_one_raster <- rasterize(class_one_shp, raster_crop)
    
    # set all raster cells = 1
    class_one_raster[!is.na(class_one_raster)] <- 1
    
    # generate outpath
    class_one_outpath <- paste('output/species_rasters/cat_1/',
                               spp,
                               sep = "")
    
    # save raster
    writeRaster(class_one_raster, 
                file = class_one_outpath,
                format = "GTiff",
                overwrite = TRUE)
    
  }
  
  # convert class two
  if(length(class_two_shp) > 0){
    
    # rasterize to extent 
    class_two_raster <- rasterize(class_two_shp, raster_crop)
    
    # set all raster cells = 1
    class_two_raster[!is.na(class_two_raster)] <- 1
    
    # generate outpath
    class_two_outpath <- paste('output/species_rasters/cat_2/',
                               spp,
                               sep = "")
    
    # save raster
    writeRaster(class_two_raster, 
                file = class_two_outpath,
                format = "GTiff",
                overwrite = TRUE)
    
  }
  
  # generate a combined (class 1 and class 2) raster
  # rasterize to extent
  combined_raster <- rasterize(snake_shape, raster_crop)
  
  # set all raster cells = 1
  combined_raster[!is.na(combined_raster)] <- 1
  
  # generate outpath
  outpath <- paste('output/species_rasters/combined_cat/',
                   spp,
                   sep = "")
  
  # save raster
  writeRaster(combined_raster, 
              file = outpath,
              format = "GTiff",
              overwrite = TRUE)
  
  
  complete <- 'yes'
  return(complete)
  
}
