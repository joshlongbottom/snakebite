# clear workspace
rm(list = ls())

# load required packages
require(raster)

# list all rasters in the directory
list <- as.list(list.files("Z:/users/joshua/Snakebite/WHO EOR rasters/", 
                           pattern = '*tif$',
                           full.names = TRUE))

# read in the global extent raster
global <- raster("Z:/users/joshua/Snakebite/raster_mask/CoastGlobal_5k.tif")

# loop through, open raster, and stack
for (i in 1:length(list)){
  
  # get path to each spp.
  path <- as.character(list[i])
  
  # read in raster
  snake_rast <- raster(path)
  
  # extend the extent to match the global raster
  snake_ext <- extend(snake_rast, global, value = NA)
  
  # inform progress
  message(paste('stacking raster', i, 'of', length(list), sep = " "))
  
  # stack the rasters
  if(i == 1){
  
    snake_stack <- snake_ext
  
    } else {
    
    snake_stack <- stack(snake_stack, snake_ext)  
    
  }
  
}

# get sum of values for each pixel in stack (basically combine all rasters
# and generate a value of number of snake spp. per pixel)
richness <- sum(snake_stack, na.rm = TRUE)

# check it's worked
plot(richness)

# write out combined output
writeRaster(richness,
            file = 'Z:/users/joshua/Snakebite/output/species_richness',
            format = 'GTiff',
            overwrite = TRUE)

