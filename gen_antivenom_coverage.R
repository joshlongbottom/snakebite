# generate antivenom efficacy layers
# clear workspace
rm(list = ls())

# load required packages
require(raster)

# read in antivenom csv
venom <- read.csv("Z:/users/joshua/Snakebite/snakebite/antivenom.csv",
                  stringsAsFactors = FALSE, 
                  na.strings = c("NA","NaN", "", " "))

# read in associated raster info dataset
snakes <- read.csv("Z:/users/joshua/Snakebite/snakebite/snake_list.csv",
                   stringsAsFactors = FALSE,
                   na.strings = c("NA","NaN", "", " "))

# read in the global extent raster
global <- raster("Z:/users/joshua/Snakebite/raster_mask/CoastGlobal_5k.tif")

# get info on the snake species covered by each antivenom
antivenom_products <- as.list(names(venom))

# loop through each antivenom product, and generate a raster showing the effective
# coverage of each product, based on the distribution of the snakes spp. that it is
# effective for

for (i in 1:length(antivenom_products)){
  
  # get antivenom product index
  # 1. get product name
  product_name <- antivenom_products[i]
  
  # 2. get list of snake spp. identification codes for that product
  effective_snakes <- venom[, colnames(venom) %in% product_name]
  
  # 3. grab paths for spp. rasters
  snake_idx <- which(snakes$id %in% effective_snakes)
  snake_sub <- snakes[snake_idx, ]
  list <- unique(snake_sub$raster_path, na.rm = TRUE)
  list <- list[!is.na(list)]
  
  # 4. loop through, open raster, and stack
  for (i in 1:length(list)){
  
  # get path to each spp.
  path <- as.character(list[i])
  
  # read in raster
  snake_rast <- raster(path)
  
  # extend the extent to match the global raster
  snake_ext <- extend(snake_rast, global, value = NA)
  
  # stack the rasters
  if(i == 1){
    
    snake_stack <- snake_ext
    
  } else {
    
    snake_stack <- stack(snake_stack, snake_ext)  
    
  }
  
  }

  # 5. get sum of values for each pixel in stack (basically combine all rasters
  # and generate a value of number of snake spp. per pixel covered by each 
  # effective antivenom)
  antivenom_coverage <- sum(snake_stack, na.rm = TRUE)

  # generate antivenom outpath
  outpath <- paste('Z:/users/joshua/Snakebite/output/antivenom_coverage/', 
                   product_name,
                   sep = '')

  # 6. write out combined output
  writeRaster(antivenom_coverage,
              file = outpath,
              format = 'GTiff',
              overwrite = TRUE)

}

