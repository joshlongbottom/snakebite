# generate antivenom efficacy layers
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreach, doMC)

# read in antivenom csv
venom <- read.csv('data/raw/antivenom.csv',
                  stringsAsFactors = FALSE, 
                  na.strings = c("NA","NaN", "", " "))

# read in associated raster info dataset
snakes <- read.csv('data/raw/snake_list_cluster.csv',
                   stringsAsFactors = FALSE)

# read in the global extent raster
global <- raster('data/raw/raster/land_sea/CoastGlobal_5k.tif')

# specify out directory
out_dir <- 'output/antivenom_coverage'

# get info on the snake species covered by each antivenom
antivenom_products <- as.list(names(venom))

# initialize the cluster (50 cores)
registerDoMC(50)

# time process
message(paste('Generating antivenom coverage surfaces', Sys.time(), sep = ' '))

# loop through each antivenom product, and generate a raster showing the effective
# coverage of each product, based on the distribution of the snakes spp. that it is
# effective for
antivenom_surface <- foreach(i = 1:length(antivenom_products)) %dopar% {
  
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
  outpath <- paste0(out_dir, 
                    product_name)

  # 6. write out combined output
  writeRaster(antivenom_coverage,
              file = outpath,
              format = 'GTiff',
              overwrite = TRUE)
  
  complete <- paste0(product_name, ' complete')
  
  return(complete)

}

message(paste('Completed generating antivenom coverage surfaces', Sys.time(), sep = ' '))
