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

# specify directories:
# out directory
out_dir <- 'output/antivenom_coverage/'
# cat 1 rasters
cat_1_ras <- 'output/species_rasters/cat_1/'
# cat 2 rasters
cat_2_ras <- 'output/species_rasters/cat_2/'
# combined rasters
combined_ras <- 'output/species_rasters/combined_cat/'

# define category ids
c1_id <- c('1', '1&2')
c2_id <- c('2', '1&2')

# get info on the snake species covered by each antivenom
antivenom_products <- as.list(names(venom))

# drop snakes in the list for which we only have species level data
snakes_to_drop <- c('Acanthophis_spp', 'Bothrops_spp', 'Crotalus_spp', 'Micrurus_spp')
snakes <- snakes[!snakes$split_spp %in% snakes_to_drop, ]

# loop through each antivenom product, and generate a raster showing the effective
# coverage of each product, based on the distribution of the snakes spp. that it is
# effective for
for(i in 1:length(antivenom_products)){
  
  # time process
  message(paste('Generating antivenom coverage surface', i, 'of', length(antivenom_products), Sys.time(), sep = ' '))
  
  # get antivenom product index
  # 1. get product name
  product_name <- antivenom_products[i]
  
  # 2. get list of snake spp. identification codes for that product
  effective_snakes <- venom[, colnames(venom) %in% product_name]
  
  # 3. grab paths for spp. rasters
  snake_idx <- which(snakes$id %in% effective_snakes)
  snake_sub <- snakes[snake_idx, ]
  # remove ones without a raster
  snake_sub <- snake_sub[!is.na(snake_sub$raster_path), ]
  
  # split by category
  cat_1_sub <- snake_sub[snake_sub$processing_cat %in% c1_id, ]
  cat_2_sub <- snake_sub[snake_sub$processing_cat %in% c2_id, ]
  
  # list individual species
  c1_list <- unique(cat_1_sub$split_spp, na.rm = TRUE)
  c2_list <- unique(cat_2_sub$split_spp, na.rm = TRUE)
  combined_list <- unique(snake_sub$split_spp, na.rm = TRUE)
  
  # initialize the cluster (50 cores)
  registerDoMC(50)
  
  # 4. for each medical class, loop through species, open rasters and stack
  if(length(c1_list) > 0){
    
    c1_antivenom_surface <- foreach(i = 1:length(c1_list)) %dopar% {
      
      # get path to each spp.
      spp <- as.character(c1_list[i])
      path <- paste0(cat_1_ras, spp, '.tif')
      
      # read in raster
      snake_rast <- raster(path)
      
      # extend the extent to match the global raster
      snake_ext <- extend(snake_rast, global, value = NA)
      
      return(snake_ext)
      
    }
    
    # stack the rasters
    c1_stack <- stack(c1_antivenom_surface)
    
  }
  
  if(length(c2_list) > 0){
    
    c2_antivenom_surface <- foreach(i = 1:length(c2_list)) %dopar% {
      
      # get path to each spp.
      spp <- as.character(c2_list[i])
      path <- paste0(cat_2_ras, spp, '.tif')
      
      # read in raster
      snake_rast <- raster(path)
      
      # extend the extent to match the global raster
      snake_ext <- extend(snake_rast, global, value = NA)
      
      return(snake_ext)
      
    }
    
    # stack the rasters
    c2_stack <- stack(c2_antivenom_surface)
    
  }
  
  if(length(combined_list) > 0){
    
    comb_antiv_surface <- foreach(i = 1:length(combined_list)) %dopar% {
      
      # get path to each spp.
      spp <- as.character(combined_list[i])
      path <- paste0(combined_ras, spp, '.tif')
      
      # read in raster
      snake_rast <- raster(path)
      
      # extend the extent to match the global raster
      snake_ext <- extend(snake_rast, global, value = NA)
      
      return(snake_ext)
      
    }
    
    # stack the rasters
    combined_stack <- stack(comb_antiv_surface)
    
  }
  
  # 5. get sum of values for each pixel in stack (basically combine all rasters and 
  # generate a value of number of snake spp. per pixel covered by each effective antivenom)
  loop_list <- as.list(c('c1_stack', 'c2_stack', 'combined_stack'))
  
  # do this in parallel (if possible - may be too computationally heavy)
  gen_surface <- foreach(i = 1:length(loop_list)) %dopar% {
    
    # check that the stack exists
    if(length(loop_list[[i]] > 0)) {
      
    antivenom_coverage <- sum(loop_list[[i]], na.rm = TRUE)
    
    # generate antivenom outpath
    outpath <- paste0(out_dir, 
                      product_name,
                      '_',
                      loop_list[[i]])
    
    # 6. write out combined output
    writeRaster(antivenom_coverage,
                file = outpath,
                format = 'GTiff',
                overwrite = TRUE)
    
    complete <- paste(product_name, loop_list[[i]], 'complete = TRUE', sep = " ")
    
    return(complete)
    
    } else{
    
    complete <- paste(product_name, loop_list[[i]], 'complete = FALSE', sep = " ")  
  
    }
  
    }  
  
  message(paste('Completed generating antivenom coverage surface', i, 'of', length(antivenom_products), Sys.time(), sep = ' '))
  
}

