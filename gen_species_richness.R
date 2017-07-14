# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreach)

# create lists of all rasters in each directory
c1_list <- as.list(list.files('output/species_rasters/cat_1', 
                              pattern = '*tif$',
                              full.names = TRUE))

c2_list <- as.list(list.files('output/species_rasters/cat_2', 
                              pattern = '*tif$',
                              full.names = TRUE))

combined_list <- as.list(list.files('output/species_rasters/combined_cat', 
                                    pattern = '*tif$',
                                    full.names = TRUE))


# read in the global extent raster
global <- raster('data/raw/raster/land_sea/CoastGlobal_5k.tif')

# loop through, open each raster and extend to the global extent
message(paste('Generating category 1 species richness surface ', Sys.time(), sep = ''))

# start with the C1 species
c1_rasters <- foreach(i = length(c1_list)) %do% {
  
  # get path to each spp.
  path <- as.character(c1_list[i])
  
  # read in raster
  snake_rast <- raster(path)
  
  # extend the extent to match the global raster
  snake_ext <- extend(snake_rast, global, value = NA)
  
  return(snake_ext)
  
}

# stack the cat 1 rasters
c1_stack <- stack(c1_rasters)  

# get sum of values for each pixel in stack (basically combine all rasters
# and generate a value of number of snake spp. per pixel)
c1_richness <- sum(c1_stack, na.rm = TRUE)

# write out the category 1 species richness output
c1_outpath <- paste('output/species_richness/', 'category_1_', Sys.Date(), sep = '')
writeRaster(c1_richness,
            file = c1_outpath,
            format = 'GTiff',
            overwrite = TRUE)

# now stack the C2 species
message(paste('Generating category 2 species richness surface ', Sys.time(), sep = ''))

c2_rasters <- foreach(i = length(c2_list)) %do% {
  
  # get path to each spp.
  path <- as.character(c2_list[i])
  
  # read in raster
  snake_rast <- raster(path)
  
  # extend the extent to match the global raster
  snake_ext <- extend(snake_rast, global, value = NA)
  
  return(snake_ext)
  
}

# stack the cat 1 rasters
c2_stack <- stack(c2_rasters)  

# get sum of values for each pixel in stack (basically combine all rasters
# and generate a value of number of snake spp. per pixel)
c2_richness <- sum(c2_stack, na.rm = TRUE)

# write out the category 1 species richness output
c2_outpath <- paste('output/species_richness/', 'category_2_', Sys.Date(), sep = '')
writeRaster(c2_richness,
            file = c2_outpath,
            format = 'GTiff',
            overwrite = TRUE)

# now generate the combined (c1 & c2 species) stack
message(paste('Generating combined species richness surface ', Sys.time(), sep = ''))

combined_class_rasters <- foreach(i = length(combined_list)) %do% {
  
  # get path to each spp.
  path <- as.character(combined_list[i])
  
  # read in raster
  snake_rast <- raster(path)
  
  # extend the extent to match the global raster
  snake_ext <- extend(snake_rast, global, value = NA)
  
  return(snake_ext)
  
}

# stack the cat 1 rasters
combined_class_rasters <- stack(combined_class_rasters)  

# get sum of values for each pixel in stack (basically combine all rasters
# and generate a value of number of snake spp. per pixel)
combined_class_richness <- sum(combined_class_rasters, na.rm = TRUE)

# write out the category 1 species richness output
combined_outpath <- paste('output/species_richness/', 'combined_categories_', Sys.Date(), sep = '')
writeRaster(combined_class_richness,
            file = combined_outpath,
            format = 'GTiff',
            overwrite = TRUE)

message(paste('Completed species richness surface processing ', Sys.time(), sep = ''))
