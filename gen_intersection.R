# generate new intersection figure
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreign, reshape2, ggplot2, scales, RColorBrewer)

# set working directory
setwd('H:/Work/Snakebite')

# Intersection of:
# -Any snake present
# -No listed antivenom
# -Lowest HAQ decile
# -Travel time of 3 hours

# load species richness surface
species_richness <- raster('output/species_richness/modified_eor_combined_categories_2017-07-25.tif')
# load antivenom naive surface
antivenom <- raster('output/antivenom_coverage/No_specific_antivenom_combined_stack.tif')
# load HAQ decile shapefile
haqi <- shapefile('HAQI/New HAQ/Final_HAQ.shp')
# load population density surface
pop_dens <- raster('rasters/population/Worldpop_GPWv4_Hybrid_201601_Global_Pop_5km_Adj_MGMatched_2015_Hybrid.tif')
# load in travel time surface
travel_surface <- raster('rasters/accessibility/accessibility_50k+_2017-07-31_aggregate_5k_2017_08_09.tif')
# load admin 0 raster
admin_0 <- raster('rasters/admin_0_updated_2017-08-01.tif')
# load admin 2 shapefile
admin_2 <- shapefile("H:/Shapefiles/FAO GAUL 2013/admin2013_2.shp")
# rasterize
admin_2 <- rasterize(admin_2, species_richness, "GAUL_CODE")

###### generate lowest HAQI raster #####
# subset shapefile to create two separate shapefiles;
# one which contains decile 1 only, and another which contains deciles one to three
haqi_decile1 <- haqi[haqi$decile <= 1, ]
haqi_decile1_3 <- haqi[haqi$decile <= 3, ]
haqi_decile1_4 <- haqi[haqi$decile <= 4, ]
haqi_decile1_5 <- haqi[haqi$decile <= 5, ]
haqi_decile1_6 <- haqi[haqi$decile <= 6, ]
haqi_decile1_7 <- haqi[haqi$decile <= 7, ]
haqi_decile1_8 <- haqi[haqi$decile <= 8, ]
haqi_decile1_9 <- haqi[haqi$decile <= 9, ]

# rasterize
haqi_d1_raster <- rasterize(haqi_decile1, species_richness, 1)
haqi_d3_raster <- rasterize(haqi_decile1_3, species_richness, 1)
haqi_d4_raster <- rasterize(haqi_decile1_4, species_richness, 1)
haqi_d5_raster <- rasterize(haqi_decile1_5, species_richness, 1)
haqi_d6_raster <- rasterize(haqi_decile1_6, species_richness, 1)
haqi_d7_raster <- rasterize(haqi_decile1_7, species_richness, 1)
haqi_d8_raster <- rasterize(haqi_decile1_8, species_richness, 1)
haqi_d9_raster <- rasterize(haqi_decile1_9, species_richness, 1)
haqi_d10_raster <- rasterize(haqi, species_richness, 1)

###### make all rasters binary, so that I can stack #####
species_richness[species_richness  >= 1] <- 1
antivenom[antivenom >= 1] <- 1

# convert travel surface
# bin accessibility values to generate mortality likelihoods
# 1. set up matrix
vals <- matrix(ncol = 3,
               c(seq(0, 6000, 60),
                 seq(60, 6060, 60),
                 seq(1, 101, 1)), 
               byrow = FALSE)

# 2. reclassify into mortality bins (suppose this is similar to /60 and rounding vals, but I want
# 61-89 minutes to contribute towards 2% mortality likelihood, opposed to 1%).
travel_surface_mod <- round(travel_surface/60, 0)
travel_surface_mod[travel_surface_mod <= 3] <- 0
travel_surface_mod[travel_surface_mod > 0] <- 1

process_list <- c('haqi_d1_raster', 'haqi_d3_raster', 'haqi_d4_raster', 'haqi_d5_raster',
                  'haqi_d6_raster', 'haqi_d7_raster', 'haqi_d8_raster', 'haqi_d9_raster',
                  'haqi_d10_raster')

d_names <- c('decile_1', 'deciles_1-3', 'deciles_1-4', 'deciles_1-5', 'deciles_1-6', 'deciles_1-7',
             'deciles_1-8', 'deciles_1-9', 'deciles_1-10')

# stack the rasters; create two stacks (one with the antivenom naive only; one without)
for(i in 1:length(process_list)){
  
  # get loop specific raster
  loop_raster <- get(process_list[i])
  decile_name <- d_names[i]
  raster_list_no_antivenom <- c(species_richness, antivenom, loop_raster, travel_surface_mod)
  raster_list <- c(species_richness, loop_raster, travel_surface_mod)
  
  # stack raster lists
  intersection_stack_no_antivenom <- stack(raster_list_no_antivenom)
  intersection_stack <- stack(raster_list)

  # get sum of values for each pixel in stack (basically combine all rasters
  # and generate a value of number of true conditions per cell, 
  # i.e. snake presence = true, antivenom = true)
  intersection_surface <- sum(intersection_stack, na.rm = TRUE)
  intersection_surface_no_antivenom <- sum(intersection_stack_no_antivenom, na.rm = TRUE)

  # convert these into binary surfaces
  intersection_binary <- intersection_surface
  intersection_binary[intersection_binary <=2 ] <- 0
  intersection_binary[intersection_binary == 3] <- 1

  intersection_binary_surface_no_antivenom <- intersection_surface_no_antivenom
  intersection_binary_surface_no_antivenom[intersection_binary_surface_no_antivenom <= 3] <- 0
  intersection_binary_surface_no_antivenom[intersection_binary_surface_no_antivenom == 4] <- 1

  # now generate most vulnerable populations based on these intersections
  # multiply the population by the binary presence/absence surface
  presence_par <- overlay(pop_dens, intersection_binary, 
                          fun = function(pop_dens, intersection_binary){
                          (pop_dens*intersection_binary)})

  presence_no_antivenom <- overlay(pop_dens, intersection_binary_surface_no_antivenom, 
                                   fun = function(pop_dens, intersection_binary_surface_no_antivenom){
                                   (pop_dens*intersection_binary_surface_no_antivenom)})

  
  # crop by this raster
  presence_par_cropped <- crop(presence_par, admin_0)
  presence_par_no_anti_cropped <- crop(presence_no_antivenom, admin_0)
  extent(admin_0) <- extent(presence_par_cropped)

  # get national PAR values
  global_par_intersection <- zonal(presence_par_cropped, admin_0, fun = 'sum', na.rm = TRUE)
  global_par_inter_no_anti <- zonal(presence_par_no_anti_cropped, admin_0, fun = 'sum', na.rm = TRUE)
  a2_par_intersection <- zonal(presence_par, admin_2, fun = 'sum', na.rm = TRUE)
  a2_par_inter_no_anti <- zonal(presence_no_antivenom, admin_2, fun = 'sum', na.rm = TRUE)

  # write out rasters
  writeRaster(presence_par, 
              file = paste0('output/intersection/15.04.2018/pixel_level_intersection_', decile_name, '_excluding_antivenom_', Sys.Date()),
              format = 'GTiff',
              overwrite = TRUE)
  writeRaster(presence_no_antivenom, 
              file = paste0('output/intersection/15.04.2018/pixel_level_intersection_', decile_name, '_including_antivenom_', Sys.Date()),
              format = 'GTiff',
              overwrite = TRUE)
  writeRaster(intersection_surface, 
              file = paste0('output/intersection/15.04.2018/intersection_surface_', decile_name, 'excluding_antivenom_', Sys.Date()),
              format = 'GTiff',
              overwrite = TRUE)
  writeRaster(intersection_surface_no_antivenom, 
              file = paste0('output/intersection/15.04.2018/intersection_surface_', decile_name, 'including_antivenom_', Sys.Date()),
              format = 'GTiff',
              overwrite = TRUE)
  write.csv(global_par_intersection,
            file = paste0('output/intersection/15.04.2018/pop_at_risk_', decile_name, '_excluding_antivenom_', Sys.Date(), '.csv'),
            row.names = FALSE)
  write.csv(global_par_inter_no_anti,
            file = paste0('output/intersection/15.04.2018/pop_at_risk_', decile_name, '_including_antivenom_', Sys.Date(), '.csv'),
            row.names = FALSE)
  write.csv(a2_par_intersection,
            file = paste0('output/intersection/15.04.2018/pop_at_risk_', decile_name, '_admin_2_excluding_antivenom_', Sys.Date(), '.csv'),
            row.names = FALSE)
  write.csv(a2_par_inter_no_anti,
            file = paste0('output/intersection/15.04.2018/pop_at_risk_', decile_name, '_admin_2_including_antivenom_', Sys.Date(), '.csv'),
            row.names = FALSE)
}
