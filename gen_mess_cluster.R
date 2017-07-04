# clear workspace
rm(list = ls ())

# set the seed to ensure reproducibility
set.seed(1)

# load packages
pacman::p_load(raster, dismo, rgeos, seegSDM, maptools, ggplot2, colorRamps, reshape2, dplyr, foreach, doMC)

# load functions
source('code/bespoke_functions_cluster.R')

# read in a list of the covariates to be used in the model
# cov_list <- read.csv('data/raw/raster/covariates_for_model_cluster.csv',
#                      stringsAsFactors = FALSE)
cov_list <- read.csv('data/raw/raster/bio_climatic_covariates_for_model_cluster.csv',
                     stringsAsFactors = FALSE)

# read in snakelist
snake_list <- read.csv('data/raw/snake_list_cluster.csv',
                       stringsAsFactors = FALSE)

# read in admin 0
admin_0 <- shapefile('data/raw/world_shapefiles/admin2013_0.shp')

# read in sea shapefile
sea_shp <- shapefile('data/raw/world_shapefiles/sea.shp')

# list species for which we have occurrence data
spp_list <- list.files('data/raw/species_data/',
                       pattern = '*_raw.csv$',
                       full.names = FALSE)

# specify outpath for distribution plots
distribution_path <- 'output/distribution_plots/'

# specify outpath for MESS geotiffs
geotif_mess <- 'output/mess_geotiff/'

# specify outpath for MESS pngs
png_mess <- 'output/mess_png/'

# specify outpath for MESS SI plots
mess_si <- 'output/mess_si/'

# specify outpath for MESS evaluation
mess_eval <- 'output/mess_evaluation/'

# loop through and create a MESS for each species, including plots showing 
# variance within environmental covariates across referenced occurrence data
for(i in 1:length(spp_list)){
  
  # get species info from spp_list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('_raw.csv', '', spp_name)
  
  # inform progress
  message(paste('processing species ', i, ' of ', length(spp_list), ' (', spp_name, ') ', Sys.time(), sep = ""))
  
  sub <- snake_list[snake_list$split_spp == spp_name, ]

  # create extent shapefile
  ext <- generate_ext(sub, admin_0)
  
  # read in EOR shapefile for each species, and dissolve if >1 polygon
  range <- prepare_eor(sub)

  # prepare covariates for analysis
  covs <- prepare_covariates(cov_list, ext, range)
  
  # get presence records for species
  locations <- load_occurrence(spp_name, range)

  if(nrow(locations) !=0){
  
  # get an index for occurrence records within species' range
  records_inside <- !is.na(over(locations, as(range, "SpatialPolygons")))
  
  records_outside <- as.data.frame(locations[records_inside == FALSE, ])
  records_inside <- as.data.frame(locations[records_inside == TRUE, ])
  
  # clean locations to remove 0,0 coordinates  
  records_outside <- records_outside[!records_outside$latitude == 0 & !records_outside$longitude == 0, ]
  
  } else {
    
  records_inside <- locations  
  records_outside <- locations  
  
  }
  
  # generate a MESS, using a prerequisite of a minimum of 5 data points within the species range
  if(nrow(records_inside) >= 5){

  # extract covariate data for all of that species' presence records, within the WHO EOR
  covs_extract <- as.data.frame(extract(covs, records_inside[c("longitude", "latitude")]))
  
  # subsample from this extracted covariate data, to generate 1000 bootstraps
  # then plot the distribution of the max and mix values for each covariate
  # to measure variation in the input data
  covariate_stats <- measure_variance(n_boot = 1000, 
                                      covs_extract,
                                      distribution_path,
                                      spp_name)

  # ### create conservative mess map
  # # inform progress
  # message(paste('generating conservative MESS', ' (', spp_name, ') ', Sys.time(), sep = ""))
  # suppressWarnings(conservative_mess_map <- mess(covs, covariate_stats, full = TRUE))
  # 
  # # make binary, interpolation/extrapolation surface for the conservative mess
  # tmp <- conservative_mess_map[['rmess']] >= 0
  # raw_mess <- conservative_mess_map[['rmess']]
  # 
  # # crop to extent
  # tmp_masked <- mask(tmp, covs[[1]])
  # raw_mess_masked <- mask(raw_mess, covs[[1]])
  # 
  # # write to disk
  # outpath <- paste(geotif_mess, spp_name, '_conservative_binary', sep = '')
  # writeRaster(tmp_masked, 
  #             file = outpath, 
  #             format = 'GTiff', 
  #             overwrite = TRUE)
  # 
  # outpath_2 <- paste(geotif_mess, spp_name, '_conservative_raw', sep = '')
  # writeRaster(raw_mess_masked, 
  #             file = outpath_2, 
  #             format = 'GTiff', 
  #             overwrite = TRUE)
  
  # generate the 1000 MESS surface 
  # run function, specify number of bootstraps, and number of reference points
  message(paste('generating bootstrapped MESS', ' (', spp_name, ') ', Sys.time(), sep = ""))
  
                                             # number of bootstraps
  bootstrapped_mess <- the_1000_mess_project(n_boot = 100,
                                             # execute in parallel?
                                             in_parallel = TRUE,
                                             # number of cores
                                             n_cores = 50,
                                             # data frame of extracted covariate values
                                             covs_extract = covs_extract,
                                             # stack of covariates
                                             covs = covs,
                                             # dataframe with species occurrence records
                                             occ_dat = records_inside,
                                             # run an evaluation plot to assess PCC?
                                             eval_plot = TRUE,
                                             # where to save the plots
                                             plot_outpath = mess_eval)
      
  # write bootstrapped mess to disk
  outpath_3 <- paste(geotif_mess, spp_name, '_bootstrapped', sep = '')
  writeRaster(bootstrapped_mess, 
              file = outpath_3, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  # convert bootstrapped MESS into a binary surface; threshold set = 95%
  binary_95 <- bootstrapped_mess
  binary_95[binary_95 < 95] <- 0
  
  # write out binary
  outpath_4 <- paste(geotif_mess, spp_name, '_binary_bootstrapped_95', sep = '')
  writeRaster(binary_95, 
              file = outpath_4, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  # convert bootstrapped MESS into a binary surface; threshold set = 90%
  binary_90 <- bootstrapped_mess
  binary_90[binary_90 < 90] <- 0
  
  # write out binary
  outpath_5 <- paste(geotif_mess, spp_name, '_binary_bootstrapped_90', sep = '')
  writeRaster(binary_90, 
              file = outpath_5, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  # convert bootstrapped MESS into a binary surface; threshold set = 75%
  binary_75 <- bootstrapped_mess
  binary_75[binary_75 < 75] <- 0
  
  # write out binary
  outpath_6 <- paste(geotif_mess, spp_name, '_binary_bootstrapped_95', sep = '')
  writeRaster(binary_75, 
              file = outpath_6, 
              format = 'GTiff', 
              overwrite = TRUE)

  # create a buffered polygon around points within the binary bootstrapped MESS
  # start by classifying as being OOR MESS +ve or OOO MESS -ve
  if(nrow(records_outside) != 0){
    
    # turn into a spatial points dataframe
    coordinates(records_outside) <- c("longitude", "latitude")
    
    # project coordinates in same projection as shape
    proj4string(records_outside) <- proj4string(range)
    
    lat <- records_outside$latitude
    lon <- records_outside$longitude
    
    lat_lon <- cbind(lon,
                     lat)
    
    # get an index for occurrence records within species' range
    vals <- extract(binary_95, lat_lon)
    
    records_outside <- as.data.frame(records_outside)
    records_outside$mess <- vals
    
    records_oor_pos <- records_outside[!(is.na(records_outside$mess)), ]
    records_oor_neg <- records_outside[is.na(records_outside$mess), ]
    
    if(nrow(records_oor_pos) != 0) {
      
      buf_lat <- records_oor_pos$latitude
      buf_lon <- records_oor_pos$longitude
      
      buffer_pts <- cbind(buf_lon,
                          buf_lat)
      
      modified_poly <- bufferMESSpositives(range,
                                           coords = buffer_pts,
                                           radius = 20,
                                           sea = sea_shp)
      
    }
    
  }
  
  ### SI plots
  # create a plot of both of the MESS outputs, specify whether to add the reference points to the plot
  # create a plotting window to plot both of the surfaces
  png_name <- paste(png_mess, spp_name, '_species_mess_maps_', Sys.Date(), '.png', sep = "")
  
  png(png_name,
      width = 450,
      height = 400,
      units = 'mm',
      res = 300)
  par(mfrow = c(2, 2))
  
  # plot the bootstrapped MESS 
  plot(bootstrapped_mess,
       main = bquote(~italic(.(title))),
       legend = TRUE,
       axes = FALSE,
       box = FALSE)
  plot(range,
       add = TRUE,
       border = 'black',
       lty = 1,
       lwd = 1.5)
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.5)
  
  title(xlab = 'Bootstrapped MESS (100 bootstraps)', line = 0)
  
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n')
  
  # plot the 95% binary MESS
  plot(binary_95,
       main = bquote(~italic(.(title))),
       legend = TRUE,
       axes = FALSE,
       box = FALSE)
  plot(range,
       add = TRUE,
       border = 'black',
       lty = 1,
       lwd = 1)
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.5)
  
  title(xlab = 'Binary bootstrapped MESS (95% threshold)', line = 0)
  
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n')
  
  # plot the 95% binary MESS with points
  plot(binary_95,
       main = bquote(~italic(.(title))),
       legend = TRUE,
       axes = FALSE,
       box = FALSE)
  plot(range,
       add = TRUE,
       border = 'black',
       lty = 1,
       lwd = 1.5)
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.5)
  
  title(xlab = 'Binary bootstrapped MESS (95% threshold)', line = 0)
  
  # add points on top
  points(records_outside$longitude, records_outside$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_inside$longitude, records_inside$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within range", "Outside range"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  
  
  dev.off()      

  }
  
  
  
}  
