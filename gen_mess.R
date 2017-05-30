# clear workspace
rm(list = ls ())

# load packages
pacman::p_load(raster, dismo, rgeos, seegSDM, maptools, ggplot2, colorRamps, reshape2, dplyr, foreach, doMC)

# load functions
source('code/bespoke_functions_cluster.R')

# read in a list of the covariates to be used in the model
cov_list <- read.csv('data/raw/raster/covariates_for_model_cluster.csv',
                     stringsAsFactors = FALSE)

# read in snakelist
snake_list <- read.csv('data/raw/snake_list_cluster.csv',
                       stringsAsFactors = FALSE)

# read in admin 0
admin_0 <- shapefile('data/raw/world_shapefiles/admin2013_0.shp')

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

  # create conservative mess map
  # inform progress
  message(paste('generating conservative MESS', ' (', spp_name, ') ', Sys.time(), sep = ""))
  suppressWarnings(conservative_mess_map <- mess(covs, covariate_stats, full = TRUE))

  # make binary, interpolation/extrapolation surface for the conservative mess
  tmp <- conservative_mess_map[['rmess']] >= 0
  raw_mess <- conservative_mess_map[['rmess']]
  
  # crop to extent
  tmp_masked <- mask(tmp, covs[[1]])
  raw_mess_masked <- mask(raw_mess, covs[[1]])
  
  # write to disk
  outpath <- paste(geotif_mess, spp_name, '_conservative_binary', sep = '')
  writeRaster(tmp_masked, 
              file = outpath, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  outpath_2 <- paste(geotif_mess, spp_name, '_conservative_raw', sep = '')
  writeRaster(raw_mess_masked, 
              file = outpath_2, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  # generate the 1000 MESS surface (uh-oh)
  # run function, specify number of bootstraps, and number of reference points
  message(paste('generating bootstrapped MESS', ' (', spp_name, ') ', Sys.time(), sep = ""))
  bootstrapped_mess <- the_1000_mess_project(n_boot = 10,
                                             in_parallel = TRUE,
                                             n_cores = 50,
                                             covs_extract = covs_extract, 
                                             covs = covs,
                                             occ_dat = records_inside,
                                             eval_plot = TRUE,
                                             plot_outpath = mess_eval)
  
  # write bootstrapped mess to disk
  outpath_3 <- paste(geotif_mess, spp_name, '_bootstrapped', sep = '')
  writeRaster(bootstrapped_mess, 
              file = outpath_3, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  # create a plot of both of the MESS outputs, specify whether to add the
  # reference points to the plot
  plots <- plot_mess(png_mess, 
                     spp_name, 
                     add_points = TRUE,
                     tmp_masked,
                     bootstrapped_mess)
  
  }
  
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
    records_outside$mess <- cellFromXY(tmp_masked, lat_lon)
    
    records_oor_pos <- as.data.frame(records_outside[!(is.na(records_outside$mess)), ])
    records_oor_neg <- as.data.frame(records_outside[is.na(records_outside$mess), ])
    
    neg_outpath <- paste('output/mess_oor_data/', 
                         spp_name, '_OOR_neg_Mess.csv', sep = '')
    
    pos_outpath <- paste('output/mess_oor_data/', 
                         spp_name, '_OOR_pos_Mess.csv', sep = '')
    
    if(nrow(records_oor_pos) != 0) {
    
      write.csv(records_oor_pos,
              pos_outpath,
              row.names = FALSE)
    }
    
    if(nrow(records_oor_neg) != 0) {
      
      write.csv(records_oor_neg,
                neg_outpath,
                row.names = FALSE)  
      
    }
  
  }
  
}  
