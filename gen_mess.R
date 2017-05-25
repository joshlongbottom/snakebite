# clear workspace
rm(list = ls ())

# load packages
pacman::p_load(raster, dismo, rgeos, seegSDM, maptools, ggplot2, colorRamps, reshape2, dplyr)

# set working directory
setwd('Z:/users/joshua/Snakebite/')

# load functions
source('snakebite/bespoke_functions.R')

# read in a list of the covariates to be used in the model
cov_list <- read.csv('rasters/covariates_for_model.csv',
                     stringsAsFactors = FALSE)

# read in snakelist
snake_list <- read.csv('snakebite/snake_list.csv',
                       stringsAsFactors = FALSE)

# read in admin 0
admin_0 <- shapefile('World shapefiles/admin2013_0.shp')

# list species for which we have occurrence data
spp_list <- list.files('output/species_occurrence_plots/species data/',
                       pattern = '*_raw.csv$',
                       full.names = FALSE)

# specify outpath for distribution plots
distribution_path <- 'output/species_mess_maps/distribution plots/'

# specify outpath for MESS geotiffs
geotif_mess <- 'output/species_mess_maps/files/'


# loop through and create stacks for each extent
# for(i in 1:length(spp_list)){
for(i in 1:1){
  
  # get species info from spp_list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('_raw.csv', '', spp_name)
  
  # inform progress
  message(paste('processing file ', i, ' of length ', length(spp_list), ' (', spp_name, ') ', sep = ""))
  
  sub <- snake_list[snake_list$split_spp == spp_name, ]

  # create extent shapefile
  ext <- generate_ext(sub, admin_0)
  
  # read in EOR shapefile for each species, and dissolve if >1 polygon
  range <- prepare_eor(sub)

  # prepare covariates for analysis
  covs <- prepare_covariates(cov_list, ext, range)
  
  # get presence records for species
  locations <- load_occurrence(spp_name, range)
  
  # get an index for occurrence records within species' range
  records_inside <- !is.na(over(locations, as(range, "SpatialPolygons")))
  
  records_outside <- as.data.frame(locations[records_inside == FALSE, ])
  records_inside <- as.data.frame(locations[records_inside == TRUE, ])
  
  # generate a MESS, using a prerequisite of a minimum of 5 data points within the species range
  if(nrow(records_inside) >= 5){

  # extract covariate data for all of that species' presence records, within the WHO EOR
  covs_extract <- as.data.frame(extract(covs, records_inside[c("longitude", "latitude")]))
  
  # subsample from this extracted covariate data, to generate 1000 bootstraps
  # then plot the distribution of the max and mix values for each covariate
  # to measure variation in the input data
  covariate_stats <- measure_variance(1000, 
                                      covs_extract,
                                      distribution_path,
                                      spp_name)

  # create conservative mess map, and time the process
  start_time <- Sys.time()
  
  suppressWarnings(conservative_mess_map <- mess(covs, covariate_stats, full = TRUE))
  
  end_time <- Sys.time()
  
  conservative_mess_time <- end_time - start_time
  
  message(paste('conservative MESS for', spp_name, 'generated in:', conservative_mess_time, 'minute(s)'))
  
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
  
  outpath_2 <- paste(geotif_mess, spp_name, 'conservative_raw', sep = '')
  writeRaster(raw_mess_masked, 
              file = outpath_2, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  # generate the 1000 MESS surface (uh-oh), and time
  start_time <- Sys.time()
  
  # run function, specify number of bootstraps
  bootstrapped_mess <- the_1000_mess_project(10, covs_extract, covs)
  
  end_time <- Sys.time()
  
  bootstrapped_mess_time <- end_time - start_time
  
  message(paste('bootstrapped MESS for', spp_name, '(10 straps) generated in:', conservative_mess_time, 'minute(s)'))
  
  # create a plotting window to plot both of the surfaces
  png_name <- paste('output/species_mess_maps/', spp_name, '_species_mess_map_', Sys.Date(), '.png', sep = "")

  png(png_name,
      width = 450,
      height = 250,
      units = 'mm',
      res = 300)
  par(mfrow = c(1, 2))
  
  # plot the conservative MESS
  title <- gsub('_', ' ', spp_name)
  
  plot(tmp_masked, 
       main = bquote(~italic(.(title))),
       legend = FALSE,
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
  
  title(xlab = 'Conservative MESS', line = 0)
  
  # add points on top
  points(records_outside$longitude, records_outside$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_inside$longitude, records_inside$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within range", "Outside range"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  # plot the 1000-MESS MESS 
  plot(bootstrapped_mess,
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
   
   
  title(xlab = 'Bootstrapped MESS', line = 0)
  
  dev.off()
  
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
    
    neg_outpath <- paste('output/species_occurrence_plots/species data oor/', 
                         spp_name, 'OOR_neg_Mess.csv', sep = '')
    
    pos_outpath <- paste('output/species_occurrence_plots/species data oor/', 
                         spp_name, 'OOR_pos_Mess.csv', sep = '')
    
    if(nrow(records_oor_pos) != 0) {
    
      write.csv(records_oor_pos,
              pos_outpath,
              row.names = FALSE)
    }
    
    if(nrow(records_oor_neg) != 0) {
      
      write.csv(records_oor_pos,
                neg_outpath,
                row.names = FALSE)  
  }
  
  }
  
  # rm(records_outside,
  #    records_inside, 
  #    records_oor_pos,
  #    records_oor_neg)
  
}  

# remind where the distribution plots have been saved (forcing to re-direct if necessary)
message(paste('distribution plots saved in the specified directory:', distribution_path, sep = " "))
