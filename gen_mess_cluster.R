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

# list species for which we have occurrence data
spp_list <- list.files('data/raw/species_data/',
                       pattern = '*_raw.csv$',
                       full.names = FALSE)

# specify outpath for distribution plots
var_path <- 'output/distribution_plots/'

# specify outpath for MESS geotiffs
geotif_mess <- 'output/mess_geotiff/'

# specify outpath for MESS pngs
png_mess <- 'output/mess_png/'

# specify outpath for MESS SI plots
mess_si <- 'output/mess_si/'

# specify outpath for MESS evaluation
mess_eval <- 'output/mess_evaluation/'

# specify outpath for modified shapefiles
new_shp_path <- 'output/modified_ranges/'

# loop through and create a MESS for each species, including plots showing 
# variance within environmental covariates across referenced occurrence data
for(i in 1:length(spp_list)){
  
  # get species info from spp_list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('_raw.csv', '', spp_name)
  
  # inform progress
  message(paste('Processing species ', i, ' of ', length(spp_list), ' (', spp_name, ') ', Sys.time(), sep = ""))
  
  sub <- snake_list[snake_list$split_spp == spp_name, ]

  # create extent shapefile
  ext <- generate_ext(sub, admin_0)
  
  # read in EOR shapefile for each species, and dissolve if >1 polygon
  range <- prepare_eor(sub)
  raw_range <- shapefile(sub$shapefile_path)

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
  
  # remove any NA extracts (if any covariate value is a NA at a particular row, remove it)
  covs_extract <- covs_extract[complete.cases(covs_extract), ]

  # generate the 1000 MESS surface 
  # run function, specify number of bootstraps, and number of reference points
  # specify if PCC and covariate variance should also be quantified using 'eval_plot' and 
  # 'variance' arguments
  message(paste('Generating bootstrapped MESS', ' (', spp_name, ') ', Sys.time(), sep = ""))
  
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
                                             # where to save the PCC plots
                                             pcc_outpath = mess_eval,
                                             # run an evaluation plot to assess variance?
                                             variance = TRUE,
                                             # where to save the variance plots
                                             var_outpath = var_path)                                 

  # convert bootstrapped MESS into a binary surface; threshold set = 95%
  binary_95 <- bootstrapped_mess
  binary_95[binary_95 < 95] <- 0
  
  # convert bootstrapped MESS into a binary surface; threshold set = 90%
  binary_90 <- bootstrapped_mess
  binary_90[binary_90 < 90] <- 0
  
  # convert bootstrapped MESS into a binary surface; threshold set = 75%
  binary_75 <- bootstrapped_mess
  binary_75[binary_75 < 75] <- 0
  
  # stack the different thresholded MESS
  mess_stack <- stack(bootstrapped_mess, binary_95, binary_90, binary_75)
  
  # write out stacked binary surface
  # bands in the following geotiff are now
  # band 1: boostrapped MESS (vals ranging from 0-100)
  # band 2: binary MESS 95% threshold
  # band 3: binary MESS 90% threshold
  # band 4: binary MESS 75% threshold
  stacked_outpath <- paste(geotif_mess, spp_name, '_stacked_bootstrapped_threshold', sep = '')
  writeRaster(mess_stack, 
              file = stacked_outpath,
              format = 'GTiff',
              overwrite = TRUE)

  }
}

# initialize the cluster (50 cores)
registerDoMC(50)

message(paste('Processing species ranges on 50 cores ', Sys.time(), sep = ""))

# run the following buffer section on all species in parallel
stage_2 <- foreach(i = 1:length(spp_list)) %dopar% {
  
  # get species info from spp_list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('_raw.csv', '', spp_name)
  
  # generate subset
  sub <- snake_list[snake_list$split_spp == spp_name, ]
  
  # create extent shapefile
  ext <- generate_ext(sub, admin_0)
  
  # read in EOR shapefile for each species, and dissolve if >1 polygon
  range <- prepare_eor(sub)
  raw_range <- shapefile(sub$shapefile_path)

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
  
  # create a buffered polygon around points within the binary bootstrapped MESS
  # start by classifying as being OOR MESS +ve or OOO MESS -ve; this requires:
  # 1. >= 5 records within the range (otherwise a MESS wasn't generated above)
  # 2. >0 records outside of the range (otherwise there's no range modification)
  if(nrow(records_inside) >= 5){
    
    # turn into a spatial points dataframe
    coordinates(records_outside) <- c("longitude", "latitude")
    
    # project coordinates in same projection as shape
    proj4string(records_outside) <- proj4string(range)
    
    # get a matrix of coordinates for use in extract function
    lat <- records_outside$latitude
    lon <- records_outside$longitude
    
    lat_lon <- cbind(lon,
                     lat)
    
    # get an index for occurrence records within species' range
    # load in the 95% surface
    mess_path <- paste(geotif_mess, spp_name, '_stacked_bootstrapped_threshold.tif', sep = '')
    
    bootstrapped_mess <- raster(mess_path, band = 1)
    binary_95 <- raster(mess_path, band = 2)
    
    vals_95threshold <- extract(binary_95, lat_lon)
    vals_95threshold <- replace(vals_95threshold, vals_95threshold == 0, NA)
    
    # attach 'vals' to the records_outside dataframe
    records_outside <- as.data.frame(records_outside)
    records_outside$mess_95 <- vals_95threshold
    
    # subset to create two dataframes
    # 1. records outside of the EOR, but within MESS interpolation space
    records_oor_pos <- records_outside[!(is.na(records_outside$mess_95)), ]
    # 2. records outside of the EOR, and in MESS extrapolation space
    records_oor_neg <- records_outside[is.na(records_outside$mess_95), ]
    
    # generate the same stats for the 90% and 75% threshold binary MESS
    binary_90 <- raster(mess_path, band = 3)
    vals_90threshold <- extract(binary_90, lat_lon)
    vals_90threshold <- replace(vals_90threshold, vals_90threshold == 0, NA)
    
    # attach 'vals' to the records_outside dataframe
    records_outside$mess_90 <- vals_90threshold
    
    # subset to create two dataframes
    # 1. records outside of the EOR, but within MESS interpolation space
    records_oor_pos_90 <- records_outside[!(is.na(records_outside$mess_90)), ]
    # 2. records outside of the EOR, and in MESS extrapolation space
    records_oor_neg_90 <- records_outside[is.na(records_outside$mess_90), ]
    
    # generate the same stats for the 90% and 75% threshold binary MESS
    binary_75 <- raster(mess_path, band = 4)
    vals_75threshold <- extract(binary_75, lat_lon)
    vals_75threshold <- replace(vals_75threshold, vals_75threshold == 0, NA)
    
    # attach 'vals' to the records_outside dataframe
    records_outside$mess_75 <- vals_75threshold
    
    # subset to create two dataframes
    # 1. records outside of the EOR, but within MESS interpolation space
    records_oor_pos_75 <- records_outside[!(is.na(records_outside$mess_75)), ]
    # 2. records outside of the EOR, and in MESS extrapolation space
    records_oor_neg_75 <- records_outside[is.na(records_outside$mess_75), ]
    
    # if there's any records which are within the MESS interpolation space
    # buffer these records, and merge them on to the current EOR shapefiles
    if(nrow(records_oor_pos) != 0) {
      
      # get matrix of coordinates
      buf_lat <- records_oor_pos$latitude
      buf_lon <- records_oor_pos$longitude
      
      buffer_pts <- cbind(buf_lon,
                          buf_lat)
      
      # run buffer function, specify a radius to buffer by (in km), and 
      # which mess layer/threshold to clip by
      modified_poly <- bufferMESSpositives(range = raw_range,
                                           coords = buffer_pts,
                                           radius = 100,
                                           mess = binary_95,
                                           admin_0,
                                           outpath = new_shp_path,
                                           spp_name = spp_name)
    
      # get final extent for a species
      new_shp_path <- paste(new_shp_path, spp_name, sep = "")
      temp_shp <- shapefile(new_shp_path)
      
      # get countries in new extent
      final_countries <- unique(temp_shp@data$COUNTRY_ID)
      
      # remove any disputed regions (iso XXX)
      final_countries <- final_countries[!(final_countries == 'XXX')]
      final_countries <- final_countries[!is.na(final_countries)]
      
      # create extent shapefile
      final_ext <- admin_0[admin_0$COUNTRY_ID %in% final_countries, ]
      
    }
    
    # generate a dataframe with stats for each species
    total_obs <- nrow(locations)
    within_eor <- nrow(records_inside)
    outside_eor <- nrow(records_outside)
    within_95 <- nrow(records_oor_pos)
    within_90 <- nrow(records_oor_pos_90)
    within_75 <- nrow(records_oor_pos_75)
    
    spec_stats <- data.frame(species = spp_name, 
                             total_occurrence_records = total_obs,
                             total_occ_within_eor = within_eor,
                             total_occ_outside_eor = outside_eor,
                             mess_positive_95 = within_95,
                             mess_positive_90 = within_90,
                             mess_positive_75 = within_75)
    
    mess_stats_path <- paste(png_mess, spp_name, '_mess_stats_', Sys.Date(), '.csv', sep = "")
    
    write.csv(spec_stats,
              mess_stats_path,
              row.names = FALSE)
    
  } else {
    
    spec_stats <- NULL
    
  }
  
  # generate SI plots - this requires a MESS to have been constructed, which depends on a spp
  # having 5 or more occurrence records, as per above.
  
  if(nrow(records_inside) >= 5){
  
  # create a plotting window
  png_name <- paste(png_mess, spp_name, '_species_mess_maps_', Sys.Date(), '.png', sep = "")
  
  png(png_name,
      width = 450,
      height = 400,
      units = 'mm',
      res = 300)
  par(mfrow = c(2, 2),
      oma = c(2, 2, 2, 2))
  
  # define species name
  title <- gsub('_', ' ', spp_name)
  
  ### plot the bootstrapped MESS 
  plot(bootstrapped_mess,
       main = "A)", adj = 0,
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
       lwd = 0.2)
  
  a_text <- paste('100 Bootstrapped MESS\n(constructed using', nrow(records_inside), 'occurrence records)', sep = " ")
  
  title(xlab = a_text, line = 2, cex.lab = 1.25)
  
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n')
  
  ### plot the 95% binary MESS
  plot(binary_95,
       main = "B)", adj = 0,
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
       lwd = 0.2)
  
  title(xlab = 'Binary bootstrapped MESS (95% threshold)', line = 1, cex.lab = 1.25)
  
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n')
  
  ### plot the 95% binary MESS with points
  plot(binary_95,
       main = "C)", adj = 0,
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
       lwd = 0.2)
  
  c_text <- paste('Binary bootstrapped MESS (95% threshold) \nwith', nrow(records_outside), 
                  'outside of range occurrence records', sep = " ")
  
  title(xlab = c_text, line = 2, cex.lab = 1.25)
  
  # add points on top
  points(records_oor_neg$longitude, records_oor_neg$latitude, pch = 4, cex = 0.5, col = 'blue', lwd = 0.5)
  points(records_oor_pos$longitude, records_oor_pos$latitude, pch = 20, cex = 0.5, col = '#D93529')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Outside range, MESS +ve", "Outside range, MESS -ve"),
         pch = c(15, 15, 20, 4),
         col = c("springgreen4","gainsboro", "#D93529", "blue"), bty = 'n')
  
  ### plot the new range shapefile, ontop of binary bootstrapped MESS
  if(nrow(records_oor_pos) != 0) {
    plot(final_ext,
         main = "D)", adj = 0,
         col = '#f2f2f2',
         border = NA)
    
    plot(modified_poly,
         add = TRUE,
         col = '#00a600',
         border = NA)
    
    plot(range,
         add = TRUE,
         col = '#c1de29',
         border = NA)
    
    plot(final_ext,
         add = TRUE,
         border = 'gray45',
         lty = 1,
         lwd = 0.2)
  
  } else {
    
    plot(ext,
         main = "D)", adj = 0,
         col = '#f2f2f2',
         border = NA)
    
    plot(range,
         add = TRUE,
         col = '#c1de29',
         border = NA)
    
    plot(ext,
         add = TRUE,
         border = 'gray45',
         lty = 1,
         lwd = 0.2)
    
  }
  
  mess_positive_title <- paste('Suggested ammended range \n(incorporating', nrow(records_oor_pos), 
                               'outside of range MESS +ve records)', sep = " ")
  
  legend('bottomleft', c("Current EOR","Proposed addition"), 
         pch = c(15, 15),
         col = c("#c1de29","#00a600"), bty = 'n')
  
  title(xlab = mess_positive_title, line = 2, cex.lab = 1.25)
  
  mtext(bquote(~italic(.(title))), side = 3, line = -1, outer = TRUE, cex = 2, font = 2)
  
  dev.off()      

  ## SI plot of the different threshold binary-bootstrapped MESS
  # create plotting window
  png_name <- paste(mess_si, spp_name, '_SI_threshold_mess_maps_', Sys.Date(), '.png', sep = "")
  
  png(png_name,
      width = 450,
      height = 200,
      units = 'mm',
      res = 300)
  par(mfrow = c(1, 2), 
      oma = c(2, 2, 2, 2))
  
  # generate species title
  title <- gsub('_', ' ', spp_name)
  
  ### plot the raw bootstrapped MESS 
  plot(bootstrapped_mess,
       main = "A)", adj = 0,
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
  
  title(xlab = '100 Bootstrapped MESS', line = 1, cex.lab = 1)
  
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n', pt.cex = 1.5)

  ### plot the 95% binary MESS with points
  binary_95[binary_95 >= 95] <- 1
  
  breakpoints <- c(0, 0.5, 1)
  colours <- c('#f2f2f2', 'gray70', 'gray70')
  
  plot(binary_95,
       breaks = breakpoints,
       col = colours,
       main = "B)", adj = 0,
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
       lwd = 0.2)
  
  c_text <- paste('MESS threshold classification of \n', nrow(records_outside), 
                  'outside of range occurrence records', sep = " ")
  
  title(xlab = c_text, line = 2, cex.lab = 1)
  
  # add points on top
  points(records_oor_neg_75$longitude, records_oor_neg_75$latitude, pch = 4, cex = 0.5, col = '#ff7f27', lwd = 0.5)
  points(records_oor_pos_75$longitude, records_oor_pos_75$latitude, pch = 21, cex = 0.5, col = 'gray58', bg = '#22b14c', lwd = 0.25)
  points(records_oor_pos_90$longitude, records_oor_pos_90$latitude, pch = 21, cex = 0.5, col = 'gray58', bg = '#E71D36', lwd = 0.25)
  points(records_oor_pos$longitude, records_oor_pos$latitude, pch = 21, cex = 0.5, col = 'gray58', bg = '#00a2e8', lwd = 0.25)
  
  legend("bottomleft", c("Within the 95%, 90% and 75% MESS thresholds", "Within the 90% and 75% MESS thresholds only", 
                         "Within the 75% MESS threshold only", "Outside of all MESS thresholds"), 
         pch = c(20, 20, 20, 4),
         col = c("#00a2e8","#E71D36", "#22b14c", "#ff7f27"), bty = 'n', cex = 0.8, pt.cex = 1)
  
  mtext(bquote(~italic(.(title))), side = 3, line = -1, outer = TRUE, cex = 2, font = 2)
  
  dev.off()
  
  }
  
  return(spec_stats)
  
  }

message(paste('Completed processing species ranges on 50 cores ', Sys.time(), sep = ""))

# merge species stats
spp_list <- list.files('data/raw/species_data/',
                       pattern = '*_raw.csv$',
                       full.names = FALSE)

temp_frame <- NULL

# loop through, read and bind
for(i in 1:length(spp_list)){

    temp <- read.csv(spp_list[i],
                     stringsAsFactors = FALSE)  
    
    temp_frame <- rbind(temp_frame,
                        temp)
    
}

# write out merged stats
mess_stats_path <- paste(png_mess, 'combined_mess_stats_', Sys.Date(), '.csv', sep = "")
write.csv(temp_frame,
          mess_stats_path,
          row.names = FALSE) 

