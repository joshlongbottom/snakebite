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

# run the following buffer section on all species in parallel
stage_2 <- foreach(i = 1:length(species_list)) %dopar% {
  
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
  
  # create a buffered polygon around points within the binary bootstrapped MESS
  # start by classifying as being OOR MESS +ve or OOO MESS -ve
  if(nrow(records_outside) != 0){
    
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
      modified_poly <- bufferMESSpositives(raw_range,
                                           coords = buffer_pts,
                                           radius = 100,
                                           mess = binary_95,
                                           admin_0,
                                           outpath = new_shp_path)
    
    }
    
  }
  
  ### SI plots
  # changed plots (07.07.2017) to remove species name from plot main title, and to
  # instead label each plot in the panel (A:D), and have one central title for all plots
  # create a plot of both of the MESS outputs, specify whether to add the reference points to the plot
  # create a plotting window to plot both of the surfaces
  message(paste('Plotting', spp_name, 'outputs ', Sys.time(), sep = " "))
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
       lwd = 1.5)
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.5)
  
  title(xlab = 'Bootstrapped MESS (100 bootstraps)', line = 0, cex.lab = 1.25)
  
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
       lwd = 0.5)
  
  title(xlab = 'Binary bootstrapped MESS (95% threshold)', line = 0, cex.lab = 1.25)
  
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
       lwd = 1.5)
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.5)
  
  title(xlab = 'Binary bootstrapped MESS (95% threshold) with occurrence records', line = 0, cex.lab = 1.25)
  
  # add points on top
  points(records_oor_neg$longitude, records_oor_neg$latitude, pch = 20, cex = 0.5, col = 'dimgray')
  points(records_inside$longitude, records_inside$latitude, pch = 20, cex = 0.5, col = 'blue')
  points(records_oor_pos$longitude, records_oor_pos$latitude, pch = 20, cex = 0.75, col = '#D93529')

  legend('bottomleft', c("Interpolation","Extrapolation", "Within range", "Outside range, MESS +ve", "Outside range, MESS -ve"), 
         pch = c(15, 15, 20, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529", "dimgray"), bty = 'n')
  
  ### plot the new range shapefile, ontop of binary bootstrapped MESS
  plot(binary_95,
       main = "D)", adj = 0,
       legend = FALSE,
       axes = FALSE,
       box = FALSE)
  
  if(nrow(records_oor_pos) != 0) {
    
    plot(modified_poly,
         add = TRUE,
         border = 'black',
         lty = 1,
         lwd = 1.5)
    
  } else {
    
    plot(range,
         add = TRUE,
         border = 'black',
         lty = 1,
         lwd = 1.5)
  }
  
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.5)
  
  title(xlab = 'Suggested ammended range (incorporating outside of range MESS +ve records)', line = 0, cex.lab = 1.25)
  
  mtext(bquote(~italic(.(title))), side = 3, line = -1, outer = TRUE, cex = 2, font = 2)
  
  dev.off()      

  ## SI plot of the different threshold binary-bootstrapped MESS
  # create plotting window
  png_name <- paste(mess_si, spp_name, '_SI_threshold_mess_maps_', Sys.Date(), '.png', sep = "")
  
  png(png_name,
      width = 450,
      height = 400,
      units = 'mm',
      res = 300)
  par(mfrow = c(2, 2), 
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
  
  title(xlab = 'Bootstrapped MESS (100 bootstraps)', line = 0, cex.lab = 1.25)
  
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n')

  ### plot the 95% binary MESS with points
  plot(binary_95,
       main = "B)", adj = 0,
       legend = FALSE,
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
  
  title(xlab = 'Binary bootstrapped MESS (95% threshold)', line = 0, cex.lab = 1.25)
  
  # add points on top
  # these are the points which are considered OOEOR MESS+ve/-ve for the 95% threshold
  points(records_oor_neg$longitude, records_oor_neg$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_oor_pos$longitude, records_oor_pos$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within the 95% MESS threshold", "Outside of the 95% MESS threshold"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  ### plot the 90% binary MESS with points
  plot(binary_90,
       main = "C)", adj = 0,
       legend = FALSE,
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
  
  title(xlab = 'Binary bootstrapped MESS (90% threshold)', line = 0, cex.lab = 1.25)
  
  # add points on top
  # these are the points which are considered OOEOR MESS+ve/-ve for the 95% threshold
  points(records_oor_neg_90$longitude, records_oor_neg_90$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_oor_pos_90$longitude, records_oor_pos_90$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within the 90% MESS threshold", "Outside of the 90% MESS threshold"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  ### plot the 75% binary MESS with points
  plot(binary_75,
       main = "D)", adj = 0,
       legend = FALSE,
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
  
  title(xlab = 'Binary bootstrapped MESS (75% threshold)', line = 0, cex.lab = 1.25)
  
  # add points on top
  # these are the points which are considered OOEOR MESS+ve/-ve for the 95% threshold
  points(records_oor_neg_75$longitude, records_oor_neg_75$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_oor_pos_75$longitude, records_oor_pos_75$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within the 75% MESS threshold", "Outside of the 75% MESS threshold"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  mtext(bquote(~italic(.(title))), side = 3, line = -1, outer = TRUE, cex = 2, font = 2)
  
  dev.off()
  
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
  
  }

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

