# clear workspace
rm(list = ls ())

# load packages
pacman::p_load(raster, dismo, rgeos, seegSDM, maptools, ggplot2, colorRamps, reshape2, dplyr)

# set working directory
setwd('Z:/users/joshua/Snakebite/')

# load function(s)
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

pdf_name <- paste('output/species_occurrence_plots/species_occurrence_plots_', Sys.Date(), '.pdf', sep = "")

pdf(pdf_name,                    
    width = 8.27,
    height = 11.29)
par(mfrow = c(3, 2))

# loop through and create stacks for each extent
# for(i in 1:length(spp_list)){
for(i in 1:1){
  
  # get country list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('_raw.csv', '', spp_name)

  sub <- snake_list[snake_list$split_spp == spp_name, ]
  
  message(paste('processing file', i, 'of length', length(spp_list), sep = " "))
  
  # create extent per species
  country_list <- as.data.frame(strsplit(sub$countries_raw, ","))
  names(country_list)[1] <- 'ISO'
  
  country_list <- as.character(unique(country_list$ISO))
  country_list <- gsub(" ", "", country_list)
  
  america_list <- c('Agkistrodon_contortrix', 'Agkistrodon_piscivorus', 'Crotalus_atrox', 
                    'Crotalus_molossus', 'Crotalus_oreganus', 'Crotalus_ruber', 'Crotalus_scutulatus',
                    'Crotalus_spp', 'Crotalus_viridis', 'Micruroides_euryxanthus', 'Micrurus_fulvius', 
                    'Micrurus_tener', 'Sistrurus_catenatus')
  
  ext <- admin_0[admin_0$COUNTRY_ID %in% country_list, ]
  
  if(spp_name %in% america_list){
    
    ext <- shapefile('World shapefiles/America.shp')
    
  }
  
  # read in EOR for each species
  range <- shapefile(sub$shapefile_path)
  
  # dissolve polygons
  n_polys <- as.numeric(nrow(range@data))
  ids <- range@data$ADMN_LEVEL  
  
  if(n_polys > 1){
    
    suppressWarnings(range <- gBuffer(range, byid = TRUE, width = 0))
    suppressWarnings(range <- unionSpatialPolygons(range, ids))
    
  }
  
  # organise covariates for analysis
  # create a list to loop through and open
  covars <- as.list(cov_list$cov_path)
            
  # create list of names to rename covs by
  names <- as.list(cov_list$cov_name)

  # loop through opening links to the rasters
  covs <- lapply(covars, raster)

  # crop raster layers by the extent of study extent
  covs <- lapply(covs, crop, ext)

  # stack the rasters
  covs <- stack(covs) 

  # loop through covariate layers, find any NAs,
  # find the cov with most NAs, and mask the whole stack by this layer
  master_mask <- masterMask(covs)

  # mask all covariates by master mask
  covs <- mask(covs, master_mask)

  # change the names of the covariates
  names(covs) <- names
  
  # get presence records for species
  dat_path <- paste('output/species_occurrence_plots/species data/', 
                    spp_name,
                    '_trimmed.csv', 
                    sep = "")
  
  locations <- read.csv(dat_path,
                        stringsAsFactors = FALSE)

  # remove NA records
  locations <- locations[!(is.na(locations$latitude)), ]
  
  # turn into a spatial points dataframe
  coordinates(locations) <- c("longitude", "latitude")
  
  # project coordinates in same projection as shape
  proj4string(locations) <- proj4string(range)
  proj4string(covs) <- proj4string(range)
  
  # get an index for occurrence records within species' range
  records_inside <- !is.na(over(locations, as(range, "SpatialPolygons")))
  
  records_outside <- locations[records_inside == FALSE, ]
  records_inside <- locations[records_inside == TRUE, ]
  records_outside <- as.data.frame(records_outside)
  records_inside <- as.data.frame(records_inside)
  
  # generate a MESS, using a prerequisite of a minimum of 5 data points within the species range
  if(nrow(records_inside) >= 5){

  # extract covariate data for all of that species' presence records, within the WHO EOR
  covs_extract <- extract(covs, records_inside[c("longitude", "latitude")])
  covs_extract_df <- as.data.frame(covs_extract)
  
  # subsample from this extracted covariate data, to generate 100 bootstraps
  # then plot the distribution of the max and mix values for each covariate
  # to measure variation in the input data
  max_cov_frame <- NULL
  min_cov_frame <- NULL
  mean_cov_frame <- NULL
  
  for(i in 1:1000){
  
    sub <- covs_extract_df[sample(1:nrow(covs_extract_df), 30, replace = TRUE), ]
  
    sub_max <- apply(sub, 2, max)
    sub_min <- apply(sub, 2, min)
    sub_mean <- apply(sub, 2, mean)
    
    max_cov_frame <- rbind(max_cov_frame,
                           sub_max)
    
    min_cov_frame <- rbind(min_cov_frame,
                           sub_min)
    
    mean_cov_frame <- rbind(mean_cov_frame,
                            sub_mean)
    
  }
  
  # max_cov_melt <- melt(max_cov_frame)
  # max_cov_melt$Var1 <- NULL
  # names(max_cov_melt) <- c('variable', 'value')
  # max_cov_msd <- melt(max_cov_frame) %>% group_by(group = Var2) %>% summarise(mean = mean(value), sd = sd(value))
  # max_cov_den <- density(max_cov_melt)
  
  max_cov_melt <- melt(max_cov_frame)
  max_cov_melt$Var1 <- NULL
  names(max_cov_melt) <- c('variable', 'value')
  max_cov_msd <- melt(max_cov_frame) %>% group_by(group = Var2) %>% summarise(mean_max = mean(value), sd_max = sd(value), median_max = median(value))
  
  min_cov_melt <- melt(min_cov_frame)
  min_cov_melt$Var1 <- NULL
  names(min_cov_melt) <- c('variable', 'value')
  min_cov_msd <- melt(min_cov_frame) %>% group_by(group = Var2) %>% summarise(mean_min = mean(value), sd_min = sd(value), median_min = median(value))
  
  mean_cov_frame <- melt(mean_cov_frame)
  mean_cov_frame$Var1 <- NULL
  names(mean_cov_frame) <- c('variable', 'value')
  
  min_cov_median <- as.data.frame(t(min_cov_msd[c(1, 4)]))
  colnames(min_cov_median) <- as.character(unlist(min_cov_median[1,]))
  min_cov_median = min_cov_median[-1, ]
  
  max_cov_median <- as.data.frame(t(max_cov_msd[c(1, 4)]))
  colnames(max_cov_median) <- as.character(unlist(max_cov_median[1,]))
  max_cov_median = max_cov_median[-1, ]
  
  covariate_stats <- rbind(min_cov_median,
                           max_cov_median)
  
  row.names(covariate_stats) <- c('1', '2')
  
  suppressWarnings(covariate_stats <- data.frame(lapply(covariate_stats, as.character), stringsAsFactors = FALSE))
  suppressWarnings(covariate_stats <- data.frame(lapply(covariate_stats, as.numeric), stringsAsFactors = FALSE))
  
  # plot the distributions of the max and minimum covariate frames
  # plot max
  suppressWarnings(ggplot(data = max_cov_melt, mapping = aes(x = value)) + 
                   geom_density(colour = 'cadetblue4', fill = 'cadetblue3') +  
                   ggtitle("Max covariate values") +
                   facet_wrap(~variable, scales = 'free'))
  
  max_path <- paste('output/species_mess_maps/distribution plots/', spp_name, '_maximum_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
  ggsave(max_path, width = 600, height = 450, units = 'mm', dpi = 300)
  
  # plot mean
  suppressWarnings(ggplot(data = mean_cov_frame, mapping = aes(x = value)) + 
                   geom_density(colour = 'cadetblue4', fill = 'cadetblue3') + 
                   ggtitle("Mean covariate values") +   
                   facet_wrap(~variable, scales = 'free')) 
  
  mean_path <- paste('output/species_mess_maps/distribution plots/', spp_name, '_mean_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
  ggsave(mean_path, width = 600, height = 450, units = 'mm', dpi = 300)
  
  # plot min
  suppressWarnings(ggplot(data = min_cov_melt, mapping = aes(x = value)) + 
                   geom_density(colour = 'cadetblue4', fill = 'cadetblue3') +
                   ggtitle("Min covariate values") +   
                   facet_wrap(~variable, scales = 'free'))

  min_path <- paste('output/species_mess_maps/distribution plots/', spp_name, '_minimum_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
  ggsave(min_path, width = 600, height = 450, units = 'mm', dpi = 300)
  
  # create mess map (mess function)
  mess_map <- mess(covs, covariate_stats, full = TRUE)

  # make binary, interpolation/extrapolation surface and save
  tmp <- mess_map[['rmess']] >= 0
  raw_mess <- mess_map[['rmess']]
  
  # crop to extent
  tmp_masked <- mask(tmp, master_mask)
  raw_mess_masked <- mask(raw_mess, master_mask)
  
  # png_name <- paste('output/species_mess_maps/', spp_name, '_species_mess_map_', Sys.Date(), '.png', sep = "")
  # 
  # png(png_name,                    
  #     width = 450,
  #     height = 250,
  #     units = 'mm',
  #     res = 300)
  # par(mfrow = c(1, 2))
  
  # plot MESS
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
  
  title(xlab = 'Binary MESS', line = 0)
  
  # add points on top
  points(records_outside$longitude, records_outside$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_inside$longitude, records_inside$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within range", "Outside range"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  # plot raw MESS
  # define a colour ramp
  # col1 <- colorRampPalette(c('turquoise3', 'slateblue4'))
  # colour_levels = 10
  # 
  # # change inf values to NA
  # is.na(raw_mess_masked) <- sapply(raw_mess_masked, is.infinite)
  # 
  # # get absolute min and max of raster values
  # max_abolute_value = max(abs(c(cellStats(raw_mess_masked, min), cellStats(raw_mess_masked, max))))
  # 
  # # set a colour sequence
  # colour_sequence = seq(-max_abolute_value, max_abolute_value, length.out = colour_levels + 1)
  # 
  # n_in_class = hist(raw_mess_masked, breaks = colour_sequence, plot = FALSE)$counts>0
  # col_to_include = min(which(n_in_class == TRUE)):max(which(n_in_class == TRUE))
  # breaks_to_include = min(which(n_in_class == TRUE)):(max(which(n_in_class == TRUE))+1)
  # 
  # plot(raw_mess_masked, 
  #      main = bquote(~italic(.(title))), 
  #      legend = TRUE,
  #      axes = FALSE, 
  #      box = FALSE,
  #      col = col1(n = colour_levels)[col_to_include],
  #      breaks = colour_sequence[breaks_to_include])

  plot(raw_mess_masked,
       main = bquote(~italic(.(title))),
       legend = FALSE,
       axes = FALSE,
       box = FALSE)
   plot(range,
      add = TRUE,
      border = 'black',
      lty = 1,
      lwd = 1)
  
  title(xlab = 'Raw MESS', line = 0)
  
  # write to disk
  outpath <- paste('output/species_mess_maps/files/', spp_name, sep = '')
  
  writeRaster(tmp_masked, 
              file = outpath, 
              format = 'GTiff', 
              overwrite = TRUE)
  
  outpath_2 <- paste('output/species_mess_maps/files/', spp_name, '_raw_MESS', sep = '')
  
  writeRaster(raw_mess_masked, 
              file = outpath_2, 
              format = 'GTiff', 
              overwrite = TRUE)
  
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
dev.off()
 
