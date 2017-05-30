# load bespoke functions

vector.is.empty <- function(data) {
  # return a logical as to whether the vector is empty (TRUE), or not (FALSE)
  return(length(data) == 0)
                   
}

prepare_covariates <- function(cov_list, ext, range) {
  # where 'cov_list' is a dataframe with covariate names and covariate file paths,
  # ext is a shapefile of the study extent, and range is the EOR shapefile - loop through the
  # covariate list, load and stack the rasters, and mask by a surface with all NA cell values
  
  # create a list of the covariates to loop through and open
  covars <- as.list(cov_list$cov_path)

  # create list of names to rename the covariates
  names <- as.list(cov_list$cov_name)

  # loop through covariate list, opening links to the rasters
  covs <- lapply(covars, raster)

  # crop all of the raster layers by the extent of study extent
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
  
  # project covariates in the same projection as the EOR map
  proj4string(covs) <- proj4string(range)
  
  # return prepped covariates
  return(covs)
  
}

prepare_eor <- function(sub) {
  # where sub is a dataframe containing a filepath for a speies, read in EOR shapefile
  # and if there's more than one polygon within the shapefile, dissolve into one single
  # polygon
  
  # read in EOR shapefile for each species
  range <- shapefile(sub$shapefile_path)
  
  # dissolve polygons
  n_polys <- as.numeric(nrow(range@data))
  ids <- range@data$ADMN_LEVEL  
  
  if(n_polys > 1){
    
    suppressWarnings(range <- gBuffer(range, byid = TRUE, width = 0))
    suppressWarnings(range <- unionSpatialPolygons(range, ids))
    
  }
  
  return(range)

}

generate_ext <- function(sub, admin_0){
  # where sub is a dataframe containing a list of ISO codes within a species' extent, and 
  # admin_0 is a shapefile with national country polygons, generate a national extent for 
  # the species
  
  # get a list of iso3 codes
  country_list <- as.data.frame(strsplit(sub$countries_raw, ","))
  names(country_list)[1] <- 'ISO'
  
  country_list <- as.character(unique(country_list$ISO))
  country_list <- gsub(" ", "", country_list)
  
  # remove any disputed regions (iso XXX)
  country_list <- country_list[!(country_list == 'XXX')]
  
  # create extent shapefile
  ext <- admin_0[admin_0$COUNTRY_ID %in% country_list, ]
  
  return(ext)
  
}

load_occurrence <- function(spp_name, range){
  # where spp_name is a vector with the name of the species (genus_species), and 
  # range is an EOR shapefile for each species, load in occurrence data, turn into a 
  # spatial points dataframe and project in the same coordinate system as the EOR
  
  # get presence records for species
  dat_path <- paste('data/raw/species_data/', 
                  spp_name,
                  '_trimmed.csv', 
                  sep = "")

  locations <- read.csv(dat_path,
                        stringsAsFactors = FALSE)

  # remove NA records
  locations <- locations[!(is.na(locations$latitude)), ]
  
  if(nrow(locations) !=0){
    
  # turn into a spatial points dataframe
  coordinates(locations) <- c("longitude", "latitude")

  # project coordinates in same projection as shape
  proj4string(locations) <- proj4string(range)
  
  }
  
  return(locations)
  
}

measure_variance <- function(n_boot, covs_extract_df, distribution_path, spp_name){
  # where n_boot is the number of bootstraps to perform on the covariate data,
  # covs_extract_df is the dataframe of covariate values for each reference point, 
  # distribution path is a vector with the directory in which plots should be saved, and
  # spp_name is the name of the species; measure the variance within each of the reference
  # data points, and plot a distribution of the maximum, mean, and medium values for each 
  # each covariate, across the specified number of boostraps (n_boot)
  max_cov_frame <- NULL
  min_cov_frame <- NULL
  mean_cov_frame <- NULL
  
  # for however many bootstraps
  for(i in 1:n_boot){
    
    # subsample the covariate data
    sub_sample <- covs_extract_df[sample(1:nrow(covs_extract_df), nrow(covs_extract_df), replace = TRUE), ]
    
    # generate a max, min, and mean value for each covariate, for each bootstrap
    sub_max <- apply(sub_sample, 2, max)
    sub_min <- apply(sub_sample, 2, min)
    sub_mean <- apply(sub_sample, 2, mean)
    
    # add these values to a dataframe
    max_cov_frame <- rbind(max_cov_frame,
                           sub_max)
    
    min_cov_frame <- rbind(min_cov_frame,
                           sub_min)
    
    mean_cov_frame <- rbind(mean_cov_frame,
                            sub_mean)
    
  }
  
  # melt the covariate max, min, mean dataframes, so they can be used to generate density plots
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

  # plot the distributions of the max and minimum covariate frames
  # plot max
  spp_italics <- gsub('_', ' ', spp_name)

  suppressWarnings(ggplot(data = max_cov_melt, mapping = aes(x = value)) + 
                     geom_density(colour = 'cadetblue4', fill = 'cadetblue3') +  
                     ggtitle(bquote(~italic(.(spp_italics))~' maximum covariate values')) +
                     facet_wrap(~variable, scales = 'free'))
  
  max_path <- paste(distribution_path, spp_name, '_maximum_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
  ggsave(max_path, width = 600, height = 450, units = 'mm', dpi = 300)
  
  # plot mean
  suppressWarnings(ggplot(data = mean_cov_frame, mapping = aes(x = value)) + 
                     geom_density(colour = 'cadetblue4', fill = 'cadetblue3') + 
                     ggtitle(bquote(~italic(.(spp_italics))~' mean covariate values')) +   
                     facet_wrap(~variable, scales = 'free')) 
  
  mean_path <- paste(distribution_path, spp_name, '_mean_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
  ggsave(mean_path, width = 600, height = 450, units = 'mm', dpi = 300)
  
  # plot min
  suppressWarnings(ggplot(data = min_cov_melt, mapping = aes(x = value)) + 
                     geom_density(colour = 'cadetblue4', fill = 'cadetblue3') +
                     ggtitle(bquote(~italic(.(spp_italics))~' minimum covariate values')) +   
                     facet_wrap(~variable, scales = 'free'))
  
  min_path <- paste(distribution_path, spp_name, '_minimum_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
  ggsave(min_path, width = 600, height = 450, units = 'mm', dpi = 300)
  
  # generate a two row dataframe, which only contains the median minimum value from all bootstraps, and the 
  # median maximum value from all bootstraps (extremely conservative)
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
  
  return(covariate_stats)
  
}

the_1000_mess_project <- function(n_boot, in_parallel, n_cores, covs_extract, covs, occ_dat, eval_plot, plot_outpath){
  # where n_boot is the number of bootstraps required, covs_extract is an extracted
  # dataframe including covariate values for each reference point, covs is a stack of covariate 
  # surfaces, occ_dat is the raw reference (occurrence data), in_parallel - if the cluster should
  # be initiated, and n_cores - how many cores to use; loop through and generate x bootstrapped MESS, 
  # and stack the binary outputs of each MESS to get a relative interpolation surface
  
  if(in_parallel){
    
    # initialize the cluster
    registerDoMC(n_cores)
    
    # for however many specified bootstraps, generate the respective number of mess
    mess_stack <- foreach(i=1:n_boot) %dopar% {
      
      # sample reference points, with replacement
      sub_sample <- covs_extract[sample(1:nrow(covs_extract), nrow(covs_extract), replace = TRUE), ]
      
      # generate a mess using the sub_sample
      suppressWarnings(mess_iteration <- mess(covs, sub_sample, full = TRUE))
      
      # convert to a binary surface
      tmp <- mess_iteration[['rmess']] >= 0
      
      # crop to extent
      tmp_masked <- mask(tmp, covs[[1]])
      
      }
      
      raster_stack <- stack(mess_stack)
    
  } else { 
    
    for(i in 1:n_boot){
      
      # sample reference points, with replacement
      sub_sample <- covs_extract[sample(1:nrow(covs_extract), nrow(covs_extract), replace = TRUE), ]
      
      # generate a mess using the sub_sample
      suppressWarnings(mess_iteration <- mess(covs, sub_sample, full = TRUE))
      
      # convert to a binary surface
      tmp <- mess_iteration[['rmess']] >= 0
      
      # crop to extent
      tmp_masked <- mask(tmp, covs[[1]])
      
      # stack with previous iteration
      if(i == 1){
        
        raster_stack <- tmp_masked
        
      } else {
        
        raster_stack <- stack(raster_stack,
                              tmp_masked)
        
      }
      
    }
    
  }
  
  # sum cell values to generate a relative interpolation surface
  monster_mess <- sum(raster_stack, na.rm = TRUE)
  
  # crop monster mess by the effective land/sea mask
  monster_mess <- mask(monster_mess, covs[[1]])
  
  if(eval_plot){
  
    ## get stats on the performance of each MESS
    # turn occ_dat into a spatial points dataframe
    coordinates(occ_dat) <- c("longitude", "latitude")
  
    # project coordinates in same projection as covariates
    proj4string(occ_dat) <- proj4string(covs[[1]])
  
    lat <- occ_dat$latitude
    lon <- occ_dat$longitude
  
    lat_lon <- cbind(lon,
                     lat)
  
    eval_stats <- NULL
  
    for(i in 1:nlayers(raster_stack)){
  
    # get an index for occurrence records correctly classified
    occ_dat$mess <- cellFromXY(raster_stack[[i]], lat_lon)
    
    # get count of number correctly classified (land within a cell value of the MESS)
    rcc <- as.data.frame(occ_dat[!(is.na(occ_dat$mess)), ])
    
    rcc_i <- (nrow(rcc) / nrow(occ_dat))
    
    eval_stats <- rbind(eval_stats,
                        rcc_i)
    
    }
  
    # plot distribution of correctly classified across all bootstraps
    eval_stats <- as.data.frame(eval_stats)
    
    eval_stats$variable <- rep('proportion correctly classified')
    
    names(eval_stats) <- c('value', 'variable')

    spp_name <- unique(occ_dat$name)
    
    spp_italics <- gsub('_', ' ', spp_name)
    
    suppressWarnings(ggplot(data = eval_stats, mapping = aes(x = value)) + 
                     geom_density(colour = 'cadetblue4', fill = 'cadetblue3') + 
                     ggtitle(bquote(~italic(.(spp_italics))~' MESS evaluation (proportion correctly classified)')) +   
                     facet_wrap(~variable, scales = 'free')) 
    
    eval_path <- paste(plot_outpath, spp_name, '_mess_evaluation_plots_', Sys.Date(), '.png', sep = "")
    ggsave(eval_path, device = "png")
    
  }
  
  return(monster_mess)
 
}

plot_mess <- function(png_mess, spp_name, add_points, tmp_masked, bootstrapped_mess){
  
  # create a plotting window to plot both of the surfaces
  png_name <- paste(png_mess, spp_name, '_species_mess_maps_', Sys.Date(), '.png', sep = "")
  
  if(add_points == TRUE){
    
  png(png_name,
      width = 450,
      height = 400,
      units = 'mm',
      res = 300)
  par(mfrow = c(2, 2))
  
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
  
  # add legend to the bottom
  legend('bottomleft', c("Interpolation","Extrapolation"), 
         pch = c(15, 15),
         col = c("springgreen4","gainsboro"), bty = 'n')
  
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
  
  # now plot again, adding the points to each plot
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
  
  # add points on top
  points(records_outside$longitude, records_outside$latitude, pch = 20, cex = 0.75, col = '#D93529')
  points(records_inside$longitude, records_inside$latitude, pch = 20, cex = 0.75, col = 'blue')
  
  legend('bottomleft', c("Interpolation","Extrapolation", "Within range", "Outside range"), 
         pch = c(15, 15, 20, 20),
         col = c("springgreen4","gainsboro", "blue", "#D93529"), bty = 'n')
  
  dev.off()
  
  } else {
    
    # just plot the two without the points
    png(png_name,
        width = 450,
        height = 200,
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
    
    # add legend to the bottom
    legend('bottomleft', c("Interpolation","Extrapolation"), 
           pch = c(15, 15),
           col = c("springgreen4","gainsboro"), bty = 'n')
    
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
  
  # create a return
  status <- 'Plots complete'
  
  return(status)
  
}
