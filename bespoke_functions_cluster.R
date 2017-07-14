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
  
  country_list_2 <- as.data.frame(strsplit(sub$countries_occ, ","))
  names(country_list_2)[1] <- 'ISO'
  
  country_list <- rbind(country_list,
                        country_list_2)
  
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

the_1000_mess_project <- function(n_boot, in_parallel, n_cores, covs_extract, 
                                  covs, occ_dat, eval_plot, pcc_outpath, variance, var_outpath){
  # where n_boot is the number of bootstraps required, covs_extract is an extracted
  # dataframe including covariate values for each reference point, covs is a stack of covariate 
  # surfaces, occ_dat is the raw reference (occurrence data), in_parallel - if the cluster should
  # be initiated, and n_cores - how many cores to use; loop through and generate x bootstrapped MESS, 
  # and stack the binary outputs of each MESS to get a relative interpolation surface
  
  if(in_parallel){
    
    # initialize the cluster
    registerDoMC(n_cores)
    
    # for however many specified bootstraps, sample reference points with replacement
    sub_sample <- foreach(i = 1:n_boot) %dopar% {
      
      temp_sample <- covs_extract[sample(1:nrow(covs_extract), nrow(covs_extract), replace = TRUE), ]
      
    }
    
    # generate the respective number of mess, using each sub_sample from above
    mess_stack <- foreach(i = 1:length(sub_sample)) %dopar% {

      # generate a mess using the sub_sample
      suppressWarnings(mess_iteration <- mess(covs, sub_sample[[i]], full = TRUE))
      
      # convert to a binary surface
      if('rmess' %in% names(mess_iteration)){
        
        tmp <- mess_iteration[['rmess']] >= 0
        
      } else {
        
        tmp <- mess_iteration[['mess']] >=0
        
      }
      
      # crop to extent
      tmp_masked <- mask(tmp, covs[[1]])
      
      }
    
    # stack the 100 separate MESS  
    raster_stack <- stack(mess_stack)
    
  } else { 
    
    for(i in 1:n_boot){
      
      # sample reference points, with replacement
      sub_sample <- covs_extract[sample(1:nrow(covs_extract), nrow(covs_extract), replace = TRUE), ]
      
      # generate a mess using the sub_sample
      suppressWarnings(mess_iteration <- mess(covs, sub_sample, full = TRUE))
      
      # convert to a binary surface
      if('rmess' %in% names(mess_iteration)){
        
        tmp <- mess_iteration[['rmess']] >= 0
        
      } else {
        
        tmp <- mess_iteration[['mess']] >=0
        
      }
      
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
    
    # subset to get a temporary MESS layer
    temp_MESS <- raster_stack[[i]]
    
    # get an index for occurrence records correctly classified
    vals <- extract(temp_MESS, lat_lon, method='bilinear')
    
    # vals of 1 suggest 100% prediction
    vals <- vals[vals == 1]
    
    rcc_i <- (length(vals) / nrow(occ_dat))
    
    eval_stats <- rbind(eval_stats,
                        rcc_i)
    
    }
  
    # plot distribution of correctly classified across all bootstraps
    eval_stats <- as.data.frame(eval_stats)
    
    title_text <- paste('(n = ', nrow(occ_dat), ' records)', sep = "")
    
    eval_stats$variable <- rep(title_text)
    
    names(eval_stats) <- c('value', 'variable')

    spp_name <- unique(occ_dat$name)
    
    spp_italics <- gsub('_', ' ', spp_name)
    
    x <- suppressWarnings(ggplot(data = eval_stats, mapping = aes(x = value, y = ..scaled..)) + 
                            geom_density(colour = 'cadetblue4', fill = 'cadetblue3') + 
                            ggtitle(bquote(~italic(.(spp_italics))~' MESS evaluation')) +
                            labs(x = "Proportion",
                                 y = "Relative density") +
                            theme(plot.title = element_text(size = 20)) +
                            theme(strip.text = element_text(size = 15)) +
                            theme(axis.title = element_text(size = 15)) +   
                            facet_wrap(~variable, scales = 'fixed')) 
    
    x + scale_x_continuous(breaks = c(seq(0,1,0.1)), limits = c(0, 1)) + scale_y_continuous(breaks = NULL) 
    
    eval_path <- paste(pcc_outpath, spp_name, '_mess_evaluation_plots_', Sys.Date(), '.png', sep = "")
    ggsave(eval_path, width = 600, height = 450, units = 'mm', dpi = 300, device = 'png')
    
  }

  if(variance){
    
    max_cov_frame <- NULL
    min_cov_frame <- NULL
    mean_cov_frame <- NULL
    
    # for however many bootstraps
    for(i in 1:length(sub_sample)){
      
      # generate a max, min, and mean value for each covariate, for each bootstrap
      sub_max <- apply(sub_sample[[i]], 2, max)
      sub_min <- apply(sub_sample[[i]], 2, min)
      sub_mean <- apply(sub_sample[[i]], 2, mean)
      
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
    names(max_cov_melt) <- c('variable', 'max_value')
    
    min_cov_melt <- melt(min_cov_frame)
    min_cov_melt$Var1 <- NULL
    names(min_cov_melt) <- c('variable', 'value')
    min_cov_msd <- melt(min_cov_frame) %>% group_by(group = Var2) %>% summarise(mean_min = mean(value), sd_min = sd(value), median_min = median(value))
    names(min_cov_melt) <- c('variable', 'min_value')
    
    mean_cov_frame <- melt(mean_cov_frame)
    mean_cov_frame$Var1 <- NULL
    names(mean_cov_frame) <- c('variable', 'mean_value')
    
    # plot the distributions of the max and minimum covariate frames
    # plot max
    spp_italics <- gsub('_', ' ', spp_name)
    
    all_cov_stats <- cbind(max_cov_melt,
                           min_cov_melt,
                           mean_cov_frame)
    
    all_cov_stats <- all_cov_stats[c(1, 2, 4, 6)]
    
    all_cov_stats$variable <- gsub('tasselled_cap_brightness_sd', 'Tasselled cap brightness (SD)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('lstday_mean', 'Daytime land surface temperature (Mean)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('lstday_sd', 'Daytime land surface temperature (SD)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('lstnight_mean', 'Night-time land surface temperature (Mean)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('lstnight_sd', 'Night-time land surface temperature (SD)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('tasselled_cap_wetness_mean', 'Tasselled cap wetness (Mean)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('tasselled_cap_wetness_sd', 'Tasselled cap wetness (SD)', all_cov_stats$variable, fixed = TRUE)
    all_cov_stats$variable <- gsub('elevation', 'Elevation', all_cov_stats$variable, fixed = TRUE)
    
    x <- suppressWarnings(ggplot(data = all_cov_stats) +
                          geom_density(aes(x = min_value, y = ..scaled..), 
                                           colour = '#ff7473', fill = '#ff7473') +
                          geom_density(aes(x = mean_value, y = ..scaled..), 
                                           colour = '#ffc952', fill = '#ffc952') +
                          geom_density(aes(x = max_value, y = ..scaled..), 
                                           colour = '#47b8e0', fill = '#47b8e0') +
                          ggtitle(bquote(~italic(.(spp_italics))~'variance (100 bootstraps)')) +
                          labs(x = "Covariate value",
                               y = "Relative density") +
                          theme(plot.title = element_text(size = 20))+
                          theme(strip.text = element_text(size = 13))+
                          theme(axis.title = element_text(size = 15))+
                          facet_wrap(~variable, scales = 'free'))
    
    x + scale_y_continuous(breaks = NULL)
    
    cov_opath <- paste(var_outpath, spp_name, '_cov_distribution_plots_', Sys.Date(), '.png', sep = "")
    ggsave(cov_opath, width = 600, height = 450, units = 'mm', dpi = 300, device = 'png')

    }
  
  return(monster_mess)
 
}

bufferMESSpositives <- function (range, coords, radius, mess, admin_0, outpath, spp_name) {
  # function to generate a buffer of a given radius around the locations of some points, and merge
  # with an EOR polygon, where:
  # `range` is a SpatialPolygons object of length 1 (polygon of current EOR)
  # `coords` is a two-column matrix of coordinates (x then y; points requiring a buffer)
  # `radius` is a single positive number giving the radius of the circle to generate around each point (in km)
  # `mess` is a binary multivariate environmental similarity surface (MESS), which have values of '1'
  # for areas of interpolation, and '0' for areas of extrapolation
  # `admin_0` is a global FAO admin 0 shapefile, with an iso3 attribute named "COUNTRY_ID"
  # `outpath` is the path to the directory where the output should be saved
  # `spp_name` is a vector giving the name of the species "Genus_species"
  
  # packages needed
  require(sp)
  require(rgeos)
  
  # check things
  stopifnot(inherits(range, 'SpatialPolygons'))
  stopifnot(inherits(coords, c('matrix', 'data.frame')))
  stopifnot(ncol(coords) == 2)
  stopifnot(radius > 0)
  
  # turn coords into SpatialPoints
  spts <- SpatialPoints(coords,
                        proj4string = CRS(proj4string(range)))

  # convert radius to decimal degrees at the equator
  radius <- radius/111.32
  
  # turn SpatialPoints in SpatialPolygons circles
  suppressWarnings(spts_buffer <- gBuffer(spts, width = radius))
  
  # clip buffered points by the MESS surface (caution that this is at a 5 km resolution)
  # first change 0 binary values to NA in the mess object
  mess[mess == 0] <- NA
  
  # then convert raster into a polygon, and merge points
  mess <- rasterToPolygons(mess, dissolve = FALSE)
  
  # change layer name for each mess
  names(mess@data) <- 'layer'
  
  # get list of unique values to merge on, use these to reduce size of polygon
  ids <- mess@data$layer
  mess <- unionSpatialPolygons(mess, ids)
  
  # fix potential self-intersection issue with buffered points
  suppressWarnings(spts_buffer <- gBuffer(spts_buffer, byid = TRUE, width = 0))
  # then clip the buffered points by the limited mess
  spts_buffer_clip <- gIntersection(spts_buffer, mess)
  
  # get ISO for each point
  polys <- spts_buffer_clip@polygons[[1]]@Polygons
  
  # create a reference for each polygon within the buffered spatialPolys object
  pl <- vector("list", length(polys))
  
  for (i in 1:length(polys)) { 
    
    pl[i] <- Polygons(list(polys[[i]]), i) 
    
  }
  
  buffer_clip_spolys <- SpatialPolygons(pl)
  
  # generate row ids
  row.ids = sapply(slot(buffer_clip_spolys, "polygons"), function(i) slot(i, "ID"))
  
  # project 
  proj4string(buffer_clip_spolys) <- CRS(proj4string(range))
  
  # get the centroids of each polygon, so we can determine which country they land in
  suppressWarnings(centroids <- getSpPPolygonsLabptSlots(buffer_clip_spolys))
  
  # change centroid coordinates to a spatial points object
  centroids <- SpatialPoints(centroids,
                             proj4string = CRS(projection(range)))
  
  # get the ISO code for all of the polygon centroids
  iso <- over(centroids, admin_0)$COUNTRY_ID
  
  # get a list of the medical classification of each country in the EOR
  class_dataframe <- as.data.frame(range)
  lookup <- match(iso, class_dataframe$COUNTRY_ID)
  med_class <- class_dataframe$Med_Class[lookup]
  
  # if an ISO isn't in the current EOR, generate a medical class:
  # if the rest of the species is considered class 2, assign class 2
  # if there's both C1 and C2, assign as being C1
  class_check <- unique(class_dataframe$Med_Class)
  
  if(length(class_check) == 2){
    
    temp_class <- 1
    
  } else {
    
    temp_class <- class_check
  
  }
  
  med_class[is.na(med_class)] <- temp_class

  # create a dataframe, to allow coersion from a SpatialPolygons object to SpatialPolygonsDataframe
  buffer_dataframe <- data.frame(FID = 1:length(buffer_clip_spolys))
  buffer_dataframe <- data.frame(COUNTRY_ID = iso,
                                 ADMN_LEVEL = 0, 
                                 Med_Class = med_class)
  
  # merge to form a spatialPolygonsDataframe
  buffer_spdf <- SpatialPolygonsDataFrame(buffer_clip_spolys, buffer_dataframe)
  buffer_spdf <- spChFIDs(buffer_spdf, paste("buffer", row.names(buffer_spdf), sep="."))
  
  variable_list <- c('COUNTRY_ID',
                     'Med_Class',
                     'ADMN_LEVEL')

  # subset both shapefiles to just have the necessary columns
  raw_sub <- range[,(names(range) %in% variable_list)]
  buf_sub <- buffer_spdf[,(names(buffer_spdf) %in% variable_list)]
  
  # merge buffer, and raw shapefile to create one shapefile which has multiple polygons (keeping med class)
  merged_range_mp <- rbind(raw_sub, buf_sub)
  merged_range_co <- aggregate(merged_range_mp, c("COUNTRY_ID", "Med_Class", "ADMN_LEVEL"))
  
  # write out the new shapefile, define path
  shp_path <- paste(outpath, spp_name, sep = "")
  shapefile(merged_range_co, shp_path, overwrite = TRUE)
  
  # combine points to create a single range polygon
  modified_range <- gUnion(range, spts_buffer_clip)

  # return expanded MESS polygon
  return (modified_range)
  
}
