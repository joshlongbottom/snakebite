#clear workspace
rm(list = ls ())

# load packages
library(raster)
library(dismo)
library(rgeos)
library(seegSDM)
library(maptools)
library(ggplot2)

# read in a list of the covariates to be used in the model
cov_list <- read.csv('Z:/users/joshua/Snakebite/rasters/covariates_for_model.csv',
                     stringsAsFactors = FALSE)

# read in snakelist
snake_list <- read.csv('Z:/users/joshua/Snakebite/snakebite/snake_list.csv',
                       stringsAsFactors = FALSE)

# read in admin 0
admin_0 <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/admin2013_0.shp')

# list species for which we have occurrence data
spp_list <- list.files('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data/',
                       pattern = '*csv$',
                       full.names = FALSE)

pdf('Z:/users/joshua/Snakebite/output/species_mess_maps/species_mess_maps_Y2017M04D26v2.pdf',                    
    width = 8.27,
    height = 11.29)
par(mfrow = c(3, 2))

# loop through and create stacks for each extent
for(i in 1:length(spp_list)){

  # get country list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('.csv', '', spp_name)

  sub <- snake_list[snake_list$split_spp == spp_name, ]
  
  message(paste('processing file', i, 'of length', length(spp_list), sep = " "))
  
  # create extent per species
  country_list <- as.data.frame(strsplit(sub$countries_raw, ","))
  names(country_list)[1] <- 'ISO'
  
  country_list <- as.character(unique(country_list$ISO))
  country_list <- gsub(" ", "", country_list)
  
  america_list <- c('Agkistrodon_contortrix', 'Agkistrodon_piscivorus')
  
  ext <- admin_0[admin_0$COUNTRY_ID %in% country_list, ]
  
  if(spp_name %in% america_list){
    
    ext <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/America.shp')
    
  }
  
  # read in EOR for each species
  range <- shapefile(sub$shapefile_path)
  
  # dissolve polygons
  n_polys <- as.numeric(nrow(range@data))
  ids <- range@data$ADMN_LEVEL  
  
  if(n_polys > 1){
    
    range <- gBuffer(range, byid = TRUE, width = 0)
    range <- unionSpatialPolygons(range, ids)
    
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
  dat_path <- paste('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data/', 
                    spp_name, 
                    '.csv', 
                    sep = "")
  
  dat <- read.csv(dat_path,
                  stringsAsFactors = FALSE)
  
  locations <- unique(dat[c('lat', 'lon')]) 
  
  # extract covariate data for all of that species' presence records
  covs_extract <- extract(covs, locations[c('lon', 'lat')])

  # create mess map (mess function)
  mess_map <- mess(covs, covs_extract, full = TRUE)

  # make binary, interpolation/extrapolation surface and save
  tmp <- mess_map[['rmess']] >= 0
  raw_mess <- mess_map[['rmess']]
  
  # crop to extent
  tmp_masked <- mask(tmp, master_mask)
  raw_mess_masked <- mask(raw_mess, master_mask)
  
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
       lwd = 0.5)
  
  title(xlab = 'Binary MESS', line = 0)
  
  # add points on top
  points(locations$lon, locations$lat, pch = 20, cex = 0.2, col = '#D93529')

  plot(raw_mess_masked, 
       main = bquote(~italic(.(title))), 
       legend = TRUE,
       axes = FALSE, 
       box = FALSE)
  plot(range,
       add = TRUE,
       border = 'black',
       lty = 1,
       lwd = 0.5)
  # plot(raw_mess_masked, 
  #      horizontal = TRUE,
  #      add = TRUE,
  #      smallPlot = c(.15, .5, .84, .86), 
  #      legend.only = TRUE)
  
  title(xlab = 'Raw MESS', line = 0)
  
  # add points on top
  points(locations$lon, locations$lat, pch = 20, cex = 0.2, col = '#D93529')
  
  # write to disk
  outpath <- paste('Z:/users/joshua/Snakebite/output/species_mess_maps/files/', spp_name, sep = '')
  writeRaster(tmp_masked, 
              file = outpath, 
              format = 'GTiff', 
              overwrite = TRUE)
}

dev.off() 
