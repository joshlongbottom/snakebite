# script to get maximal extent of repository data
# clear workspace
rm(list = ls())

pacman::p_load(raster)

# list species
spp_list <- list.files('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data',
                       pattern = '*_trimmed.csv$',
                       full.names = TRUE)

# read in admin 0 shapefile
national <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/admin2013_0.shp')

occ_countries <- NULL

for(i in 1:length(spp_list)){
  
  message(paste('processing file', i, 'of', length(spp_list)), sep = " ")
  
  # get presence records for species
  locations <- read.csv(spp_list[i],
                        stringsAsFactors = FALSE)

  # remove NA records
  locations <- locations[!(is.na(locations$latitude)), ]

  if(nrow(locations) !=0){
  
  # turn into a spatial points dataframe
  lat <- locations$latitude
  lon <- locations$longitude
  
  pts <- cbind(lon,
               lat)
  
  pts_sp <- SpatialPoints(pts,
                          proj4string = CRS(projection(national)))

  # get the ISO codes for all records
  codes <- over(pts_sp, national)$COUNTRY_ID
  
  # get unique country ISO codes
  iso <- unique(codes)
  
  # remove NAs and unassigned countries (XXX)
  iso <- iso[!is.na(iso)]
  iso <- iso[!(iso == 'XXX')]
  
  # get string of all countries for snake distribution
  countries <- paste(iso, collapse = ", ")

  spp_name <- gsub("Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data/", '', spp_list[i])
  spp_name <- gsub("_trimmed.csv", '', spp_name)
  
  # create a dataframe of each country within a species' range
  spec_dist <- as.data.frame(cbind(spp_name, countries))
  
  if(i == 1){
    
    occ_countries <- spec_dist 
    
  } else {
    
    occ_countries <- rbind(occ_countries,
                           spec_dist)
  
    }
  
  
  }
  
}

write.csv(occ_countries,
          'Z:/users/joshua/Snakebite/output/occurrence_data_iso.csv',
          row.names = FALSE)