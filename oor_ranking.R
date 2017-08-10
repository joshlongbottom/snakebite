# clear workspace
rm(list = ls ())

# load packages
pacman::p_load(raster, dismo, rgeos, sp, spdep)

# load functions
source('Z:/users/joshua/Snakebite/snakebite/bespoke_functions_cluster.R')

# read in snakelist
snake_list <- read.csv('Z:/users/joshua/Snakebite/snakebite/snake_list.csv',
                       stringsAsFactors = FALSE)

# list species for which we have occurrence data
spp_list <- list.files('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data/',
                       pattern = '*_raw.csv$',
                       full.names = FALSE)

combined <- NULL

for(i in 1:length(spp_list)){
  
  # get species info from spp_list
  spp_name <- as.character(spp_list[i])
  spp_name <- gsub('_raw.csv', '', spp_name)
  
  # inform progress
  message(paste('Processing species ', i, ' of ', length(spp_list), ' (', spp_name, ') ', Sys.time(), sep = ""))
  
  sub <- snake_list[snake_list$split_spp == spp_name, ]
  
  if(spp_name == 'Lachesis_stenophrys'){
    
    sub <- snake_list[snake_list$split_spp == 'Lachesis_stenophyrs', ]
    
  }
  
  data_out <- NULL
  
  # read in EOR shapefile for each species, and dissolve if >1 polygon
  raw_range <- shapefile(sub$shapefile_path)
  
  # get presence records for species
  dat_path <- paste0('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data/', 
                    spp_name,
                    '_trimmed.csv')
  
  locations <- read.csv(dat_path,
                        stringsAsFactors = FALSE)
  
  # remove NA records
  locations <- locations[!(is.na(locations$latitude)), ]
  
  if(nrow(locations) !=0){
    
    # turn into a spatial points dataframe
    coordinates(locations) <- c("longitude", "latitude")
    
    # project coordinates in same projection as shape
    proj4string(locations) <- proj4string(raw_range)
    
  
  
  distances <- apply(gDistance(locations, raw_range,byid = TRUE), 2, min)
  
  in_range <- distances[distances == 0]
  
  data_out <- data.frame(species = spp_name,
                         number_of_records = length(distances),
                         number_in_range = length(in_range),
                         number_out_of_range = (length(distances) - length(in_range)),
                         cumulative_distance = sum(distances))
  
  }
  
  if(!(vector.is.empty(data_out))){
    
    combined <- rbind(combined,
                      data_out)
    
  }
}

outpath <- paste0('Z:/users/joshua/Snakebite/output/species_occurrence_plots/oor_stats/combined_', Sys.Date(), '.csv')

write.csv(combined,
          outpath,
          row.names = FALSE)