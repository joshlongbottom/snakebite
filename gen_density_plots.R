# density plots (investigation in to a ppp)
# clear workspace
rm(list = ls())

# load required packages
require(spatstat)
require(dismo)

# load function
vector.is.empty <- function(x) return(length(x) == 0)

# read in snake species list
species_sheet <- read.csv('Z:/users/joshua/Snakebite/snakebite/snake_list.csv',
                          stringsAsFactors = FALSE)

# get a list of the unique snake species 
species_list <- unique(species_sheet$split_spp)
species_list <- as.list(species_list[!(species_list == "")])

# for every species in the list, grab occurrence records from gbif, and 
# generate a ppp using data

# loop through and populate pdf with species occurrence plots
for(i in 1:length(species_list)){
  
  # specify species
  species <- species_list[i]
  
  # subset species sheet (to obtain file path to shapefile)
  sub <- species_sheet[species_sheet$split_spp == species, ]
  
  # define species name    
  title <- sub$split_spp    
  
  # split the species into 'genus' and 'spp'
  split <- as.data.frame(strsplit(title, "_"))
  genus <- as.character(split[rownames(split) == '1', ])
  
  if(nrow(split) == 2){
    
    spp <- as.character(split[rownames(split) == '2', ])
  
    }
  
  if(nrow(split) == 3){
    
    temp_1 <- as.character(split[rownames(split) == '2', ])
    temp_2 <- as.character(split[rownames(split) == '3', ])
    spp <- paste(temp_1, temp_2)
  
    }
  
  # get gbif presence for that species
  # if the shapefile is for a group of species, only get genus results from gbif
  if(spp != 'spp'){
    
    spp_data <- gbif(genus, spp, end = 100000)
    
  } else {
    
    spp_data <- gbif(genus, '', end = 100000)  
    
  }
  
  # if spp_data is not an empty dataframe get unique lat/longs for the species
  if(vector.is.empty(spp_data) == FALSE) {
    
    if('lat' %in% colnames(spp_data)){
      
      locations <- unique(spp_data[c('lat', 'lon')])
      
      # if there are records, remove any NAs
      locations <- locations[!is.na(locations$lat), ]
      
      locations <- locations[!is.na(locations$lon), ]
      
    } else {
      
      locations <- NULL
      
    }
  }
  
  # if there are records, gen an observation window
  # define shapefile for that species
  shape_path <- sub$shapefile_path
  
  # read in shapefile
  shape <- shapefile(shape_path)
  extents <- extent(shape)
  latitudes <- c(extents@ymin, extents@ymax)
  longitudes <- c(extents@xmin, extents@xmax)
  
  gbif_lat <- c(min(locations$lat), max(locations$lat))
  gbif_lon <- c(min(locations$lon), max(locations$lon))
  
  # create observation windows
  # shapefile window 
  spp_window <- owin(longitudes, latitudes)
  # gbif window
  gbif_window <- owin(gbif_lon, gbif_lat)
  
  # turn observation points into a point pattern process
  snake_ppp <- ppp(locations$lon, locations$lat, window = gbif_window)
  
  # define titles for the plots
  title_dens <- paste(genus, spp, 'density (default)', sep = " ")
  title_dens_sigma <- paste(genus, spp, 'density (sigma)', sep = " ")
  
  # generate density plots
  # default
  snake_dens <- density(snake_ppp)
  
  # standard deviations
  snake_dens_sigma <- density(snake_ppp, sigma = bw.diggle(snake_ppp))
  
  genus <- tolower(genus)
  
  # define output path
  path <- paste('Z:/users/joshua/Snakebite/density_plots/', genus, '_', spp, '_density.pdf', sep = '')
  
  # start a plotting window
  pdf(path,                    
      width = 8.27,
      height = 11.29)
  
  par(mfrow = c(2,2))
  
  plot(shape, main = title)
  
  plot(shape, main = title)
  points(locations$lon, locations$lat, pch = 20, cex = 0.5, col = "red")
  
  plot(snake_dens, main = title_dens)
  plot(shape, add = TRUE)
  
  plot(snake_dens_sigma, main = title_dens_sigma)
  plot(shape, add = TRUE)
  
  dev.off()  
  
}

