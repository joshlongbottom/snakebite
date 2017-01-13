# script to grab all available gbif data per species, and plot this on top of WHO EOR maps
# clear workspace
rm(list = ls())

# load required packages
require(dismo)
require(raster)

# load function
vector.is.empty <- function(x) return(length(x) == 0)

# read in snake species list
species_sheet <- read.csv('Z:/users/joshua/Snakebite/snakebite/snake_list.csv',
                          stringsAsFactors = FALSE)

# get a list of the unique snake species 
species_list <- unique(species_sheet$split_spp)
species_list <- as.list(species_list[!(species_list == "")])

# for every species in the list, grab occurrence records from gbif, and plot ontop
# of WHO EOR converted shapefile

# start a plotting window
pdf('Z:/users/joshua/Snakebite/species_occurrence_plots.pdf',                    
    width = 8.27,
    height = 11.29)
par(mfrow = c(4, 3))

# loop through and populate pdf with species occurrence plots
for(i in 1:length(species_list)){
  
    # specify species
    species <- species_list[i]
    
    # subset species sheet (to obtain file path to shapefile)
    sub <- species_sheet[species_sheet$split_spp == species, ]
    
    # define shapefile for that species
    shape_path <- sub$shapefile_path
    
    # read in shapefile
    shape <- shapefile(shape_path)
    
    # define title for the plot
    title <- as.character(species)
    
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
    
    # if there are records, plot them on top of shapefile
    if(vector.is.empty(locations) == FALSE){
      
    # plot shapefile
    plot(shape, main = title)
    
    # plot coordinates
    points(locations$lon, locations$lat, pch = 20, cex = 0.5, col = "red")
    
    } else{
      
    # just plot the shapefile
    plot(shape, main = title)  
    
      }
  }
  
dev.off()  
