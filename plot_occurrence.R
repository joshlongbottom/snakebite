# script to grab all available gbif data per species, and plot this on top of WHO EOR maps
# clear workspace
rm(list = ls())

# load required packages
require(dismo)
require(raster)
require(maptools)
require(rgeos)

# load function
vector.is.empty <- function(x) return(length(x) == 0)

# read in snake species list
species_sheet <- read.csv('Z:/users/joshua/Snakebite/snakebite/snake_list.csv',
                          stringsAsFactors = FALSE)

# read in shapefiles
admin0 <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/admin2013_0.shp')
africa <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/Africa.shp')
america <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/America.shp')
latin_america <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/Latin_America.shp')
america_mex <- shapefile('Z:/users/joshua/Snakebite/World shapefiles/USA_MEX.shp')

# get a list of the unique snake species 
species_list <- unique(species_sheet$split_spp)
species_list <- as.list(species_list[!(species_list == "")])

# for every species in the list, grab occurrence records from gbif, and plot ontop
# of WHO EOR converted shapefile
# start a plotting window
pdf('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species_occurrence_plots_Y2017M04D18.pdf',                    
    width = 8.27,
    height = 11.29)
par(mfrow = c(3, 2))

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
    
    # inform progress
    message(paste('generating occurrence plot', i, 'of', length(species_list), sep = " "))
    
    # get unique country ISO codes
    iso <- unique(shape@data$COUNTRY_ID)
        # remove NAs and unassigned countries (XXX)
        iso <- iso[!is.na(iso)]
        iso <- iso[!(iso == 'XXX')]
        # get string of all countries for snake distribution
        countries <- paste(iso, collapse = ", ")
        countries_raw <- countries
        
    # create a shapefile of countries for each species
    country_shp <- admin0[admin0@data$COUNTRY_ID %in% iso, ]    
    
    if('USA' %in% iso){
      
      country_shp <- america
    }
    
    africa_list <- c('Bitis_arietans', 'Bitis_gabonica', 'Bitis_nasicornis', 'Bitis_rhinoceros', 'Atractaspis_irregularis', 
                     'Dendroaspis_polylepis', 'Dispholidus_typus', 'Naja_melanoleuca', 'Naja_nigricollis', 'Atractaspis_bibronii',
                     'Dendroaspis_angusticeps', 'Dendroaspis_jamesoni', 'Dendroaspis_polylepis', 'Dendroaspis_viridis', 'Echis_leucogaster',
                     'Echis_ocellatus', 'Echis_pyramidum', 'Hemachatus_haemachatus', 'Naja_anchietae', 'Naja_annulata', 
                     'Naja_annulifera', 'Naja_haje', 'Naja_katiensis', 'Naja_mossambica', 'Naja_nivea', 'Pseudohaje_goldii',
                     'Thelotornis_capensis', 'Thelotornis_kirklandii', 'Thelotornis_mossambicanus')
    
    if(species %in% africa_list){
      
      country_shp <- africa
    }
    
    latin_america_list <- c('Bothrocophias_hyoprora', 'Bothrocophias_microphthalmus', 'Bothrops_alternatus', 
                            'Bothrops_ammodytoides', 'Bothrops_attrox', 'Bothrops_bilineatus', 'Bothrops_brazili',
                            'Bothrops_diporus', 'Bothrops_jararaca', 'Bothrops_jararacussu', 'Bothrops_mattogrossensis',
                            'Bothrops_moojeni', 'Bothrops_neuwiedi', 'Bothrops_pubescens', 'Bothrops_spp','Bothrops_taeniatus',
                            'Crotalus_durissus', 'Lachesis_muta', 'Micrurus_corallinus', 'Micrurus_lemniscatus', 'Micrurus_spixii')
    
    if(species %in% latin_america_list){
      
      country_shp <- latin_america
    }
     
    usa_mex <- c('Crotalus_atrox', 'Crotalus_molossus', 'Crotalus_oreganus', 'Crotalus_ruber', 'Crotalus_scutulatus',
                 'Crotalus_spp', 'Crotalus_viridis', 'Micruroides_euryxanthus', 'Micrurus_fulvius', 'Micrurus_tener', 
                 'Sistrurus_catenatus')
    
    if(species %in% usa_mex){
      
      country_shp <- america_mex
    }
       
    if(nchar(countries) > 34){
      
        countries <- 'More than 7 countries'
    
        }    
    
    # if Russia is in ISO, sub out
    if('RUS' %in% iso){
      
      country_shp <- country_shp[!country_shp@data$COUNTRY_ID == 'RUS', ] 
    
      }
    
    # dissolve polygons
    n_polys <- as.numeric(nrow(shape@data))
    ids <- shape@data$ADMN_LEVEL  
    
    if(n_polys > 1){
      
      shape <- gBuffer(shape, byid = TRUE, width = 0)
      shape <- unionSpatialPolygons(shape, ids)
    
    }
    
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
    
    # define title for the plot
    title <- paste(genus, spp, sep = " ")
    species_n <- paste(genus, spp, sep = "_")
    
    # create a dataframe of each country within a species' range
    spec_dist <- as.data.frame(cbind(species_n, countries_raw))

    if(i == 1){
      
      spp_countries <- spec_dist
      
    } else {
      
      spp_countries <- rbind(spp_countries,
                             spec_dist)
    }
    
    # get gbif presence for that species
    # if the shapefile is for a group of species, don't grab results from gbif
    if(spp != 'spp'){
    
      spp_data <- gbif(genus, spp, end = 100000)
    
      } else {
      
      spp_data <- NA  
    
    }
    
    # if spp_data is not an empty dataframe get unique lat/longs for the species
    if((vector.is.empty(spp_data) == FALSE) & 'lat' %in% colnames(spp_data)) {
      
      locations <- unique(spp_data[c('lat', 'lon', 'year')])
    
      # if there are records, remove any NAs
      locations <- locations[!is.na(locations$lat), ]
      locations <- locations[!is.na(locations$lon), ]
    
      # write out the data points
      dat_path <- paste('Z:/users/joshua/Snakebite/output/species_occurrence_plots/species data/', species_n, '.csv', sep = "")
      write.csv(spp_data,
                dat_path)
      
      na_dates <- locations[is.na(locations$year), ]
      non_na <- locations[!is.na(locations$year), ]
      post_08 <- non_na[non_na$year >= 2008, ]
      pre_08 <- non_na[non_na$year < 2008, ]
      
      } else {
        
      locations <- NULL
          
      }

    # add number of data points to plot as a y axis
    if(vector.is.empty(locations) == FALSE) {
      
      row_length <- nrow(locations)
    
    } else {
      
      row_length <- 0
      
    }
    
    # records text
    records <- paste(row_length, 'records obtained from GBIF', sep = " ")
    
    # plot ranges, starting with background
    plot(country_shp, 
         col = 'grey',
         border = 'white',
         main = bquote(~italic(.(title))), 
         lwd = 0.3)  
    
    title(ylab = records, xlab = countries, line = 0)
    
    # plot shapefile
    plot(shape,
         add = TRUE,
         col = '#4775B5',
         border = '#686868',
         lty = 1,
         lwd = 0.3)
    plot(country_shp,
         add = TRUE,
         border = 'white',
         lwd = 0.3)
    
    # if there are records, plot them on top of shapefile
    if(vector.is.empty(locations) == FALSE){
      
      if(vector.is.empty(na_dates) == FALSE){
        points(na_dates$lon, na_dates$lat, pch = 20, cex = 0.2, col = '#FFFFBF')
        
      }
      
      if(vector.is.empty(pre_08) == FALSE){
      points(pre_08$lon, pre_08$lat, pch = 20, cex = 0.2, col = '#F59D6E')
    
      } 
      
      if(vector.is.empty(post_08) == FALSE){
        points(post_08$lon, post_08$lat, pch = 20, cex = 0.2, col = '#D93529')
        
      }
        }
  
    # add legend to plot
    legend('bottomleft', c("No date","Pre-2008", "Post-2008"), pch = 20,
          col = c("#FFFFBF","#F59D6E", "#D93529"), bty = 'n')

    # clear environment a little
    rm(spp_data,
       locations)
  }
  
dev.off()  
