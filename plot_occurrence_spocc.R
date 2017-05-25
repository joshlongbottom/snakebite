# script to grab all available gbif data per species, and plot this on top of WHO EOR maps
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(dismo, raster, maptools, rgeos, plyr, spocc, scrubr, rvertnet, dplyr)

# set working directory
setwd('Z:/users/joshua/Snakebite/')

# load function(s)
source('snakebite/bespoke_functions.R')

# read in snake species list
species_sheet <- read.csv('snakebite/snake_list.csv',
                          stringsAsFactors = FALSE)

# limit species list to those for which I have digitised the range map
species_sheet <- species_sheet[!(species_sheet$shapefile_path == ""), ]

# read in continent shapefiles
admin0 <- shapefile('World shapefiles/admin2013_0.shp')
africa <- shapefile('World shapefiles/Africa.shp')
latin_america <- shapefile('World shapefiles/Latin_America.shp')
america_mex <- shapefile('World shapefiles/USA_MEX.shp')

# get a list of the unique snake species 
species_list <- unique(species_sheet$split_spp)
species_list <- as.list(species_list[!(species_list == "")])

# define lists for species which require a bespoke background plot
africa_list <- c('Bitis_arietans', 'Bitis_gabonica', 'Bitis_nasicornis', 'Atractaspis_irregularis', 
                 'Dendroaspis_polylepis', 'Dispholidus_typus', 'Naja_melanoleuca', 'Naja_nigricollis', 'Atractaspis_bibronii',
                 'Dendroaspis_polylepis', 'Naja_haje', 'Thelotornis_capensis')

latin_america_list <- c('Bothrocophias_hyoprora', 'Bothrocophias_microphthalmus', 'Bothrops_attrox', 'Bothrops_mattogrossensis',
                        'Bothrops_spp','Bothrops_taeniatus',
                        'Crotalus_durissus', 'Lachesis_muta', 'Micrurus_lemniscatus')

usa_mex <- c('Crotalus_atrox', 'Crotalus_molossus', 'Crotalus_oreganus', 'Crotalus_ruber', 'Crotalus_scutulatus',
             'Crotalus_spp', 'Crotalus_viridis', 'Micruroides_euryxanthus', 'Micrurus_fulvius', 'Micrurus_tener', 
             'Sistrurus_catenatus')

# for every species in the list, grab occurrence records from gbif, and plot ontop
# of WHO EOR converted shapefile

# start a plotting window
pdf_name <- paste('output/species_occurrence_plots/species_occurrence_plots_', Sys.Date(), '.pdf', sep = "")

pdf(pdf_name,                    
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
    
    # read in digitised range shapefile
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
    

    if(species %in% africa_list){
      
      country_shp <- africa
    }
    
    if(species %in% latin_america_list){
      
      country_shp <- latin_america
    }
    
    if(species %in% usa_mex){
      
      country_shp <- america_mex
    }
       
    if(nchar(countries) > 34){
      
        countries <- 'More than 7 countries'
    
        }    
    
    # dissolve polygons
    n_polys <- as.numeric(nrow(shape@data))
    ids <- shape@data$ADMN_LEVEL  
    
    if(n_polys > 1){
      
      suppressWarnings(shape <- gBuffer(shape, byid = TRUE, width = 0))
      suppressWarnings(shape <- unionSpatialPolygons(shape, ids))
    
    }
    
    # define species name    
    title <- sub$split_spp    
   
    # title <- gsub('_', ' ', title)
            
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
    
    if(!('spp' %in% title)){
      
      # get bison, ecoengine, idigbio records
      suppressWarnings(dat <- occ(title, from = c('gbif', 'bison', 'ecoengine', 
                                                  'idigbio', 'inat', 'vertnet'), 
                                  limit = 100000))

      # convert each into a separate dataframe
      dat_gbif <- as.data.frame(dat$gbif$data)
      dat_gbif$source <- rep('gbif', nrow(dat_gbif))
      names(dat_gbif) <- tolower(names(dat_gbif))
      
      dat_bison <- as.data.frame(dat$bison$data)
      dat_bison$source <- rep('bison', nrow(dat_bison))
      names(dat_bison) <- tolower(names(dat_bison))
      
      dat_ecoengine <- as.data.frame(dat$ecoengine$data)
      dat_ecoengine$source <- rep('ecoengine', nrow(dat_ecoengine))
      names(dat_ecoengine) <- tolower(names(dat_ecoengine))
      
      dat_idigbio <- as.data.frame(dat$idigbio$data)
      dat_idigbio$source <- rep('idigbio', nrow(dat_idigbio))
      names(dat_idigbio) <- tolower(names(dat_idigbio))
      
      dat_inat <- as.data.frame(dat$inat$data)
      dat_inat$source <- rep('inaturalist', nrow(dat_inat))
      names(dat_inat) <- tolower(names(dat_inat))
      
      dat_vertnet <- as.data.frame(dat$vertnet$data)
      dat_vertnet$source <- rep('vertnet', nrow(dat_vertnet))
      names(dat_vertnet) <- tolower(names(dat_vertnet))
      
      # merge into one dataframe
      dat_all <- rbind.fill(dat_gbif,
                           dat_bison,
                           dat_ecoengine,
                           dat_idigbio,
                           dat_inat,
                           dat_vertnet)
      
      # correct dataframe names
      sub_n <- gsub(' ', '_', title)
      sub_n <- tolower(sub_n)
      
      names(dat_all) <- gsub(sub_n, '', names(dat_all))
      names(dat_all) <- gsub('[.]', '', names(dat_all))
      
      # remove any listed elements in dat_all dataframe
      suppressWarnings(dat_all <- data.frame(lapply(dat_all, as.character), stringsAsFactors = FALSE))
  
      # convert to dataframe (only retains 'species name', 'longitude', 'latitude', 'prov', 'date', 'key')
      dat <- occ2df(dat)
      
    # if dat is not an empty dataframe get unique lat/longs for the species, across all repositories
    if(nrow(dat) > 0) {
      
      # clean coordinates by dropping incomplete
      suppressMessages(dat <- coord_incomplete(dat, lat = NULL, lon = NULL, drop = TRUE))
      
      if('date' %in% colnames(dat)){
      
        # split collection date string to get collection year
        dat$year <- substring(dat$date, 1, 4)
      
        # flag and drop duplicate lat/long/year combinations
        dat$duplicated <- duplicated(dat[c("longitude", "longitude", "year")])
        dat <- dat[dat$duplicated == FALSE, ]
      
        } else {
        
        # flag and drop duplicate lat/long combinations
        dat$duplicated <- duplicated(dat[c("longitude", "longitude")])
        dat <- dat[dat$duplicated == FALSE, ]
        
        }
      
      locations <- dat
      
      # write out the data points
      dat_path <- paste('output/species_occurrence_plots/species data/', species_n, '_raw.csv', sep = "")
      write.csv(dat_all,
                dat_path,
                na = "",
                row.names = FALSE)
      
      dat_path_trim <- paste('output/species_occurrence_plots/species data/', species_n, '_trimmed.csv', sep = "")
      write.csv(dat,
                dat_path_trim,
                row.names = FALSE)
      
      # split data by collection year
      if('year' %in% colnames(locations)){
      
        na_dates <- locations[is.na(locations$year), ]
        non_na <- locations[!is.na(locations$year), ]
        post_08 <- non_na[non_na$year >= 2008, ]
        pre_08 <- non_na[non_na$year < 2008, ]
      
        } else {
        
        na_dates <- locations
        non_na <- NULL
        post_08 <- NULL
        pre_08 <- NULL
        
      }
      
    } else {
      
      locations <- NULL
      
    }
    }  

    # add number of data points to plot as a y axis
    if(vector.is.empty(locations) == FALSE) {
      
      row_length <- nrow(locations)
    
    } else {
      
      row_length <- 0
      
    }
    
    # records text
    records <- paste(row_length, 'records obtained from online repositories', sep = " ")
    
    # plot ranges, starting with background
    plot(country_shp, 
         col = 'grey',
         border = 'white',
         main = bquote(~italic(.(title))), 
         lwd = 0.2)  
    
    title(ylab = records, xlab = countries, line = 0)
    
    # plot shapefile
    plot(shape,
         add = TRUE,
         col = '#4775B5',
         border = '#686868',
         lty = 1,
         lwd = 0.2)
    plot(country_shp,
         add = TRUE,
         border = 'white',
         lwd = 0.2)
    
    # if there are records, plot them on top of shapefile
    if(vector.is.empty(locations) == FALSE){
      
      if(vector.is.empty(na_dates) == FALSE){
        
        points(na_dates$longitude, na_dates$latitude, pch = 20, cex = 0.2, col = 'lightgoldenrod')
        
      }
      
      if(vector.is.empty(pre_08) == FALSE){
      
        points(pre_08$longitude, pre_08$latitude, pch = 20, cex = 0.2, col = '#F59D6E')
    
      } 
      
      if(vector.is.empty(post_08) == FALSE){
        
        points(post_08$longitude, post_08$latitude, pch = 20, cex = 0.2, col = '#D93529')
        
      }
        }
  
    # add legend to plot
    legend('bottomleft', c("No date","Pre-2008", "Post-2008"), pch = 20,
          col = c("lightgoldenrod","#F59D6E", "#D93529"), bty = 'n')

    # clear environment a little
    suppressWarnings(rm(spp_data, locations, na_dates, non_na, post_08, pre_08, spp, genus, dat_bison,
                        dat_gbif, dat_idigbio, dat_inat, dat_trim, dat_ecoengine, occ_trim, raw_occ, 
                        vert_trim, vert_date_c, dat, vert_latitude, vert_longitude, vert_record, vert_uncert,
                        vert_year, vert_year_idx, latitude, longitude, dat_2, date_c, split, col_year, 
                        col_year_idx, record, uncert, uncert_idx, dat_path, dat_path_trim, species_n,
                        sub_n, title))
    
  }
  
dev.off()  

# write out updated countries list
countries_path <- paste('output/species_occurrence_plots/species_country_list_', Sys.Date(), '.csv', sep = "")

write.csv(spp_countries,
          countries_path,
          row.names = FALSE)
