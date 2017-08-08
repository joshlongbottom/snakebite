# clear workspace
rm(list = ls())

# load required functions
pacman::p_load(raster, maptools)

# read in functions file
source('code/bespoke_functions_cluster.R')

# read in admin 0
admin_0 <- shapefile('data/raw/world_shapefiles/admin2013_0.shp')

# read in the global extent raster
global <- raster('data/raw/raster/land_sea/CoastGlobal_5k.tif')

# specify outpath for EOR pngs
png_eor <- 'output/eor_png/'

# list species which have a MESS (>5 records)
png_plots <- list.files('output/mess_png/',
                        pattern = '*png$',
                        full.names = FALSE)

# remove dates
date_v <- '2017-07-19'
png_plots <- gsub(date_v, '', png_plots)
png_plots <- gsub('_species_mess_maps_.png', '', png_plots, fixed = TRUE)

# list all species
snake_list <- read.csv('data/raw/snake_list_cluster.csv',
                       stringsAsFactors = FALSE,
                       na.strings = c("NA","NaN", " ", ""))

sub_list <- snake_list[!is.na(snake_list$shapefile_path),]

individual_list <- unique(sub_list$split_spp)

# now, list species which have not already been plot, but for which we have a range
to_plot <- individual_list[!(individual_list %in% png_plots)]

# loop through these species, and generate a plot of the un-altered EOR
for(i in 1:length(to_plot)){
  
  # get species name
  spp_name <- as.character(to_plot[i])
  
  # get sub list
  sub <- snake_list[snake_list$split_spp == spp_name, ]
  
  # create extent shapefile
  ext <- generate_ext(sub, admin_0)
  
  # crop the global extent by the country extent
  glob_ext <- crop(global, ext)
  
  # read in EOR shapefile for each species, and dissolve if >1 polygon
  range <- prepare_eor(sub)
  
  # create a plotting window
  png_name <- paste(png_eor, spp_name, '_species_eor_maps_', Sys.Date(), '.png', sep = "")
  
  species_sub <- gsub('_', ' ', spp_name)
  
  png(png_name,
      width = 450,
      height = 400,
      units = 'mm',
      res = 300)
  par(oma = c(2, 2, 2, 2))
  
  breakpoints <- c(0, 0.5, 1)
  colours <- c('white', '#f2f2f2', '#f2f2f2')
  
  plot(glob_ext,
       breaks = breakpoints,
       col = colours,
       legend = FALSE,
       axes = FALSE,
       box = FALSE)

  plot(range,
       add = TRUE,
       col = '#c1de29',
       border = NA)
  
  plot(ext,
       add = TRUE,
       border = 'gray45',
       lty = 1,
       lwd = 0.2)
  
  legend('bottomleft', c("Current EOR"), 
         pch = c(15),
         col = c("#c1de29"), bty = 'n', pt.cex = 2)
  
  mtext(bquote(~italic(.(species_sub))), side = 3, line = -1, outer = TRUE, cex = 2, font = 2)
  
  dev.off()
  
}
  
  