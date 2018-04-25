# compare hospitals
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreign, reshape2, ggplot2, scales, RColorBrewer, seegSDM, rowr)

# set working directory
setwd('H:/Work/Snakebite')

# load in travel time surface
travel_surface <- raster('rasters/accessibility/accessibility_50k+_2017-07-31_aggregate_5k_2017_08_09.tif')

# load in hospital network
hospitals <- read.csv('africa hospitals/Ouma_Okiro_Snow_Africa_Hospitals_Data.tab.csv', 
                      stringsAsFactors = FALSE)

# get coordinates
hospital_locs <- hospitals[c('Long', 'Lat')]

# get travel times at each of these points
hospital_vals <- extract(travel_surface, hospital_locs)

# convert to hours
hospital_vals <- hospital_vals/60

# add back into dataset
hospitals$distance <- hospital_vals

breaks <- seq(0, 143, 0.5)

hospital_cut <- cut(hospital_vals, breaks, right=FALSE)

hospital_freq <- table(hospital_cut)

# randomly sample some points in Africa
# read in Africa shapefile
africa <- shapefile('World shapefiles/Africa.shp')
africa_raster <- rasterize(africa, travel_surface)

# generate samples (100 iterations)
for(i in 1:100){

  africa_samples <- bgSample(africa_raster, n = 4908)

  # generate values
  random_distances <- extract(travel_surface, africa_samples)
  # convert to hours
  random_distance_vals <- random_distances/60
  
  if(i == 1){
    
    random_samples <- random_distance_vals
  
  } else {
      
    random_samples <- cbind.fill(random_samples, 
                                 random_distance_vals, fill = 0)
    }
  
}

random_samples <- as.data.frame(random_samples)

for(i in 1:100){
  
  breaks <- seq(0, 143, 0.5)
  
  random_sample <- random_samples[c(i)]
  random_sample <- unlist(random_sample)
  random_cut <- cut(random_sample, breaks, right=FALSE)
  
  random_freq <- table(random_cut)

  cumfreq0 <- c(0, cumsum(random_freq))
  
  if(i == 1){
  
    plot(breaks, cumfreq0,
       main ="Cumulative Frequency",
       xlab ="Duration (Hours)",
       ylab ="Cumulative points", type = 'l')

  } else {
    
    lines(breaks, cumfreq0, col = 'black')

  }
}

abline(a = 4908, b = 0, lty = 2, col = 'red')
cumfreq0 <- c(0, cumsum(hospital_freq))
lines(breaks, cumfreq0, col = 'blue')

# loop through and plot based on country
hosp_duplicate <- hospitals
countries <- sort(unique(hosp_duplicate$Country))

pdf('africa hospitals/comparison_plots.pdf',                    
    width = 8.27,
    height = 11.29)
par(mfrow = c(4, 3))

br <- seq(0, 143, by = 0.5)

for (s in seq_along(countries)) { 
  
  sub <- hosp_duplicate[which(hosp_duplicate$Country == countries[s]), ]
  
  # set title for the plot
  title <- unique(sub$Country) 
  
  max_m <- round(max((sub$distance), na.rm = TRUE))
  
  hist(sub$distance,
       col = grey(0.4),
       xlim = c(0, max_m),
       border = 'white',
       breaks = br,
       xlab = 'Duration (Hours)',
       main = title)
} 

dev.off()
