# script to generate population at risk estimates
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreign, reshape2, ggplot2, scales, RColorBrewer)

# load species richness surface
species_richness <- raster('Z:/users/joshua/Snakebite/output/species_richness/modified_eor_combined_categories_2017-07-25.tif')
c1_species_richness <- raster('Z:/users/joshua/Snakebite/output/species_richness/modified_eor_category_1_2017-07-25.tif')
c2_species_richness <- raster('Z:/users/joshua/Snakebite/output/species_richness/modified_eor_category_2_2017-07-25.tif')

# load antivenom naive surface
antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom_combined_stack.tif')
c1_antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom_c1_stack.tif')
c2_antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom_c2_stack.tif')

# define lists
par_list <- c('species_richness', 'c1_species_richness', 'c2_species_richness', 'antivenom', 'c1_antivenom', 'c2_antivenom')
title_vector <- c('Population at risk of exposure to one or more medically important snake species \nper HAQI decile',
                  'Population at risk of exposure to one or more Category 1 snake species \nper HAQI decile',
                  'Population at risk of exposure to one or more Category 2 snake species \nper HAQI decile',
                  'Population at risk of exposure to one or more medically important snake species \nwith no effective therapy, per HAQI decile',
                  'Population at risk of exposure to one or more Category 1 snake species \nwith no effective therapy, per HAQI decile',
                  'Population at risk of exposure to one or more Category 2 snake species \nwith no effective therapy, per HAQI decile')
outpath_vector <- c('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_spp',
                    'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c1_spp',
                    'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c2_spp',
                    'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_therapy_naive_spp',
                    'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c1_therapy_naive_spp',
                    'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c2_therapy_naive_spp')

# load population density surface
pop_dens <- raster('Z:/users/joshua/Snakebite/rasters/population/Worldpop_GPWv4_Hybrid_201601_Global_Pop_5km_Adj_MGMatched_2015_Hybrid.tif')

# load admin 0, and admin 1 raster
admin_0 <- raster('Z:/users/joshua/Snakebite/rasters/admin_0_updated_2017-08-01.tif')
admin_1 <- raster('Z:/users/joshua/Snakebite/rasters/admin_1.tif')

# load in accessibility surface
# accessibility <- raster('Z:/users/joshua/Snakebite/rasters/accessibility/accessibility_50k+_2017-01-05_final.tif')
accessibility <- raster('Z:/users/joshua/Snakebite/rasters/accessibility/accessibility_50k+_2017-01-05_aggregate_5k_2017_02_08.tif')
# # change -9999 to NA
# accessibility <- reclassify(accessibility, c(-10000, -1, NA))
# # resample to 5k
# accessibility <- resample(accessibility, species_richness, method = 'ngb')
# writeRaster(accessibility, 
#             file = 'Z:/users/joshua/Snakebite/rasters/accessibility/accessibility_50k+_2017-01-05_aggregate_5k_2017_02_08',
#             format = 'GTiff',
#             overwrite = TRUE)

# read in admin 0 shapefile dbf
countries <- read.dbf('Z:/users/joshua/Snakebite/World shapefiles/merged_admin0.dbf',
                      as.is = TRUE)

a1_dbf <- read.dbf('Z:/users/joshua/admin2013/admin2013_1.dbf',
                   as.is = TRUE)

# extend admin 0 and accessibility to the same extent as species richness surface
# get all extents
raster_list <- c(species_richness,
                 c1_species_richness,
                 c2_species_richness,
                 antivenom,
                 c1_antivenom,
                 c2_antivenom,
                 accessibility,
                 pop_dens,
                 admin_0,
                 admin_1)

# loop through and grab extents
extents <- t(sapply(raster_list, function (x) as.vector(extent(x))))

# get the smallest extent of all layers
ext <- extent(c(max(extents[, 1]),
                min(extents[, 2]),
                max(extents[, 3]),
                min(extents[, 4])))

# crop all layers by this minimal extent
species_richness <- crop(species_richness, ext)
c1_species_richness <- crop(c1_species_richness, ext)
c2_species_richness <- crop(c2_species_richness, ext)
antivenom <- crop(antivenom, ext)
c1_antivenom <- crop(c1_antivenom, ext)
c2_antivenom <- crop(c2_antivenom, ext)
accessibility <- crop(accessibility, ext)
pop_dens <- crop(pop_dens, ext)
admin_0 <- crop(admin_0, ext)
admin_1 <- crop(admin_1, ext)

# temp fix for weird admin_0 issue
# set extent equal to species richness
extent(admin_0) <- extent(species_richness)
extent(admin_1) <- extent(species_richness)

# read in HAQI data
haqi <- read.csv('Z:/users/joshua/Snakebite/HAQ_extract.csv',
                 stringsAsFactors = FALSE)

#### get global populations per country/HAQI, in order to calculate % of pop at risk ####
global_pop <- zonal(pop_dens, admin_0, fun = 'sum', na.rm = TRUE)
global_pop <- as.data.frame(global_pop)

# create an matching index
match_idx <- match(global_pop$zone, countries$GAUL_CODE)

# append iso and country name
global_pop$iso <- countries$COUNTRY_ID[match_idx]
global_pop$name <- countries$name[match_idx]

# merge HAQI with PAR
match_idx <- match(global_pop$iso, haqi$COUNTRY_ID)
global_pop$haqi <- haqi$haqi_2015[match_idx]

# add deciles
global_pop$decile[global_pop$haqi < 42.9] <- 1 
global_pop$decile[(global_pop$haqi >= 42.9) & (global_pop$haqi <= 47) ] <- 2
global_pop$decile[(global_pop$haqi > 47) & (global_pop$haqi <= 51.3) ] <- 3
global_pop$decile[(global_pop$haqi > 51.3) & (global_pop$haqi <= 59) ] <- 4
global_pop$decile[(global_pop$haqi > 59) & (global_pop$haqi <= 63.4) ] <- 5
global_pop$decile[(global_pop$haqi > 63.4) & (global_pop$haqi <= 69.7) ] <- 6
global_pop$decile[(global_pop$haqi > 69.7) & (global_pop$haqi <= 74.4) ] <- 7
global_pop$decile[(global_pop$haqi > 74.4) & (global_pop$haqi <= 79.4) ] <- 8
global_pop$decile[(global_pop$haqi > 79.4) & (global_pop$haqi <= 86.3) ] <- 9
global_pop$decile[global_pop$haqi > 86.3 ] <- 10

# get a total population per decile
decile_pop <- do.call(rbind,lapply(split(global_pop, global_pop$decile),function(df) sum(df$sum)))

decile_population <- data.frame(decile = rep(NA, length(decile_pop)),
                                pop = rep(NA, length(decile_pop)))

decile_population$decile <- row.names(decile_pop)
decile_population$pop <- decile_pop

##### loop through and generate population at risk estimates for: ####
# 1. Population living in areas suitable for one or more species (any medical classification)
# 2. Population living in areas suitable for one or more Category 1 species
# 3. Population living in areas suitable for one or more Category 2 species
# 4. Population living in areas suitable for one or more species for which no effective antivenom exists (any med class)
# 5. Population living in areas suitable for one or more Cat 1 species for which no effective antivenom exists
# 6. Population living in areas suitable for one or more Cat 2 species for which no effective antivenom exists

for(i in 1:length(par_list)){

# convert into a binary surface (`1` = presence of 1 or more species, `0` = absence)
species_presence <- get(par_list[[i]])
species_presence[species_presence < 1 ] <- 0
species_presence[species_presence  >= 1] <- 1
  
# multiply the population by the binary presence/absence surface
presence_par <- overlay(pop_dens, species_presence, fun = function(pop_dens, species_presence){
                                                                  (pop_dens*species_presence)})

# convert this to a national estimate using zonal()
national_par <- zonal(presence_par, admin_0, fun = 'sum', na.rm = TRUE)
national_par <- as.data.frame(national_par)

# match this to get location names
# create an matching index
match_idx <- match(national_par$zone, countries$GAUL_CODE)

# append iso and country name
national_par$iso <- countries$COUNTRY_ID[match_idx]
national_par$name <- countries$name[match_idx]

# correct West Bank and Gaza to PSE
national_par$iso[national_par$name == 'West Bank'] <- 'PSE'
national_par$iso[national_par$name == 'Gaza Strip'] <- 'PSE'

# aggregate based on iso code
aggregated_par <- do.call(rbind,lapply(split(national_par, national_par$iso),function(df) sum(df$sum)))

new_df <- data.frame(iso = rep(NA, length(aggregated_par)),
                     par = rep(NA, length(aggregated_par)))

new_df$iso <- row.names(aggregated_par)
new_df$par <- aggregated_par

national_par <- new_df

# merge HAQI with PAR
match_idx <- match(national_par$iso, haqi$COUNTRY_ID)
national_par$haqi <- haqi$haqi_2015[match_idx]

# add deciles
national_par$decile[national_par$haqi < 42.9] <- 1 
national_par$decile[(national_par$haqi >= 42.9) & (national_par$haqi <= 47) ] <- 2
national_par$decile[(national_par$haqi > 47) & (national_par$haqi <= 51.3) ] <- 3
national_par$decile[(national_par$haqi > 51.3) & (national_par$haqi <= 59) ] <- 4
national_par$decile[(national_par$haqi > 59) & (national_par$haqi <= 63.4) ] <- 5
national_par$decile[(national_par$haqi > 63.4) & (national_par$haqi <= 69.7) ] <- 6
national_par$decile[(national_par$haqi > 69.7) & (national_par$haqi <= 74.4) ] <- 7
national_par$decile[(national_par$haqi > 74.4) & (national_par$haqi <= 79.4) ] <- 8
national_par$decile[(national_par$haqi > 79.4) & (national_par$haqi <= 86.3) ] <- 9
national_par$decile[national_par$haqi > 86.3 ] <- 10

# get a total population at risk per decile
decile_par <- do.call(rbind,lapply(split(national_par, national_par$decile),function(df) sum(df$par)))

decile_nat_par <- data.frame(decile = rep(NA, length(decile_par)),
                             par = rep(NA, length(decile_par)))

decile_nat_par$decile <- row.names(decile_par)
decile_nat_par$par <- decile_par

# gen % of decile at risk
combined <- cbind(decile_nat_par,
                  decile_population)

combined[3] <- NULL
combined$percent_par <- (combined$par/combined$pop)*100
combined$remaining <- 100-combined$percent_par

# new dataframe
decile_plot_1 <- data.frame(decile = rep(NA, 10),
                            variable = rep(NA, 10),
                            value = rep(NA, 10))
decile_plot_2 <- data.frame(decile = rep(NA, 10),
                            variable = rep(NA, 10),
                            value = rep(NA, 10))

decile_plot_1$decile <- combined$decile
decile_plot_1$variable <- rep('Exposure to one or more species')
decile_plot_1$value <- combined$percent_par
decile_plot_2$decile <- combined$decile
decile_plot_2$variable <- rep('No risk of exposure')
decile_plot_2$value <- combined$remaining

combined <- rbind(decile_plot_1,
                  decile_plot_2)

combined$decile <- as.numeric(combined$decile)

# generate plot title
plot_title <- title_vector[[i]]
# generate plot outpath
plot_outpath <- paste0(outpath_vector[[i]], '_', Sys.Date(), '_hist.png')
csv_outpath <- paste0(outpath_vector[[i]], '_', Sys.Date(), '.csv')
geotiff_outpath <- paste0(outpath_vector[[i]], '_', Sys.Date())

# plot stacked barplot for population at risk
ggplot(combined, 
       aes(x = decile, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  scale_x_continuous(breaks = c(seq(0,10,1)))+
  scale_y_continuous(labels = percent) +
  labs(x = "Decile",
       y = "Population (%)")+
  theme(legend.title = element_blank(),
        # panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle(plot_title)

ggsave(plot_outpath, 
       dpi = 300, device = 'png')

# write out par of exposure to 1 or more snake species
# first write out the csv
write.csv(national_par,
          csv_outpath,
          row.names = FALSE)

# then the raster
writeRaster(presence_par, 
            file = geotiff_outpath,
            format = 'GTiff',
            overwrite = TRUE)

}

#### generate 'snake-human exposure events' risk surface ####
# this is the number of unique exposure events per person, per pixel. i.e. if there are 5 snakes in a pixel
# with 50 individuals, there are 250 possible snake-human exposure events (each person could be exposed to up to 
# 5 species)
exposure_events_par <- overlay(pop_dens, species_richness, fun = function(pop_dens, species_richness){
                                                                         (pop_dens*species_richness)})

# round up the values to integers
exposure_events_par <- round(exposure_events_par)

# convert this to a national estimate using zonal()
national_exposure_par <- zonal(exposure_events_par, admin_0, fun = 'sum', na.rm = TRUE)
national_exposure_par <- as.data.frame(national_exposure_par)

# match this to get location names
# create an matching index
match_idx <- match(national_exposure_par$zone, countries$GAUL_CODE)

# append iso and country name
national_exposure_par$iso <- countries$COUNTRY_ID[match_idx]
national_exposure_par$name <- countries$name[match_idx]

# merge HAQI with PAR
match_idx <- match(national_exposure_par$iso, haqi$COUNTRY_ID)
national_exposure_par$haqi <- haqi$haqi_2015[match_idx]

# merge with total population
match_idx <- match(national_exposure_par$iso, global_pop$iso)
national_exposure_par$population <- global_pop$sum[match_idx]

# add deciles
national_exposure_par$decile[national_exposure_par$haqi < 42.9] <- 1 
national_exposure_par$decile[(national_exposure_par$haqi >= 42.9) & (national_exposure_par$haqi <= 47) ] <- 2
national_exposure_par$decile[(national_exposure_par$haqi > 47) & (national_exposure_par$haqi <= 51.3) ] <- 3
national_exposure_par$decile[(national_exposure_par$haqi > 51.3) & (national_exposure_par$haqi <= 59) ] <- 4
national_exposure_par$decile[(national_exposure_par$haqi > 59) & (national_exposure_par$haqi <= 63.4) ] <- 5
national_exposure_par$decile[(national_exposure_par$haqi > 63.4) & (national_exposure_par$haqi <= 69.7) ] <- 6
national_exposure_par$decile[(national_exposure_par$haqi > 69.7) & (national_exposure_par$haqi <= 74.4) ] <- 7
national_exposure_par$decile[(national_exposure_par$haqi > 74.4) & (national_exposure_par$haqi <= 79.4) ] <- 8
national_exposure_par$decile[(national_exposure_par$haqi > 79.4) & (national_exposure_par$haqi <= 86.3) ] <- 9
national_exposure_par$decile[national_exposure_par$haqi > 86.3 ] <- 10

# generate average exposure per person, per country
national_exposure_par$average_exposure <- national_exposure_par$sum/national_exposure_par$population

# get a total population at risk per decile
decile_par <- do.call(rbind,lapply(split(national_exposure_par, national_exposure_par$decile),function(df) sum(df$sum)))

decile_nat_par <- data.frame(decile = rep(NA, length(decile_par)),
                             exposures = rep(NA, length(decile_par)))

decile_nat_par$decile <- row.names(decile_par)
decile_nat_par$exposures <- decile_par

# gen % of decile at risk
combined <- cbind(decile_nat_par,
                  decile_population)

combined[3] <- NULL
combined$average_exposure <- (combined$exposures/combined$pop)

# plot
plot(combined$decile, combined$average_exposure,
     xlab = 'Decile',
     ylab = 'Average snake-human exposures per person',
     xaxt = 'n')
axis(1, xaxp = c(0, 10, 10))

# write out par of exposure to 1 or more snake species
# first write out the csv
csv_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/snake_human_exposure_events', '_', Sys.Date(), '.csv')
write.csv(national_exposure_par,
          csv_outpath,
          row.names = FALSE)

# then the raster
geotiff_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/snake_human_exposure_events', '_', Sys.Date())
writeRaster(exposure_events_par, 
            file = geotiff_outpath,
            format = 'GTiff',
            overwrite = TRUE)

# generate admin 1, snake-human exposure events
# use 'exposure_events_par' from above; but admin 1 raster
# convert this to a national estimate using zonal()
admin_1_exposure_par <- zonal(exposure_events_par, admin_1, fun = 'sum', na.rm = TRUE)
admin_1_exposure_par <- as.data.frame(admin_1_exposure_par)

# match this to get location names
# create an matching index
match_idx <- match(admin_1_exposure_par$zone, a1_dbf$GAUL_CODE)

# append iso and country name
admin_1_exposure_par$iso <- a1_dbf$COUNTRY_ID[match_idx]
admin_1_exposure_par$name <- a1_dbf$name[match_idx]

# merge HAQI with PAR
match_idx <- match(admin_1_exposure_par$iso, haqi$COUNTRY_ID)
admin_1_exposure_par$haqi <- haqi$haqi_2015[match_idx]

# generate admin 1 population estimates
admin_1_pop <- zonal(pop_dens, admin_1, fun = 'sum', na.rm = TRUE)
admin_1_pop <- as.data.frame(admin_1_pop)

# create an matching index
match_idx <- match(admin_1_pop$zone, a1_dbf$GAUL_CODE)

# append iso and country name
admin_1_pop$iso <- a1_dbf$COUNTRY_ID[match_idx]
admin_1_pop$name <- a1_dbf$name[match_idx]

# merge HAQI with PAR
match_idx <- match(admin_1_pop$iso, haqi$COUNTRY_ID)
admin_1_pop$haqi <- haqi$haqi_2015[match_idx]

# merge with total population
match_idx <- match(admin_1_exposure_par$zone, admin_1_pop$zone)
admin_1_exposure_par$population <- admin_1_pop$sum[match_idx]

# generate average exposure per person, per admin 1
admin_1_exposure_par$average_exposure <- admin_1_exposure_par$sum/admin_1_exposure_par$population

# write out par of exposure to 1 or more snake species
# first write out the csv
csv_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/snake_human_exposure_events_admin_1_', Sys.Date(), '.csv')
write.csv(admin_1_exposure_par,
          csv_outpath,
          row.names = FALSE)

#### distance based mortality ####
# bin accessibility values to generate mortality likelihoods
# 1. set up matrix
vals <- matrix(ncol = 3,
               c(seq(0, 6000, 60),
                 seq(60, 6060, 60),
                 seq(0, 100, 1)), 
               byrow = FALSE)

# 2. reclassify into mortality bins (suppose this is similar to /60 and rounding vals, but I want
# 61-89 minutes to contribute towards 2% mortality likelihood, opposed to 1%).
accessibility_mortality <- reclassify(accessibility, vals)
accessibility_mortality <- reclassify(accessibility_mortality, c(101, 1000000000, 100))

# write out the new 'distance based mortality' raster
acc_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/distance_based_mortality_raw_percentage', '_', Sys.Date())
writeRaster(accessibility_mortality, 
            file = acc_outpath,
            format = 'GTiff',
            overwrite = TRUE)

# loop through and classify proportion of population within each time bin
for(i in 1:25){
  
  # message to inform progress
  message(paste0("Processsing distance '", i, "'"))
  
  # generate a binary time/distance surface
  if(i != 25){
    
    j <- i-1
    temp <- reclassify(accessibility_mortality, c(0, j, j))
    temp <- reclassify(accessibility_mortality, c(j, 101, NA))
    
  } else {
    
    temp <- accessibility_mortality
    
  }
  
  # mask population dens by this
  pop_mask <- pop_dens
  pop_mask <- mask(pop_mask, temp)
  
  # get zonal statistics
  national_pop_dist <- zonal(pop_mask, admin_0, fun = 'sum', na.rm = TRUE)
  national_pop_dist <- as.data.frame(national_pop_dist)
  
  # rename dataframe
  # create string for distance
  distance_n <- paste0('pop_within_', i, '_hours')
  names(national_pop_dist) <- c('zone',
                                distance_n)
  
  # bind dataframes
  if(i == 1){
    
    combined_frame <- national_pop_dist
    
  } else {
    
    combined_frame <- merge(combined_frame, national_pop_dist)
    
  }
  
}

# convert these raw distance-based populations into proportion of total population
# merge with global pop estimates
# first, rename dataframe
names(global_pop) <- c("zone",
                       "total_population",
                       "iso",
                       "name",
                       "haqi",
                       "decile")

# create matching index
match_idx1 <- match(combined_frame$zone, global_pop$zone)

# sub in total population values
combined_frame$total_population <- global_pop$total_population[match_idx1]

# generate % fields
combined_frame$pop_within_1_hours_percent <- (combined_frame$pop_within_1_hours/combined_frame$total_population)*100
combined_frame$pop_within_2_hours_percent <- (combined_frame$pop_within_2_hours/combined_frame$total_population)*100
combined_frame$pop_within_3_hours_percent <- (combined_frame$pop_within_3_hours/combined_frame$total_population)*100
combined_frame$pop_within_4_hours_percent <- (combined_frame$pop_within_4_hours/combined_frame$total_population)*100
combined_frame$pop_within_5_hours_percent <- (combined_frame$pop_within_5_hours/combined_frame$total_population)*100
combined_frame$pop_within_6_hours_percent <- (combined_frame$pop_within_6_hours/combined_frame$total_population)*100
combined_frame$pop_within_7_hours_percent <- (combined_frame$pop_within_7_hours/combined_frame$total_population)*100
combined_frame$pop_within_8_hours_percent <- (combined_frame$pop_within_8_hours/combined_frame$total_population)*100
combined_frame$pop_within_9_hours_percent <- (combined_frame$pop_within_9_hours/combined_frame$total_population)*100
combined_frame$pop_within_10_hours_percent <- (combined_frame$pop_within_10_hours/combined_frame$total_population)*100
combined_frame$pop_within_11_hours_percent <- (combined_frame$pop_within_11_hours/combined_frame$total_population)*100
combined_frame$pop_within_12_hours_percent <- (combined_frame$pop_within_12_hours/combined_frame$total_population)*100
combined_frame$pop_within_13_hours_percent <- (combined_frame$pop_within_13_hours/combined_frame$total_population)*100
combined_frame$pop_within_14_hours_percent <- (combined_frame$pop_within_14_hours/combined_frame$total_population)*100
combined_frame$pop_within_15_hours_percent <- (combined_frame$pop_within_15_hours/combined_frame$total_population)*100
combined_frame$pop_within_16_hours_percent <- (combined_frame$pop_within_16_hours/combined_frame$total_population)*100
combined_frame$pop_within_17_hours_percent <- (combined_frame$pop_within_17_hours/combined_frame$total_population)*100
combined_frame$pop_within_18_hours_percent <- (combined_frame$pop_within_18_hours/combined_frame$total_population)*100
combined_frame$pop_within_19_hours_percent <- (combined_frame$pop_within_19_hours/combined_frame$total_population)*100
combined_frame$pop_within_20_hours_percent <- (combined_frame$pop_within_20_hours/combined_frame$total_population)*100
combined_frame$pop_within_21_hours_percent <- (combined_frame$pop_within_21_hours/combined_frame$total_population)*100
combined_frame$pop_within_22_hours_percent <- (combined_frame$pop_within_22_hours/combined_frame$total_population)*100
combined_frame$pop_within_23_hours_percent <- (combined_frame$pop_within_23_hours/combined_frame$total_population)*100
combined_frame$pop_within_24_hours_percent <- (combined_frame$pop_within_24_hours/combined_frame$total_population)*100
combined_frame$pop_within_25_hours_percent <- (combined_frame$pop_within_25_hours/combined_frame$total_population)*100

# sub in ISO code & decile
combined_frame$iso <- global_pop$iso[match_idx1]
combined_frame$decile <- global_pop$decile[match_idx1]

# subset to have a raw population dataframe, and a % dataframe
raw_pop_distance <- combined_frame[c(53:54, 1:27)]
percent_pop_distance <- combined_frame[c(53:54, 1, 27:52)]

# write this dataframe to disk
combined_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/proportion_of_whole_population_within_x_distance_', Sys.Date(), '.csv')
write.csv(combined_frame,
          combined_outpath,
          row.names = FALSE)

# remove countries without a HAQi
percent_pop_distance <- percent_pop_distance[!is.na(percent_pop_distance$decile), ]

## for each decile, loop through and generate heatmaps
for(i in 1:10){
  
  # subset to get decile
  decile_frame <- percent_pop_distance[percent_pop_distance$decile == i, ]

  # drop some variables prior to reshaping
  decile_frame$zone <- NULL
  decile_frame$total_population <- NULL
  decile_frame$decile <- NULL
  
  melt_percentage_pop <- melt(decile_frame, id.vars = c('iso'))
  
  # rename variables in melted dataframe
  melt_percentage_pop$variable <- gsub('pop_within_1_hours_percent', '< 1', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_2_hours_percent', '< 2', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_3_hours_percent', '< 3', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_4_hours_percent', '< 4', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_5_hours_percent', '< 5', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_6_hours_percent', '< 6', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_7_hours_percent', '< 7', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_8_hours_percent', '< 8', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_9_hours_percent', '< 9', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_10_hours_percent', '< 10', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_11_hours_percent', '< 11', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_12_hours_percent', '< 12', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_13_hours_percent', '< 13', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_14_hours_percent', '< 14', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_15_hours_percent', '< 15', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_16_hours_percent', '< 16', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_17_hours_percent', '< 17', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_18_hours_percent', '< 18', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_19_hours_percent', '< 19', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_20_hours_percent', '< 20', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_21_hours_percent', '< 21', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_22_hours_percent', '< 22', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_23_hours_percent', '< 23', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_24_hours_percent', '< 24', melt_percentage_pop$variable, fixed = TRUE)
  melt_percentage_pop$variable <- gsub('pop_within_25_hours_percent', '24 or more', melt_percentage_pop$variable, fixed = TRUE)
  
  # round percentages to 2dp
  melt_percentage_pop$value <- as.numeric(melt_percentage_pop$value)
  melt_percentage_pop$value <- round(melt_percentage_pop$value, digits = 2)
  
  # plot data
  # define colours
  colours <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(20)
  
  # change variable to an ordered factor...
  melt_percentage_pop$variable <- factor(melt_percentage_pop$variable, c("< 1", "< 2", "< 3", "< 4",
                                                                         "< 5", "< 6", "< 7", "< 8",
                                                                         "< 9", "< 10", "< 11", "< 12",
                                                                         "< 13", "< 14", "< 15", "< 16",
                                                                         "< 17", "< 18", "< 19", "< 20",
                                                                         "< 21", "< 22", "< 23", "< 24", 
                                                                         "24 or more"))
  
  # define title text
  title_text <- paste0('HAQI Decile ', i)
  
  # create plot
  p <- ggplot(melt_percentage_pop, aes(x = iso, y = variable)) +
       geom_tile(aes(fill = cut(value, seq(0, 100, 5),
                             include.lowest=TRUE)), 
                 colour="white", 
                 size = 0.1) + 
       labs(x = 'Country',
            y = 'Hours from closest city with population \u2265 50,000') +
       ggtitle(title_text) +
       scale_fill_manual(values = colours,
                         labels = c("0-5",
                                    "5-10",
                                    "10-15",
                                    "15-20",
                                    "20-25",
                                    "25-30",
                                    "30-35",
                                    "35-40",
                                    "40-45",
                                    "45-50",
                                    "50-55",
                                    "55-60",
                                    "60-65",
                                    "65-70",
                                    "70-75",
                                    "75-80",
                                    "80-85",
                                    "85-90",
                                    "90-95",
                                    "95-100"),
                        name = "Proportion of population (%)",
                        drop = FALSE)
  
  # rotate x axis labels and remove axis ticks
  p + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.ticks = element_blank())
  
  # save output
  # generate outpath
  distance_opath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/decile_',i, '_proportion_of_whole_population_per_distance_', Sys.Date(), '.png')
  ggsave(distance_opath, width = 400, height = 350, units = 'mm', dpi = 300, device = 'png')
  
}

#### now generate these estimates for the proportion of the PAR ####
# load PAR surface from above
# any exposure
exposure_any <- raster('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_spp_2017-08-01.tif')
exposure_any <- crop(exposure_any, admin_0)

# any c1 exposure
exposure_c1 <- raster('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c1_spp_2017-08-01.tif')
exposure_c1 <- crop(exposure_c1, admin_0)

# any c1 exposure
exposure_c2 <- raster('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c2_spp_2017-08-01.tif')
exposure_c2 <- crop(exposure_c2, admin_0)

# list metrics
process_list <- c('exposure_any', 'exposure_c1', 'exposure_c2')

# cycle through and generate stats
for(i in 1:length(process_list)){
  
  # get appropriate raster
  process_raster <- get(process_list[[i]])
  
  # inform progress
  message(paste0('Processing raster ', i, ' of ', length(process_list)))
  
  # generate national estimates using zonal
  exposure_pop <- zonal(process_raster, admin_0, fun = 'sum', na.rm = TRUE)
  exposure_pop <- as.data.frame(exposure_pop)
  
  # create an matching index
  match_idx <- match(exposure_pop$zone, countries$GAUL_CODE)
  
  # append iso and country name
  exposure_pop$iso <- countries$COUNTRY_ID[match_idx]
  exposure_pop$name <- countries$name[match_idx]
  
  # merge HAQI with PAR
  match_idx <- match(exposure_pop$iso, haqi$COUNTRY_ID)
  exposure_pop$haqi <- haqi$haqi_2015[match_idx]
  
  # add deciles
  exposure_pop$decile[exposure_pop$haqi < 42.9] <- 1 
  exposure_pop$decile[(exposure_pop$haqi >= 42.9) & (exposure_pop$haqi <= 47) ] <- 2
  exposure_pop$decile[(exposure_pop$haqi > 47) & (exposure_pop$haqi <= 51.3) ] <- 3
  exposure_pop$decile[(exposure_pop$haqi > 51.3) & (exposure_pop$haqi <= 59) ] <- 4
  exposure_pop$decile[(exposure_pop$haqi > 59) & (exposure_pop$haqi <= 63.4) ] <- 5
  exposure_pop$decile[(exposure_pop$haqi > 63.4) & (exposure_pop$haqi <= 69.7) ] <- 6
  exposure_pop$decile[(exposure_pop$haqi > 69.7) & (exposure_pop$haqi <= 74.4) ] <- 7
  exposure_pop$decile[(exposure_pop$haqi > 74.4) & (exposure_pop$haqi <= 79.4) ] <- 8
  exposure_pop$decile[(exposure_pop$haqi > 79.4) & (exposure_pop$haqi <= 86.3) ] <- 9
  exposure_pop$decile[exposure_pop$haqi > 86.3 ] <- 10
  
  # get a total population per decile
  exposure_decile_pop <- do.call(rbind,lapply(split(exposure_pop, exposure_pop$decile),function(df) sum(df$sum)))
  
  exposure_decile_population <- data.frame(decile = rep(NA, length(exposure_decile_pop)),
                                           pop = rep(NA, length(exposure_decile_pop)))
  
  exposure_decile_population$decile <- row.names(exposure_decile_pop)
  exposure_decile_population$pop <- exposure_decile_pop
  
  # loop through and classify proportion of population within each time bin
  for(n in 1:25){
    
    # message to inform progress
    message(paste0("Processsing distance '", n, "'"))
    
    # generate a binary time/distance surface
    if(n != 25){
      
      j <- n-1
      temp <- reclassify(accessibility_mortality, c(0, j, j))
      temp <- reclassify(accessibility_mortality, c(j, 101, NA))
      
    } else {
      
      temp <- accessibility_mortality
      
    }
    
    # mask exposure surface by this
    pop_mask <- process_raster
    pop_mask <- mask(pop_mask, temp)
    
    # get zonal statistics
    national_pop_dist <- zonal(pop_mask, admin_0, fun = 'sum', na.rm = TRUE)
    national_pop_dist <- as.data.frame(national_pop_dist)
    
    # rename dataframe
    # create string for distance
    distance_n <- paste0('pop_within_', n, '_hours')
    names(national_pop_dist) <- c('zone',
                                  distance_n)
    
    # bind dataframes
    if(n == 1){
      
      combined_frame <- national_pop_dist
      
    } else {
      
      combined_frame <- merge(combined_frame, national_pop_dist)
      
    }
    
  }
  
  # convert these raw distance-based populations into proportion of total population
  # merge with global pop estimates
  # first, rename dataframe
  names(exposure_pop) <- c("zone",
                           "total_population",
                           "iso",
                           "name",
                           "haqi",
                           "decile")
  
  # create matching index
  match_idx1 <- match(combined_frame$zone, exposure_pop$zone)
  
  # sub in total population values, ISO code & decile
  combined_frame$total_population <- exposure_pop$total_population[match_idx1]
  combined_frame$iso <- exposure_pop$iso[match_idx1]
  combined_frame$decile <- exposure_pop$decile[match_idx1]
  
  # generate % fields
  combined_frame$pop_within_1_hours_percent <- (combined_frame$pop_within_1_hours/combined_frame$total_population)*100
  combined_frame$pop_within_2_hours_percent <- (combined_frame$pop_within_2_hours/combined_frame$total_population)*100
  combined_frame$pop_within_3_hours_percent <- (combined_frame$pop_within_3_hours/combined_frame$total_population)*100
  combined_frame$pop_within_4_hours_percent <- (combined_frame$pop_within_4_hours/combined_frame$total_population)*100
  combined_frame$pop_within_5_hours_percent <- (combined_frame$pop_within_5_hours/combined_frame$total_population)*100
  combined_frame$pop_within_6_hours_percent <- (combined_frame$pop_within_6_hours/combined_frame$total_population)*100
  combined_frame$pop_within_7_hours_percent <- (combined_frame$pop_within_7_hours/combined_frame$total_population)*100
  combined_frame$pop_within_8_hours_percent <- (combined_frame$pop_within_8_hours/combined_frame$total_population)*100
  combined_frame$pop_within_9_hours_percent <- (combined_frame$pop_within_9_hours/combined_frame$total_population)*100
  combined_frame$pop_within_10_hours_percent <- (combined_frame$pop_within_10_hours/combined_frame$total_population)*100
  combined_frame$pop_within_11_hours_percent <- (combined_frame$pop_within_11_hours/combined_frame$total_population)*100
  combined_frame$pop_within_12_hours_percent <- (combined_frame$pop_within_12_hours/combined_frame$total_population)*100
  combined_frame$pop_within_13_hours_percent <- (combined_frame$pop_within_13_hours/combined_frame$total_population)*100
  combined_frame$pop_within_14_hours_percent <- (combined_frame$pop_within_14_hours/combined_frame$total_population)*100
  combined_frame$pop_within_15_hours_percent <- (combined_frame$pop_within_15_hours/combined_frame$total_population)*100
  combined_frame$pop_within_16_hours_percent <- (combined_frame$pop_within_16_hours/combined_frame$total_population)*100
  combined_frame$pop_within_17_hours_percent <- (combined_frame$pop_within_17_hours/combined_frame$total_population)*100
  combined_frame$pop_within_18_hours_percent <- (combined_frame$pop_within_18_hours/combined_frame$total_population)*100
  combined_frame$pop_within_19_hours_percent <- (combined_frame$pop_within_19_hours/combined_frame$total_population)*100
  combined_frame$pop_within_20_hours_percent <- (combined_frame$pop_within_20_hours/combined_frame$total_population)*100
  combined_frame$pop_within_21_hours_percent <- (combined_frame$pop_within_21_hours/combined_frame$total_population)*100
  combined_frame$pop_within_22_hours_percent <- (combined_frame$pop_within_22_hours/combined_frame$total_population)*100
  combined_frame$pop_within_23_hours_percent <- (combined_frame$pop_within_23_hours/combined_frame$total_population)*100
  combined_frame$pop_within_24_hours_percent <- (combined_frame$pop_within_24_hours/combined_frame$total_population)*100
  combined_frame$pop_within_25_hours_percent <- (combined_frame$pop_within_25_hours/combined_frame$total_population)*100
  
  # drop countries with 0 PAR of exposure
  combined_frame <- combined_frame[!(combined_frame$total_population == 0), ]
  
  # subset to have a raw population dataframe, and a % dataframe
  raw_pop_distance <- combined_frame[c(28:29, 1:27)]
  percent_pop_distance <- combined_frame[c(28:29, 1, 27, 30:54)]
  
  # loop vector
  loop_vector <- process_list[i]
  
  # write this dataframe to disk
  combined_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/proportion_of_', loop_vector, '_PAR_within_x_distance_', Sys.Date(), '.csv')
  write.csv(combined_frame,
            combined_outpath,
            row.names = FALSE)
  
  # remove countries without a HAQi
  percent_pop_distance <- percent_pop_distance[!is.na(percent_pop_distance$decile), ]
  
  ## for each decile, loop through and generate heatmaps
  for(d in 1:10){
    
    # subset to get decile
    decile_frame <- percent_pop_distance[percent_pop_distance$decile == d, ]
    
    # drop some variables prior to reshaping
    decile_frame$zone <- NULL
    decile_frame$total_population <- NULL
    decile_frame$decile <- NULL
    
    melt_percentage_pop <- melt(decile_frame, id.vars = c('iso'))
    
    # rename variables in melted dataframe
    melt_percentage_pop$variable <- gsub('pop_within_1_hours_percent', '< 1', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_2_hours_percent', '< 2', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_3_hours_percent', '< 3', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_4_hours_percent', '< 4', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_5_hours_percent', '< 5', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_6_hours_percent', '< 6', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_7_hours_percent', '< 7', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_8_hours_percent', '< 8', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_9_hours_percent', '< 9', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_10_hours_percent', '< 10', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_11_hours_percent', '< 11', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_12_hours_percent', '< 12', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_13_hours_percent', '< 13', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_14_hours_percent', '< 14', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_15_hours_percent', '< 15', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_16_hours_percent', '< 16', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_17_hours_percent', '< 17', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_18_hours_percent', '< 18', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_19_hours_percent', '< 19', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_20_hours_percent', '< 20', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_21_hours_percent', '< 21', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_22_hours_percent', '< 22', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_23_hours_percent', '< 23', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_24_hours_percent', '< 24', melt_percentage_pop$variable, fixed = TRUE)
    melt_percentage_pop$variable <- gsub('pop_within_25_hours_percent', '24 or more', melt_percentage_pop$variable, fixed = TRUE)
    
    # round percentages to 2dp
    melt_percentage_pop$value <- as.numeric(melt_percentage_pop$value)
    melt_percentage_pop$value <- round(melt_percentage_pop$value, digits = 2)
    
    # plot data
    # define colours
    colours <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(20)
    
    # change variable to an ordered factor...
    melt_percentage_pop$variable <- factor(melt_percentage_pop$variable, c("< 1", "< 2", "< 3", "< 4",
                                                                           "< 5", "< 6", "< 7", "< 8",
                                                                           "< 9", "< 10", "< 11", "< 12",
                                                                           "< 13", "< 14", "< 15", "< 16",
                                                                           "< 17", "< 18", "< 19", "< 20",
                                                                           "< 21", "< 22", "< 23", "< 24", 
                                                                           "24 or more"))
    
    # define title text
    if(loop_vector == "exposure_c2"){
    
      title_text <- paste0('HAQI Decile ', d, ' - Exposure to one or more category 2 species')
    
    } else {
      
      if(loop_vector == "exposure_c1"){
        
        title_text <- paste0('HAQI Decile ', d, ' - Exposure to one or more category 1 species')
        
      } else {
        
        if(loop_vector == "exposure_any"){
          
          title_text <- paste0('HAQI Decile ', d, ' - Exposure to one or more medically important species')
          
        }
      }
    }
    
    # create plot
    p <- ggplot(melt_percentage_pop, aes(x = iso, y = variable)) +
      geom_tile(aes(fill = cut(value, seq(0, 100, 5),
                               include.lowest=TRUE)), 
                colour="white", 
                size = 0.1) + 
      labs(x = 'Country',
           y = 'Hours from closest city with population \u2265 50,000') +
      ggtitle(title_text) +
      scale_fill_manual(values = colours,
                        labels = c("0-5",
                                   "5-10",
                                   "10-15",
                                   "15-20",
                                   "20-25",
                                   "25-30",
                                   "30-35",
                                   "35-40",
                                   "40-45",
                                   "45-50",
                                   "50-55",
                                   "55-60",
                                   "60-65",
                                   "65-70",
                                   "70-75",
                                   "75-80",
                                   "80-85",
                                   "85-90",
                                   "90-95",
                                   "95-100"),
                        name = "Proportion of population (%)",
                        drop = FALSE)
    
    # rotate x axis labels and remove axis ticks
    p + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
              axis.ticks = element_blank())
    
    # save output
    # generate outpath
    distance_opath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/decile_',d, '_proportion_of_', loop_vector, '_PAR_per_distance_', Sys.Date(), '.png')
    ggsave(distance_opath, width = 400, height = 350, units = 'mm', dpi = 300, device = 'png')
    
  }
  
  rm(exposure_pop,
     combined_frame)
  
}

