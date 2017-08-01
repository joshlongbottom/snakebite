# script to generate population at risk estimates
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreign, reshape2, ggplot2, scales)

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

# load admin 0 raster
admin_0 <- raster('Z:/users/joshua/Snakebite/rasters/admin_0_updated_2017-08-01.tif')

# load in accessibility surface
accessibility <- raster('Z:/users/joshua/Snakebite/rasters/accessibility/accessibility_50k+_2017-01-05_final.tif')

# read in admin 0 shapefile dbf
countries <- read.dbf('Z:/users/joshua/Snakebite/World shapefiles/merged_admin0.dbf',
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
                 admin_0)

# loop through and grab extents
extents <- t(sapply(raster_list, function (x) as.vector(extent(x))))

# get the smallest extent squares of all layers
ext <- extent(c(max(extents[, 1]),
                min(extents[, 2]),
                max(extents[, 3]),
                min(extents[, 4])))

# crop all layers by this
species_richness <- crop(species_richness, ext)
c1_species_richness <- crop(c1_species_richness, ext)
c2_species_richness <- crop(c2_species_richness, ext)
antivenom <- crop(antivenom, ext)
c1_antivenom <- crop(c1_antivenom, ext)
c2_antivenom <- crop(c2_antivenom, ext)
accessibility <- crop(accessibility, ext)
pop_dens <- crop(pop_dens, ext)
admin_0 <- crop(admin_0, ext)

# aggregate accessibility to 5km resolution
accessibility <- aggregate(accessibility, fact = 5, fun = mean)

# read in HAQI data
haqi <- read.csv('Z:/users/joshua/Snakebite/HAQ_extract.csv',
                 stringsAsFactors = FALSE)

#### get global populations per country/HAQI, in order to calculate % of pop at risk
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

# loop through and generate population at risk estimates for
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

### generate 'snake-human exposure events' risk surface
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

### distance based mortality
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
accessibility_mortality <- reclassify(accessibility_mortality, c(101, 1000000000, 100, -10000, -1, NA))

# write out the new 'distance based mortality' raster
acc_outpath <- paste0('Z:/users/joshua/Snakebite/output/population_at_risk/distance_based_mortality_raw_percentage', '_', Sys.Date())
writeRaster(accessibility_mortality, 
            file = acc_outpath,
            format = 'GTiff',
            overwrite = TRUE)



