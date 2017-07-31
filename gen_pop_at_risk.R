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

# load population density surface
pop_dens <- raster('Z:/users/joshua/Snakebite/rasters/population/Worldpop_GPWv4_Hybrid_201601_Global_Pop_5km_Adj_MGMatched_2015_Hybrid.tif')

# load admin 0 raster
admin_0 <- raster('Z:/users/joshua/Cx.tritaeniorhynchus/Culex Tritaeniorhynchus/data/clean/raster/polygon_rasters/admin0_raster.flt')

# read in admin 0 shapefile dbf
countries <- read.dbf("Z:/users/joshua/admin2013/admin2013_0.dbf",
                      as.is = TRUE)

# extend to the same extent as species richness surface
admin_0 <- extend(admin_0, species_richness, value = NA)

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

# correct West Bank and Gaza to PSE
global_pop$iso[global_pop$name == 'West Bank'] <- 'PSE'
global_pop$iso[global_pop$name == 'Gaza Strip'] <- 'PSE'

# aggregate based on iso code
aggregated_pop <- do.call(rbind,lapply(split(global_pop, global_pop$iso),function(df) sum(df$sum)))

new_df <- data.frame(iso = rep(NA, length(aggregated_pop)),
                     pop = rep(NA, length(aggregated_pop)))

new_df$iso <- row.names(aggregated_pop)
new_df$pop <- aggregated_pop

global_pop <- new_df

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
decile_pop <- do.call(rbind,lapply(split(global_pop, global_pop$decile),function(df) sum(df$pop)))

decile_population <- data.frame(decile = rep(NA, length(decile_pop)),
                                pop = rep(NA, length(decile_pop)))

decile_population$decile <- row.names(decile_pop)
decile_population$pop <- decile_pop

### generate population living in areas suitable for 1 or more med spp
# convert species richness surface into a binary surface (`1` = presence of 1 or more venomous
# snake species, `NoData` = absence)
species_presence <- species_richness
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
  ggtitle('Population at risk of exposure to one or more medically important snake species \nper HAQI decile')

ggsave('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_spp_hist.png', 
        dpi = 300, device = 'png')

# write out par of exposure to 1 or more snake species
# first write out the csv
write.csv(national_par,
          'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_spp.csv',
          row.names = FALSE)

# then the raster
writeRaster(presence_par, 
            file = 'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_spp',
            format = 'GTiff',
            overwrite = TRUE)

### generate population living in areas suitable for 1 or more medically important spp for which
### no effective therapy exists
# convert antivenom naive surface into a binary surface (`1` = presence of 1 or more therapy
# naive snake species, `NoData` = absence)
antiv_n_naive <- antivenom
antiv_n_naive[antiv_n_naive < 1 ] <- 0
antiv_n_naive[antiv_n_naive  >= 1] <- 1
  
# multiply the population by the binary presence/absence antivenom surface
naive_par <- overlay(pop_dens, antiv_n_naive, fun = function(pop_dens, antiv_n_naive){
                                                            (pop_dens*antiv_n_naive)})

# convert this to a national estimate using zonal()
national_naive_par <- zonal(naive_par, admin_0, fun = 'sum', na.rm = TRUE)
national_naive_par <- as.data.frame(national_naive_par)

# create an matching index
rm(match_idx)
match_idx <- match(national_naive_par$zone, countries$GAUL_CODE)

# append iso and country name
national_naive_par$iso <- countries$COUNTRY_ID[match_idx]
national_naive_par$name <- countries$name[match_idx]

# correct West Bank and Gaza to PSE
national_naive_par$iso[national_naive_par$name == 'West Bank'] <- 'PSE'
national_naive_par$iso[national_naive_par$name == 'Gaza Strip'] <- 'PSE'

# aggregate based on iso code
aggregated_par <- do.call(rbind,lapply(split(national_naive_par, national_naive_par$iso),function(df) sum(df$sum)))

new_df <- data.frame(iso = rep(NA, length(aggregated_par)),
                     par = rep(NA, length(aggregated_par)))

new_df$iso <- row.names(aggregated_par)
new_df$par <- aggregated_par

national_naive_par <- new_df

# merge HAQI with PAR
match_idx <- match(national_naive_par$iso, haqi$COUNTRY_ID)
national_naive_par$haqi <- haqi$haqi_2015[match_idx]

# add deciles
national_naive_par$decile[national_naive_par$haqi < 42.9] <- 1 
national_naive_par$decile[(national_naive_par$haqi >= 42.9) & (national_naive_par$haqi <= 47) ] <- 2
national_naive_par$decile[(national_naive_par$haqi > 47) & (national_naive_par$haqi <= 51.3) ] <- 3
national_naive_par$decile[(national_naive_par$haqi > 51.3) & (national_naive_par$haqi <= 59) ] <- 4
national_naive_par$decile[(national_naive_par$haqi > 59) & (national_naive_par$haqi <= 63.4) ] <- 5
national_naive_par$decile[(national_naive_par$haqi > 63.4) & (national_naive_par$haqi <= 69.7) ] <- 6
national_naive_par$decile[(national_naive_par$haqi > 69.7) & (national_naive_par$haqi <= 74.4) ] <- 7
national_naive_par$decile[(national_naive_par$haqi > 74.4) & (national_naive_par$haqi <= 79.4) ] <- 8
national_naive_par$decile[(national_naive_par$haqi > 79.4) & (national_naive_par$haqi <= 86.3) ] <- 9
national_naive_par$decile[national_naive_par$haqi > 86.3 ] <- 10

# get a total population at risk per decile
decile_naive_par <- do.call(rbind,lapply(split(national_naive_par, national_naive_par$decile),function(df) sum(df$par)))

naive_decile_par <- data.frame(decile = rep(NA, length(decile_naive_par)),
                               par = rep(NA, length(decile_naive_par)))

naive_decile_par$decile <- row.names(decile_naive_par)
naive_decile_par$par <- decile_naive_par

# gen % of decile at risk
combined <- cbind(naive_decile_par,
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
        ggtitle('Population at risk of exposure to one or more medically important snake species \nwith no effective therapy, per HAQI decile')

ggsave('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_therapy_naive_spp_hist.png', 
       dpi = 300, device = 'png')

# write out par of exposure to 1 or more therapy naive snake species
# first write out the csv
write.csv(national_naive_par,
          'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_therapy_naive_spp.csv',
          row.names = FALSE)

# then the raster
writeRaster(naive_par, 
            file = 'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_therapy_naive_spp',
            format = 'GTiff',
            overwrite = TRUE)

### no effective therapy exists (cat 1 only)
# convert antivenom naive surface into a binary surface (`1` = presence of 1 or more therapy
# naive snake species, `NoData` = absence)
antiv_c1_naive <- c1_antivenom
antiv_c1_naive[antiv_c1_naive < 1 ] <- 0
antiv_c1_naive[antiv_c1_naive  >= 1] <- 1

# multiply the population by the binary presence/absence antivenom surface
c1_naive_par <- overlay(pop_dens, antiv_c1_naive, fun = function(pop_dens, antiv_c1_naive){
                                                                (pop_dens*antiv_c1_naive)})

# convert this to a national estimate using zonal()
national_c1_naive_par <- zonal(c1_naive_par, admin_0, fun = 'sum', na.rm = TRUE)
national_c1_naive_par <- as.data.frame(national_c1_naive_par)

# create an matching index
rm(match_idx)
match_idx <- match(national_c1_naive_par$zone, countries$GAUL_CODE)

# append iso and country name
national_c1_naive_par$iso <- countries$COUNTRY_ID[match_idx]
national_c1_naive_par$name <- countries$name[match_idx]

# correct West Bank and Gaza to PSE
national_c1_naive_par$iso[national_c1_naive_par$name == 'West Bank'] <- 'PSE'
national_c1_naive_par$iso[national_c1_naive_par$name == 'Gaza Strip'] <- 'PSE'

# aggregate based on iso code
aggregated_par <- do.call(rbind,lapply(split(national_c1_naive_par, national_c1_naive_par$iso),function(df) sum(df$sum)))

new_df <- data.frame(iso = rep(NA, length(aggregated_par)),
                     par = rep(NA, length(aggregated_par)))

new_df$iso <- row.names(aggregated_par)
new_df$par <- aggregated_par

national_c1_naive_par <- new_df

# merge HAQI with PAR
match_idx <- match(national_c1_naive_par$iso, haqi$COUNTRY_ID)
national_c1_naive_par$haqi <- haqi$haqi_2015[match_idx]

# add deciles
national_c1_naive_par$decile[national_c1_naive_par$haqi < 42.9] <- 1 
national_c1_naive_par$decile[(national_c1_naive_par$haqi >= 42.9) & (national_c1_naive_par$haqi <= 47) ] <- 2
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 47) & (national_c1_naive_par$haqi <= 51.3) ] <- 3
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 51.3) & (national_c1_naive_par$haqi <= 59) ] <- 4
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 59) & (national_c1_naive_par$haqi <= 63.4) ] <- 5
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 63.4) & (national_c1_naive_par$haqi <= 69.7) ] <- 6
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 69.7) & (national_c1_naive_par$haqi <= 74.4) ] <- 7
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 74.4) & (national_c1_naive_par$haqi <= 79.4) ] <- 8
national_c1_naive_par$decile[(national_c1_naive_par$haqi > 79.4) & (national_c1_naive_par$haqi <= 86.3) ] <- 9
national_c1_naive_par$decile[national_c1_naive_par$haqi > 86.3 ] <- 10

# get a total population at risk per decile
decile_naive_par <- do.call(rbind,lapply(split(national_naive_par, national_naive_par$decile),function(df) sum(df$par)))

naive_decile_par <- data.frame(decile = rep(NA, length(decile_naive_par)),
                               par = rep(NA, length(decile_naive_par)))

naive_decile_par$decile <- row.names(decile_naive_par)
naive_decile_par$par <- decile_naive_par

# gen % of decile at risk
combined <- cbind(naive_decile_par,
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
  ggtitle('Population at risk of exposure to one or more medically important snake species \nwith no effective therapy, per HAQI decile')

ggsave('Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_therapy_naive_spp_hist.png', 
       dpi = 300, device = 'png')

# write out par of exposure to 1 or more therapy naive snake species
# first write out the csv
write.csv(national_c1_naive_par,
          'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c1_therapy_naive_spp.csv',
          row.names = FALSE)

# then the raster
writeRaster(c1_naive_par, 
            file = 'Z:/users/joshua/Snakebite/output/population_at_risk/exposure_to_one_or_more_c1_therapy_naive_spp',
            format = 'GTiff',
            overwrite = TRUE)

### generate 'snake-human exposure events' risk surface
exposure_events_par <- overlay(pop_dens, species_presence, fun = function(pop_dens, species_presence){
                                                                         (pop_dens*species_presence)})

exposure_events_par <- round(exposure_events_par)

# convert this to a national estimate using zonal()
national_exposure_par <- zonal(exposure_events_par, admin_0, fun = 'sum', na.rm = TRUE)
national_exposure_par <- as.data.frame(national_exposure_par)

# match this to get location names
# read in admin 0 shapefile dbf
countries <- foreign::read.dbf("Z:/users/joshua/admin2013/admin2013_0.dbf")

# create an matching index
match_idx <- match(national_exposure_par$zone, countries$GAUL_CODE)

# append iso and country name
national_exposure_par$iso <- countries$COUNTRY_ID[match_idx]
national_exposure_par$name <- countries$name[match_idx]

# write out par of exposure to 1 or more snake species
# first write out the csv
write.csv(national_exposure_par,
          'Z:/users/joshua/Snakebite/output/population_at_risk/snake_human_exposure_events.csv',
          row.names = FALSE)

# then the raster
writeRaster(exposure_events_par, 
            file = 'Z:/users/joshua/Snakebite/output/population_at_risk/snake_human_exposure_events',
            format = 'GTiff',
            overwrite = TRUE)



