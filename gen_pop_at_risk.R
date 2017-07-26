# script to generate population at risk estimates
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster)

# load species richness surface
species_richness <- raster('Z:/users/joshua/Snakebite/output/species_richness/modified_eor_combined_categories_2017-07-25.tif')

# load antivenom naive surface
antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom_combined_stack.tif')
c1_antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom_c1_stack.tif')
c2_antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom_c2_stack.tif')

# load population density surface
pop_dens <- raster('Z:/users/joshua/Snakebite/rasters/population/Worldpop_GPWv4_Hybrid_201601_Global_Pop_5km_Adj_MGMatched_2015_Hybrid.tif')

# load admin 0 raster
admin_0 <- raster('Z:/users/joshua/Cx.tritaeniorhynchus/Culex Tritaeniorhynchus/data/clean/raster/polygon_rasters/admin0_raster.flt')

# extend to the same extent as species richness surface
admin_0 <- extend(admin_0, species_richness, value = NA)

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
# read in admin 0 shapefile dbf
countries <- foreign::read.dbf("Z:/users/joshua/admin2013/admin2013_0.dbf")

# create an matching index
match_idx <- match(national_par$zone, countries$GAUL_CODE)

# append iso and country name
national_par$iso <- countries$COUNTRY_ID[match_idx]
national_par$name <- countries$name[match_idx]

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



