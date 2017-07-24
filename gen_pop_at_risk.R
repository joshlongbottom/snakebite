# script to generate population at risk estimates
# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster, foreach, doMC)

# load species richness surface
species_richness <- raster('Z:/users/joshua/Snakebite/output/species_richness.tif')

# load antivenom naive surface
antivenom <- raster('Z:/users/joshua/Snakebite/output/antivenom_coverage/No_specific_antivenom.tif')

# load population density surface
pop_dens <- raster('Z:/mastergrids/Other_Global_Covariates/Population/Worldpop_GPWv4_Hybrid_201601/Global_Pop_1km_Adj_MGMatched_2015_Hybrid.tif')

# load admin 0 raster
admin_0 <- brick('Z:/users/joshua/Cx.tritaeniorhynchus/Culex Tritaeniorhynchus/data/clean/raster/polygon_rasters/admin0_raster.flt')

# rasterize spp richness raster to be at the same spatial res as the 1km pop
species_richness <- disaggregate(species_richness, fact = c(5,5))

# rasterize antivenom surface at the same spatial res as the 1km pop
antivenom <- disaggregate(antivenom, fact = c(5,5))

# extend world pop to the same extent as the spp richness raster
pop_dens <- extend(pop_dens, species_richness, value = NA)

# get a count of the maximum number of species in the richness surface
spp_vals <- values(species_richness)
n_spp <- max(spp_vals)

# initialize the cluster (50 cores)
registerDoMC(50)

### generate population living in areas suitable for x or more med spp
par_stack <- foreach(i = 1:n_spp) %dopar% {
  
  # convert species richness surface into a binary surface (`1` = presence of `i` quantity of venomous
  # snake species, `NoData` = absence)
  species_presence <- species_richness
  species_presence[species_presence < i ] <- 0
  species_presence[species_presence  >= i] <- 1
  
  # multiply the population by the binary presence/absence surface
  presence_par <- overlay(pop_dens, species_presence, fun = function(pop_dens, species_presence){
                                                                    (pop_dens*species_presence)})
  
  return(presence_par)
  
}

names(par_stack) <- paste0("par_exposure_", 1:n_spp)
list2env(par_stack, envir = .GlobalEnv)

# get a count of the maximum number of species in the antivenom naive surface
naive_spp <- values(antivenom)
n_naive <- max(naive_spp)

# initialize the cluster (50 cores)
registerDoMC(50)

### generate population living in areas suitable for x or more med spp
antiv_par_stack <- foreach(i = 1:n_naive) %dopar% {
  
  # convert antivenom naive surface into a binary surface (`1` = presence of `i` quantity of therapy
  # naive snake species, `NoData` = absence)
  antiv_n_naive <- antivenom
  antiv_n_naive[antiv_n_naive < i ] <- 0
  antiv_n_naive[antiv_n_naive  >= i] <- 1
  
  # multiply the population by the binary presence/absence antivenom surface
  presence_par <- overlay(pop_dens, species_presence, fun = function(pop_dens, species_presence){
                                                                    (pop_dens*species_presence)})
  
  return(presence_par)
  
}

names(antiv_par_stack) <- paste0("antivenom_naive_", 1:n_spp)
list2env(antiv_par_stack, envir = .GlobalEnv)

# generate global PAR totals for surface
# loop through each PAR raster, and get a sum of the pixel values to generate a global
# PAR estimate
par_totals <- for(i in 1:n_spp) {
  
  # get surface 
  obj_name <- paste0("par_exposure_", i)
  
  # get cell values across surface
  par_exp_vals <- values(obj_name)

  # sum cell values across surface
  par_totals <- sum(par_exp_vals)
  
  return(par_totals)
  
}

# generate global PAR totals for antivenom naive
anti_par_totals <- for(i in 1:n_naive) {
  
  # get surface 
  obj_name <- paste0("antivenom_naive_", i)
  
  # get cell values across surface
  par_naive_vals <- values(obj_name)
  
  # sum cell values across surface
  par_totals <- sum(par_naive_vals)
  
  return(par_totals)
  
}

# generate national estimates

