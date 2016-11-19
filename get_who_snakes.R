# script to download snake expert opinion maps from WHO
# clear workspace
rm(list = ls())

require(utils)

# get snakes (list of snakes species on WHO website, as csv)
# list was compiled from copying source data
snakes <- read.csv("Z:/users/joshua/Snakebite/snakebite/snake_list.csv",
                   stringsAsFactors = FALSE)

# list spp
species_list <- as.list(unique(snakes$species))

# loop through and grab pdfs from website
for(i in 1:length(species_list)){
  
  # get species name
  spp_name <- unlist(species_list[i])
  
  # get species url
  spp_url <- paste("http://apps.who.int/bloodproducts/snakeantivenoms/database/Images/SnakesDistribution/Large/map_",
                   spp_name,
                   ".pdf",
                   sep = "")
  
  # specify outpath
  outpath <- paste("Z:/users/joshua/Snakebite/WHO EOR maps/",
                   spp_name,
                   ".pdf",
                   sep = "")
  # download file
  download.file(spp_url,
                outpath,
                mode = "wb")
  
}