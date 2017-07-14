# script to download snake expert opinion maps from WHO
# clear workspace
rm(list = ls())

pacman::p_load(utils, pdftools)

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
  outpath <- paste("Z:/users/joshua/Snakebite/WHO EOR PDFs/",
                   spp_name,
                   ".pdf",
                   sep = "")
  # download file
  download.file(spp_url,
                outpath,
                mode = "wb")
  
}

# list the downloaded pdf files
list <- as.list(list.files("Z:/users/joshua/Snakebite/WHO EOR PDFs/", pattern = '.pdf$'))

# loop through and convert these pdf files into png
for(i in 1:length(list)){
  
  # get file name
  file_name <- unlist(list[i])
  
  # get species name
  spp <- gsub('.pdf$', '', file_name)
  
  # gen path to pdf
  paths <- paste("Z:/users/joshua/Snakebite/WHO EOR PDFs/",
                 file_name,
                 sep = "")
  
  # open file and render to png
  bitmap <- pdf_render_page(paths, page = 1, dpi = 300)
  
  # generate outpath
  outpath <- paste("Z:/users/joshua/Snakebite/WHO EOR PNG/",
                   spp,
                   ".png",
                   sep = "")
  
  # save bitmap image
  png::writePNG(bitmap, outpath)
  
}
