# clear workspace
rm(list = ls())

# load libraries
require(pdftools)

# list files
list <- as.list(list.files("Z:/users/joshua/Snakebite/WHO EOR maps/"))

# loop through and convert pdf files into png
for(i in 1:length(list)){
  
  # get file name
  file_name <- unlist(list[i])
  
  # get species name
  spp <- gsub('.pdf$', '', file_name)
  
  # gen path to pdf
  paths <- paste("Z:/users/joshua/Snakebite/WHO EOR maps/",
                 file_name,
                 sep = "")
  
  # open file and render to png
  bitmap <- pdf_render_page(paths, page = 1)
  
  # generate outpath
  outpath <- paste("Z:/users/joshua/Snakebite/WHO EOR maps/",
                   spp,
                   ".png",
                   sep = "")
  
  # save bitmap image
  png::writePNG(bitmap, outpath)
  
}
