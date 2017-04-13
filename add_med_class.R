# script to add medical class to shapefiles
# clear workspace
rm(list = ls())

# load required packages
require(foreign)

# set working directory
setwd('Z:/users/joshua/Snakebite/WHO EOR shapefiles/')

shape_list <- list.files('Z:/users/joshua/Snakebite/WHO EOR shapefiles/',
                         pattern = '*dbf$')

class_list <- read.csv('Z:/users/joshua/Snakebite/snakebite/snake_list.csv',
                       stringsAsFactors = FALSE)

# loop through and add medical class to attributes
for(i in 1:length(shape_list)){
  
  # get file name
  file_name <- unlist(shape_list[i])
  
  # get species med cat
  # sub out dbf extension
  spp <- gsub('.dbf$', '', file_name)
  
  # inform progress
  message(paste('Adding attribute to', spp, '[file', i,  'of', length(shape_list), ']', sep = " "))
  
  # subset to get class
  sub <- class_list[class_list$split_spp == spp, ]
  cat <- as.numeric(unique(sub$majority_cat))
  
  # read in dbf
  snake_shape <- read.dbf(file_name, as.is = TRUE)
  
  # add new attribute
  snake_shape$Med_Class <- rep(cat, nrow(snake_shape))
  
  # write out dbf
  write.dbf(snake_shape, file_name)
  
}
