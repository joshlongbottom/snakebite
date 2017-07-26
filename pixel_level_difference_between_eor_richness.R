# clear workspace
rm(list = ls())

# load required packages
require(raster)

# generate pixel level difference surface between original EOR and modified EOR
# read in the original EOR
original <- raster('Z:/users/joshua/Snakebite/output/species_richness/original_eor_combined_categories_2017-07-25.tif')

# read in modified EOR
modified <- raster('Z:/users/joshua/Snakebite/output/species_richness/modified_eor_combined_categories_2017-07-25.tif')

# calculate pixel level difference
difference <- modified - original

# check pixel values
table(values(difference))

# one pixel has a value of `-1`, change this to `NA`
difference[difference == -1] <- NA

# write out pixel level difference
writeRaster(difference,
            file = 'Z:/users/joshua/Snakebite/output/species_richness/pixel_level_difference_combined_categories_2017-07-26.tif',
            format = 'GTiff',
            overwrite = TRUE)

