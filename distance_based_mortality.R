# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(raster)

# load in accessibility surface
accessibility <- raster('Z:/users/joshua/Snakebite/rasters/accessibility/accessibility_50k+_2017-01-05_final.tif')

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
writeRaster(accessibility_mortality, 
            file = 'Z:/users/joshua/Snakebite/output/population_at_risk/distance_based_mortality_raw_percentage',
            format = 'GTiff',
            overwrite = TRUE)

