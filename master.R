# master script
# script to loop through each separate file, and generate surfaces

# convert shapes to raster
source('Z:/users/joshua/Snakebite/snakebite/convert_shp_to_raster.R')

# update plots
source('Z:/users/joshua/Snakebite/snakebite/plot_occurrence.R')

# generate species richness
source('Z:/users/joshua/Snakebite/snakebite/gen_species_richness.R')

# generate antivenom coverage
source('Z:/users/joshua/Snakebite/snakebite/gen_antivenom_coverage.R')
