# script to generate SI plots
# this script depends on all previous plots being archived within the 'archive' folder
# in each sub-directory

# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(png, grid, ggplot2, gridExtra)

# SI figure 1 file:
# this SI figure contains plots of:
# A) Boostrapped MESS 
# B) 95% threshold MESS
# C) 95% threshold MESS with occurrence records placed on top
# D) Recommended modified range

# list the SI figure 1 files
species_list <- list.files('Z:/users/joshua/Snakebite/output/species_mess_maps/mess_png', 
                           pattern = '.*[.]png',
                           full.names = TRUE)

# specify filepath to save combined output
mess_si_name <- paste('Z:/users/joshua/Snakebite/output/species_mess_maps/merged_si_files/Species_MESS_SI_', Sys.Date(), '.pdf', sep = "")

# open all png files
plots <- lapply(species_list,
                function(x){
                img <- as.raster(readPNG(x))
                rasterGrob(img, interpolate = FALSE)
            
              }
)

# merge into one single PDF
ggsave(mess_si_name, 
       width = 8.27,
       height = 11.29,
       marrangeGrob(grobs = plots, nrow = 2, ncol = 1, top = NULL))

rm(plots)
gc()

# SI figure 2 file:
# this SI figure contains plots of:
# A) Boostrapped MESS 
# B) 95% threshold MESS with occurrence records placed on top
# C) 90% threshold MESS with occurrence records placed on top
# D) 75% threshold MESS with occurrence records placed on top

# list the SI figure 2 files
threshold_list <- list.files('Z:/users/joshua/Snakebite/output/species_mess_maps/mess_si/', 
                           pattern = '.*[.]png',
                           full.names = TRUE)

# specify filepath to save combined output
threshold_si_name <- paste('Z:/users/joshua/Snakebite/output/species_mess_maps/merged_si_files/Species_threshold_SI_', Sys.Date(), '.pdf', sep = "")

# open all png files
plots <- lapply(threshold_list,
                function(x){
                img <- as.raster(readPNG(x))
                rasterGrob(img, interpolate = FALSE)
                  
                }
)

# merge into one single PDF
ggsave(threshold_si_name, 
       width = 8.27,
       height = 11.29,
       marrangeGrob(grobs = plots, nrow = 2, ncol = 1, top = NULL))

rm(plots)
gc()

# SI figure 3 file:
# this SI figure contains plots of:
# 1) Maximum covariate values (across 1000 bootstraps of the reference data)
# 2) Mean covariate values (across 1000 bootstraps of the reference data)
# 3) Minimum covariate values (across 1000 bootstraps of the reference data)

covariate_plot_list <- list.files('Z:/users/joshua/Snakebite/output/species_mess_maps/distribution_plots/', 
                           pattern = '.*[.]png',
                           full.names = TRUE)

# specify filepath to save combined output
covariate_si_name <- paste('Z:/users/joshua/Snakebite/output/species_mess_maps/merged_si_files/Species_covariate_variance_SI_', Sys.Date(), '.pdf', sep = "")

plots <- lapply(covariate_plot_list,
                function(x){
                img <- as.raster(readPNG(x))
                rasterGrob(img, interpolate = FALSE)
                  
                }
)

ggsave(covariate_si_name, 
       width = 8.27,
       height = 11.29,
       marrangeGrob(grobs = plots, nrow = 3, ncol = 1, top = NULL))

rm(plots)
gc()

# SI figure 4 file:
# this SI figure contains plots of:
# 1) Proportion correctly classified

pcc_plot_list <- list.files('Z:/users/joshua/Snakebite/output/species_mess_maps/mess_evaluation/', 
                                  pattern = '.*[.]png',
                                  full.names = TRUE)

# specify filepath to save combined output
pcc_si_name <- paste('Z:/users/joshua/Snakebite/output/species_mess_maps/merged_si_files/Species_PCC_SI_', Sys.Date(), '.pdf', sep = "")

plots <- lapply(pcc_plot_list,
                function(x){
                  img <- as.raster(readPNG(x))
                  rasterGrob(img, interpolate = FALSE)
                  
                }
)

ggsave(pcc_si_name, 
       width = 8.27,
       height = 11.29,
       marrangeGrob(grobs = plots, nrow = 2, ncol = 1, top = NULL))

rm(plots)
gc()
