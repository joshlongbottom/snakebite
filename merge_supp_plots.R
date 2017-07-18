# script to generate SI plots
# this script depends on all previous plots being archived within the 'archive' folder
# in each sub-directory

# clear workspace
rm(list = ls())

# load required packages
pacman::p_load(png, grid, ggplot2, gridExtra, parallel)

# SI figure 1 file:
# this SI figure contains plots of:
# A) Boostrapped MESS 
# B) 95% threshold MESS
# C) 95% threshold MESS with occurrence records placed on top
# D) Recommended modified range

# list the SI figure 1 files
species_list <- list.files('output/mess_png/', 
                           pattern = '.*[.]png',
                           full.names = TRUE)

# specify filepath to save combined output
mess_si_name <- paste('output/merged_si/Species_MESS_SI_', Sys.Date(), '.pdf', sep = "")

# open all png files
plots <- mclapply(species_list,
                  function(x){
                  img <- as.raster(readPNG(x))
                  rasterGrob(img, interpolate = FALSE)
                  },
                  mc.cores = 50)

# merge into one single PDF
ggsave(mess_si_name, 
       width = 11.69,
       height = 8.27,
       marrangeGrob(grobs = plots, nrow = 1, ncol = 1, top = NULL))

rm(plots)
gc()

# SI figure 2 file:
# this SI figure contains plots of:
# A) Boostrapped MESS 
# B) 95% threshold MESS with occurrence records placed on top
# C) 90% threshold MESS with occurrence records placed on top
# D) 75% threshold MESS with occurrence records placed on top

# list the SI figure 2 files
threshold_list <- list.files('output/mess_si/', 
                           pattern = '.*[.]png',
                           full.names = TRUE)

# specify filepath to save combined output
threshold_si_name <- paste('output/merged_si/Species_threshold_SI_', Sys.Date(), '.pdf', sep = "")

# open all png files
plots <- mclapply(threshold_list,
                  function(x){
                  img <- as.raster(readPNG(x))
                  rasterGrob(img, interpolate = FALSE)
                  },
                  mc.cores = 50)

# merge into one single PDF
ggsave(threshold_si_name, 
       width = 11.69,
       height = 8.27,
       marrangeGrob(grobs = plots, nrow = 2, ncol = 1, top = NULL))

rm(plots)
gc()

# SI figure 3 file:
# this SI figure contains plots of:
# 1) Maximum covariate values (across 1000 bootstraps of the reference data)
# 2) Mean covariate values (across 1000 bootstraps of the reference data)
# 3) Minimum covariate values (across 1000 bootstraps of the reference data)

covariate_plot_list <- list.files('output/distribution_plots/', 
                           pattern = '.*[.]png',
                           full.names = TRUE)

# specify filepath to save combined output
covariate_si_name <- paste('output/merged_si/Species_covariate_variance_SI_', Sys.Date(), '.pdf', sep = "")

plots <- mclapply(covariate_plot_list,
                  function(x){
                  img <- as.raster(readPNG(x))
                  rasterGrob(img, interpolate = FALSE)
                  },
                  mc.cores = 50)

ggsave(covariate_si_name, 
       width = 11.69,
       height = 8.27,
       marrangeGrob(grobs = plots, nrow = 1, ncol = 1, top = NULL))

rm(plots)
gc()

# SI figure 4 file:
# this SI figure contains plots of:
# 1) Proportion correctly classified

pcc_plot_list <- list.files('output/mess_evaluation/', 
                                  pattern = '.*[.]png',
                                  full.names = TRUE)

# specify filepath to save combined output
pcc_si_name <- paste('output/merged_si/Species_PCC_SI_', Sys.Date(), '.pdf', sep = "")

plots <- mclapply(pcc_plot_list,
                  function(x){
                  img <- as.raster(readPNG(x))
                  rasterGrob(img, interpolate = FALSE)
                  },
                  mc.cores = 50)

ggsave(pcc_si_name, 
       width = 8.27,
       height = 11.69,
       marrangeGrob(grobs = plots, nrow = 2, ncol = 1, top = NULL))

rm(plots)
gc()
