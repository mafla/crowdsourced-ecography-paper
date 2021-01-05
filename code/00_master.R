rm(list=ls())
setwd("~/work/FLORKART/Mahecha-etal-CrowdSourced/FloraIncDimRed-master")

# list of packages we need for the analytics
packages2use = c(
  # fordata wrangling
  "tidyverse",
  # for some plots
  "reshape2",
  "gridExtra",
  # for ecological distance computations
  "vegan",
  # for dimensionality reduction
  "dimRed",
  "igraph",
  # for cca
  "yacca",
  # for clustering
  "apcluster",
  # for spatial machine learning
  "CAST",
  "caret",
  "randomForest",
  # for PCA (not that this has also a CCA, so we need to take care to yacca::cca)
  "FactoMineR",
  "factoextra",
  # For parallel computing
  "parallel",
  # for geo processing
  "rgdal",
  "sp",
  "rgeos",
  # for plotting
  "corrplot",
  "pals",
  "wordcloud",
   # for producing and saving maps
  "mapview",
  "htmlwidgets",
  "webshot",
  "leaflet")
  
for(p in packages2use){
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}


# noes for parallel computations
n_nodes = min(40, detectCores()-1)

# wrapper function for generating the maps
cat("====== load extra code =======\n")
source("f_map2mtb.R")
# function for computing k-NN graph and derive geodesic distance matrix
source("f_d2dgeo.R")
# function to truncate a vector by quantiles
source("f_trunc_quant.R")

##-------------
# read data
cat("====== run 01_read_data.R =======\n")
source("01_read_data.R")

##-------------
# set to FALSE if you have run this code before and only want to work on plots etc...
# this would then try to load preestiated intermediate results
calc_from_scratch = FALSE

##-------------
# comparison of florkart and flora incognita
cat("====== run 02_compare_fk_fi.R =======\n")
source("02_compare_fk_fi.R")

##-------------
# dimensionality reduction
cat("====== run 03_dimred.R =======\n")
source("03_dimred.R")

##-------------
# canonical correlation analysis
cat("====== run 04_cca.R =======\n")
source("04_cca.R")

##-------------
# prediction of reduced dimensions and canonical variates
#considering spatial cross validation issues
## predictors are population, climate and soil data - the latter PCA preprocessed to consider redundancies
cat("====== run 05_predict.R =======\n")
source("05_predict.R")
