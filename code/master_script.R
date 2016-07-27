# script to run all analysis

# clear workspace
rm(list = ls())

# ~~~~~~~~~~~~
# load packages
library(seegSDM)
library(raster)
library(snowfall)
library(dismo)
library(rgdal)

# ~~~~~~~~~~~~~
# load functions file
source('code/functions_parasite.R')

# ~~~~~~~~~~~~
# prepare model inputs

# organise covariate layers for MBS fitting
source('code/sort_rasters_MBS.R')

# organise covariate layers for SE Asia prediction
source('code/sort_rasters_parasite.R')

# organise occurrence data 
source('code/sorting_polygons_MBS.R')

# ~~~~~~~~~~~~
# run models, calculate stats, create MESS map and plot outputs
source('code/run_parasite_MBS.R')

