# organise covariates for South East Asia parasite models - just need covs_current to make prediction

# clear working space
rm(list=ls())

# list the covariates of interest
covars <- c('~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/tasselled cap brightness/TCB_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/tasselled cap brightness/TCB_SD_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//worldpop_gpwv4_mosaic_export_5k_MG_Reallocated_filled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/enhanced vegetation index/EVI_Fixed_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/enhanced vegetation index/EVI_Fixed_SD_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/tasselled cap wetness/TCW_Mean_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/tasselled cap wetness/TCW_SD_5km_Mean_DEFLATE.ABRAID_Extent.Gapfilled.tif',
            '~/Dropbox/DPHIL/mapping/covariates//SRTM elevation/mod_dem_5k_CLEAN.flt',
            '~/Dropbox/DPHIL/mapping/covariates//access.asc',
            '~/Dropbox/DPHIL/mapping/covariates//tempaucpf_5k.tif',
            '~/Dropbox/DPHIL/mapping/covariates/forest_extent_SEAsia/ifl_01to12_5km/ifl2012.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates/forest_extent_SEAsia/plant_01to12_5km/plt2012.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates/MODIS satellite data/landcover 2001-12/2012.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif')

# give them readable names
names <- c('TCB_mean', 
           'TCB_SD',
           'human_pop', 
           'EVI_mean',
           'EVI_SD',
           'TCW_mean',
           'TCW_SD',
           'SRTM_elevation',
           'urban_access',
           'Pf_temp',
           'forest_intact',
           'forest_disturbed',
           'open_shrublands',
           'woody_savannas',
           'savannas',
           'grasslands',
           'permanent_wetlands',
           'croplands',
           'cropland_natural_vegetation_mosaic')

# loop through opening links to the rasters
covs <- lapply(covars, raster)

# load SEAsia human pop raster to get the extent
extent_template <- raster('~/Dropbox/DPHIL/mapping/covariates/SEAsia_pop5k/SEAsia_pop5k.tif')
ext <- extent(extent_template)

# get all extents
#extents <- t(sapply(covs, function (x) as.vector(extent(x))))

# get the smallest extent squares of all layers
#ext <- extent(c(max(extents[, 1]),
#                min(extents[, 2]),
#                max(extents[, 3]),
#                min(extents[, 4])))

# crop all layers by this
covs <- lapply(covs, crop, ext)

# stack the rasters
covs <- stack(covs) 

cbind(names(covs), names)

warning('check the names match up!')

# give them nicer names
names(covs) <- names

# set all NAs in Pf_temp layer to 0
covs[['Pf_temp']][is.na(covs[['Pf_temp']])] <- 0

# load annual reservoir and vector covs
annual_covs <- brick('data/raw/raster/covariates/annual_reservoirs/annual_reservoirs_covs.grd')

# select the most current distribution for fascicularis, nemestrina and the leucosphyrus group
annual_covs <- annual_covs[[c('fascicularis_2012', 'nemestrina_2012', 'leucosphyrus_group_2012')]]

# rename layers
names(annual_covs)[names(annual_covs)=='fascicularis_2012'] <- 'fascicularis'
names(annual_covs)[names(annual_covs)=='nemestrina_2012'] <- 'nemestrina'
names(annual_covs)[names(annual_covs)=='leucosphyrus_group_2012'] <- 'leucosphyrus_group'

# add to covariate stack
covs <- addLayer(covs, annual_covs)

# mask this layer by the human population layer
#covs[[11]] <- mask(covs[[11]], covs[[3]])

# load P falciparum temperature suitability layer
#Pf_temp <- raster('covariates/tempaucpf_5k.tif')

# crop by extent
#Pf_temp <- crop(Pf_temp, ext)

# set all NAs to 0
#Pf_temp[is.na(Pf_temp)] <- 0

# mask by SE Asia population layer
#Pf_temp <- mask(Pf_temp, covs[[3]])

# add to covs stack
#covs <- addLayer(covs, Pf_temp)

# rename
#names(covs)[names(covs)=='tempaucpf_5k'] <- 'Pf_temp'

# load annual species covariates
#reservoir_covs <- brick('data/raw/raster/covariates/annual_reservoirs/annual_reservoirs_covs.grd')

# drop leonina, dirus_complex and leucosphyrus group layers
#reservoir_covs <- dropLayer(reservoir_covs, c(13:24, 37:48, 61:72))

# drop temporal layers
#reservoir_covs <- dropLayer(reservoir_covs, c(1:11, 13:23, 25:35))

# crop all layers by the extent
#reservoir_covs <- crop(reservoir_covs, ext)

# add reservoir covs to other covs
#covs <- addLayer(covs, reservoir_covs)

# load admin0 raster file
#admin0_ras <- raster('data/raw/raster/polygon_rasters/Admin0_Raster.flt')

# crop by extent
#study_extent <- crop(admin0_ras, ext)

# change all values that are not NA to 1
#study_extent[!(is.na(study_extent))] <- 1

# create master mask
master_mask <- masterMask(covs)

#mask all covariates by master mask
covs <- mask(covs, master_mask)

# mask study extent by master mask
#study_extent <- mask(study_extent, master_mask)

# save study extent
writeRaster(master_mask,
            file = 'data/clean/raster/SEAsia_extent',
            overwrite=TRUE)

# output covs as multiband .grd files
writeRaster(covs,
            file = 'data/clean/raster/SEAsia_covs',
            overwrite = TRUE)

