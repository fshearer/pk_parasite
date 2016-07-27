# organise covariates for Malaysia/Brunei/Singapore knowlesi model

# clear working space
rm(list=ls())

# load admin shapefiles
admin0 <- shapefile('data/raw/raster/admin/admin2013_0.shp')
admin1 <- shapefile('data/raw/raster/admin/admin2013_1.shp')

# get Admin1 gaul codes for Malaysia, Brunei and Singapore
country_codes  <- c('MYS', 'BRN', 'SGP')
codes <- admin1$GAUL_CODE[admin1$COUNTRY_ID %in% country_codes]

# get study extent
mbs <- admin0[admin0$COUNTRY_ID %in% country_codes,]
ext <- extent(mbs)

# save Malaysia extent file
shapefile(mbs, 'output/MBS/mbs_extent.shp', overwrite = TRUE)

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
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2001.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2002.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2003.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2004.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2005.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2006.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2007.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2008.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2009.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2010.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2011.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/ifl_01to12_5km/ifl2012.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2001.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2002.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2003.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2004.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2005.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2006.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2007.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2008.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2009.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2010.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2011.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//forest_extent_SEAsia/plant_01to12_5km/plt2012.fractional.5km.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class07_Open_Shrublands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class08_Woody_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class09_Savannas.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class10_Grasslands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class11_Permanent_Wetlands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class12_Croplands.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2001.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2002.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2003.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2004.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2005.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2006.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2007.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2008.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2009.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2010.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2011.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif',
            '~/Dropbox/DPHIL/mapping/covariates//MODIS satellite data/landcover 2001-12/2012.Class14_Cropland_Natural_Vegetation_Mosaic.5km.Percentage.tif')
            
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
           'forest_intact_2001',
            'forest_intact_2002',
            'forest_intact_2003',
            'forest_intact_2004',
            'forest_intact_2005',
            'forest_intact_2006',
            'forest_intact_2007',
            'forest_intact_2008',
            'forest_intact_2009',
            'forest_intact_2010',
            'forest_intact_2011',
            'forest_intact_2012',
            'forest_disturbed_2001',
            'forest_disturbed_2002',
            'forest_disturbed_2003',
            'forest_disturbed_2004',
            'forest_disturbed_2005',
            'forest_disturbed_2006',
            'forest_disturbed_2007',
            'forest_disturbed_2008',
            'forest_disturbed_2009',
            'forest_disturbed_2010',
            'forest_disturbed_2011',
            'forest_disturbed_2012',
            'open_shrublands_2001',
            'open_shrublands_2002',
            'open_shrublands_2003',
            'open_shrublands_2004',
            'open_shrublands_2005',
            'open_shrublands_2006',
            'open_shrublands_2007',
            'open_shrublands_2008',
            'open_shrublands_2009',
            'open_shrublands_2010',
            'open_shrublands_2011',
            'open_shrublands_2012',
            'woody_savannas_2001',
            'woody_savannas_2002',
            'woody_savannas_2003',
            'woody_savannas_2004',
            'woody_savannas_2005',
            'woody_savannas_2006',
            'woody_savannas_2007',
            'woody_savannas_2008',
            'woody_savannas_2009',
            'woody_savannas_2010',
            'woody_savannas_2011',
            'woody_savannas_2012',
            'savannas_2001',
            'savannas_2002',
            'savannas_2003',
            'savannas_2004',
            'savannas_2005',
            'savannas_2006',
            'savannas_2007',
            'savannas_2008',
            'savannas_2009',
            'savannas_2010',
            'savannas_2011',
            'savannas_2012',
            'grasslands_2001',
            'grasslands_2002',
            'grasslands_2003',
            'grasslands_2004',
            'grasslands_2005',
            'grasslands_2006',
            'grasslands_2007',
            'grasslands_2008',
            'grasslands_2009',
            'grasslands_2010',
            'grasslands_2011',
            'grasslands_2012',
            'permanent_wetlands_2001',
            'permanent_wetlands_2002',
            'permanent_wetlands_2003',
            'permanent_wetlands_2004',
            'permanent_wetlands_2005',
            'permanent_wetlands_2006',
            'permanent_wetlands_2007',
            'permanent_wetlands_2008',
            'permanent_wetlands_2009',
            'permanent_wetlands_2010',
            'permanent_wetlands_2011',
            'permanent_wetlands_2012',
            'croplands_2001',
            'croplands_2002',
            'croplands_2003',
            'croplands_2004',
            'croplands_2005',
            'croplands_2006',
            'croplands_2007',
            'croplands_2008',
            'croplands_2009',
            'croplands_2010',
            'croplands_2011',
            'croplands_2012',
            'cropland_natural_vegetation_mosaic_2001',
            'cropland_natural_vegetation_mosaic_2002',
            'cropland_natural_vegetation_mosaic_2003',
            'cropland_natural_vegetation_mosaic_2004',
            'cropland_natural_vegetation_mosaic_2005',
            'cropland_natural_vegetation_mosaic_2006',
            'cropland_natural_vegetation_mosaic_2007',
            'cropland_natural_vegetation_mosaic_2008',
            'cropland_natural_vegetation_mosaic_2009',
            'cropland_natural_vegetation_mosaic_2010',
            'cropland_natural_vegetation_mosaic_2011',
            'cropland_natural_vegetation_mosaic_2012')

# loop through opening links to the rasters
covs <- lapply(covars, raster)

# crop all layers by malaysia extent
covs <- lapply(covs, crop, ext)

# stack the rasters
covs <- stack(covs) 

cbind(names(covs), names)

warning('check the names match up!')

# give them nicer names
names(covs) <- names

#load annual species covariates
reservoir_covs <- brick('data/raw/raster/covariates/annual_reservoirs/annual_reservoirs_covs.grd')

# crop all layers by malaysia extent
reservoir_covs <- crop(reservoir_covs, ext)

# add reservoir covs to other covs
covs <- addLayer(covs, reservoir_covs)

# create study mask raster layer with 1s for cells within Malaysia/Brunei/Singapore and NAs elsewhere
admin1_ras <- raster('data/raw/raster/polygon_rasters/Admin1_Raster.flt')
mbs_ext <- crop(admin1_ras, ext)
mbs_ext[! getValues(mbs_ext) %in% codes] <- NA
mbs_mask <- mbs_ext * 0 + 1

# add Malaysia extent layer and admin1 layer to covs and create master mask 
covs <- addLayer(covs, mbs_ext)
covs <- addLayer(covs, mbs_mask)
master_mask <- masterMask(covs)

# mask all covariates by master mask
covs <- mask(covs, master_mask)

# output malaysia admin1 file
writeRaster(covs[[155]],
            file = 'data/clean/raster/mbs_admin1',
            overwrite = TRUE)

# output malaysia extent file
writeRaster(covs[[156]],
            file = 'data/clean/raster/mbs_mask',
            overwrite = TRUE)

# then drop Malaysia extent and mask layers
covs <- dropLayer(covs, 155:156)

# create subsets of temporal, non-temporal and the most contemporary covariates for running models 
covs_current <- subset(covs, c(1:10, 22, 34, 46, 58, 70, 82, 94, 106, 118, 130, 142, 154))
covs_temporal <- subset(covs, 11:154)
covs_nontemporal <- subset(covs, 1:10)

# output them as multiband .grd files
writeRaster(covs_current,
            file = 'data/clean/raster/mbs_raster_current',
            overwrite = TRUE)

writeRaster(covs_temporal,
             file = 'data/clean/raster/mbs_raster_temporal',
             overwrite = TRUE)
 
writeRaster(covs_nontemporal,
             file = 'data/clean/raster/mbs_raster_nontemporal',
             overwrite = TRUE)
