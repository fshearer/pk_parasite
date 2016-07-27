# prepare parasite occurrence data for Malaysia/Brunei/Singapore model

# clear workspace
rm(list = ls())

# ~~~~~~~~~~~~
# load data

# raw occurrence data
dat <- read.csv('data/raw/occurrence/model parasite data OCT2015.csv',
                stringsAsFactors = FALSE)


# ~~~~~~~~~~~~
# sort data

# remove unnecessary columns
v <- c('No.gen..tested',
       'NGT.calc.',
       'Community.tested',
       'samples.collected.to.look.for.Pk..1..or.opportunistic..2..or.previously.confirmed.Pk..3.',
       'Unconfirmed..1.',
       'Country_id',
       'Site_name',
       'Subnational_area',
       'Centroid.Latitude',
       'Centroid.Longitude',
       'EndNote.ID',
       'Abridged_citation',
       'Other_source')

dat <- dat[,-which((names(dat)%in% v))] 

# remove extra columns un-named columns in dataset 
dat <- subset(dat, select= -c(12:25))

# subset malaysia infection data
countries <- c('Malaysia', 'Brunei', 'Singapore')
dat_mbs <- dat[dat$Country %in% countries,]

write.csv(dat_mbs,
          file = 'data/clean/occurrence/parasite_data_mbs.csv',
          row.names = FALSE)


# subset point data
dat_point <- subset(dat_mbs, dat_mbs$Geometry_type!='polygon')

# subset polygon data
poly1_dat <- subset(dat_mbs, dat_mbs$Admin_level=='para_ras1')
poly2_dat <- subset(dat_mbs, dat_mbs$Admin_level=='para_ras2')
poly3_dat <- subset(dat_mbs, dat_mbs$Admin_level=='Admin0') # no points
poly4_dat <- subset(dat_mbs, dat_mbs$Admin_level=='Admin1')
poly5_dat <- subset(dat_mbs, dat_mbs$Admin_level=='Admin2')

# load polygon files
poly1 <- brick('data/raw/raster/polygon_rasters/para_ras1.tif')
poly2 <- brick('data/raw/raster/polygon_rasters//para_ras2.tif')
#admin0 <- brick('data/raw/raster/polygon_rasters/Admin0_Raster.flt')
admin1 <- brick('data/raw/raster/polygon_rasters/Admin1_Raster.flt')
admin2 <- brick('data/raw/raster/polygon_rasters/Admin2_Raster.flt')

# load a 5km template raster of the study area extent
template <- raster('data/clean/raster/mbs_mask.grd')

# load world pop raster as some polygons are outside the SEAsia extent
human_pop <- raster('~/Dropbox/DPHIL/mapping/covariates/worldpop_gpwv4_mosaic_export_5k_MG_Reallocated_filled.tif')

# crop admin0, admin1 and admin2 raster to this extent
#poly3 <- crop(admin0, template)
poly4 <- crop(admin1, template)
poly5 <- crop(admin2, template)

# ~~~~~~~~~~~


# create empty data frame to store info for polygons over 1000 pixels  
big_polys <- dat_point[,colnames(dat_point)][0,]


# sample random points from within polygons on a 5km grid
for (polygon in c('poly1_dat', 'poly2_dat','poly4_dat', 'poly5_dat')){
  
  # get data subset for polygon
  poly_dat <- get(polygon)
  
  for (i in 1:nrow(poly_dat)){
    
    message(paste('processing point', i, 'of', nrow(poly_dat)))
    
    # get polygon code from dataset
    if (polygon %in% c('poly1_dat', 'poly2_dat')) {
      polycode <- poly_dat$Polygon_code[i]
      
    } else {
      polycode <- poly_dat$Gaul_code[i]
    }  
    
    # make sure polycode is numeric
    polycode <- as.numeric(polycode)
    
    # get raster for this polygon
    substring <- substr(polygon, 0, 5)
    raster <- get(substring)
    
    # obtain index for pixels in raster with this identifier (NB this does not exclude pixels containing NAs)
    idx <- which(getValues(raster)==polycode)
    
    # get coordinates of the cells there
    pts <- xyFromCell(raster, idx)
  
    n <- nrow(pts)
    
    # only sample from polygons containing less than 1000 pixels
    if (!(n > 1000)){
      
      # set pixels in this raster without this identifier to NA
      #raster[raster != polycode] <- NA  
      
      # create sampling raster
      # crop and mask human pop raster by polygon raster
      #sampling_ras <- crop(human_pop, raster)
      #sampling_ras <- mask(sampling_ras, raster)  
      
      # sample possible pixels 100 times, weighted by human pop raster
      #pts<- bgSample(sampling_ras, n = 100, prob= TRUE, replace = TRUE, spatial=FALSE)
      
      #if there's more than 100 of them, pick 100 at random
      if (n>100){
        pts <- pts[sample(1:n, size = 100, replace=TRUE),]
        n <- 100
      } 
      
      # pull out the info for that record
      info_mat <- poly_dat[rep(i, nrow(pts)),]
      
      # stick the coordinates in
      info_mat[,c('Latitude', 'Longitude')] <- pts[,2:1]
      
      # check that new points are within covariate mask
      info_mat <- insideMask(info_mat, template)
      
      # append it to point data
      dat_point <- rbind(dat_point, info_mat)
      
    } else {
      
      # pull out info for that record
      info_mat <- poly_dat[i,]
      
      # make a dataframe of codes for polygons containing more than 1000 pixels
      big_polys <- rbind(big_polys, info_mat)
      
    }
    
  }
  
}

# create dataframe of points and polygons not within big_polys dataset
big_polys_idx <- which(dat_mbs$ID %in% big_polys$ID)
dat_mbs_polys_excluded <- dat_mbs[-big_polys_idx,]

data_years <- dat_mbs_polys_excluded$Year

outside_range <- length(which((data_years > 2012) 
                              | (data_years < 2001))) / length(data_years)*100


# output the resulting tables
write.csv(dat_mbs_polys_excluded,
          file = 'data/clean/occurrence/parasite_data_mbs_bigpolys_excluded.csv',
          row.names = FALSE)

write.csv(dat_point,
          file = 'data/clean/occurrence/polygon_data_mbs.csv',
          row.names = FALSE)

write.csv(big_polys,
          file='data/raw/occurrence/excluded_polygons_mbs.csv',
          row.names = FALSE)




