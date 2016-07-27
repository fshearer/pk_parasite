# run knowlesi parasite malaysia/brunei/singapore model
# includes human, monkey and vector infection data
# pseudo-absences for monkey infection data were sampled from monkey and mammal presence points
# pseudo-absences for vector infection data were sampled from human pop layer
# uses a joint distribution to assess difference in distribution of parasite in humans, monkeys and vectors 

# clear workspace
rm(list = ls()) 

# set the RNG seed
set.seed(1)

# ~~~~~~~~~~~~~~~~~

# set output path
outpath <- 'output/MBS/'

# ~~~~~~~~~~~~~~~~~
# load raster data

# raster covariates
covs_current <- brick('data/clean/raster/mbs_raster_current.grd')
covs_temporal <- brick('data/clean/raster/mbs_raster_temporal.grd')
covs_nontemporal <- brick('data/clean/raster/mbs_raster_nontemporal.grd')

# drop correlated covariates and reservoir and vector species
covs_current <- dropLayer(covs_current, c('EVI_mean', 'EVI_SD', 'TCB_mean'))
covs_nontemporal <- dropLayer(covs_nontemporal, c('EVI_mean', 'EVI_SD', 'TCB_mean'))

# rename the temporal layers (remove year)
names(covs_current)[names(covs_current)=='forest_intact_2012']<-'forest_intact'
names(covs_current)[names(covs_current)=='forest_disturbed_2012']<-'forest_disturbed'
names(covs_current)[names(covs_current)=='open_shrublands_2012']<-'open_shrublands'
names(covs_current)[names(covs_current)=='woody_savannas_2012']<-'woody_savannas'
names(covs_current)[names(covs_current)=='savannas_2012']<-'savannas'
names(covs_current)[names(covs_current)=='grasslands_2012']<-'grasslands'
names(covs_current)[names(covs_current)=='croplands_2012']<- 'croplands'
names(covs_current)[names(covs_current)=='cropland_natural_vegetation_mosaic_2012']<-'cropland_natural_vegetation_mosaic'
names(covs_current)[names(covs_current)=='permanent_wetlands_2012']<- 'permanent_wetlands'
names(covs_current)[names(covs_current)=='fascicularis_2012']<- 'fascicularis'
names(covs_current)[names(covs_current)=='nemestrina_2012']<- 'nemestrina'
names(covs_current)[names(covs_current)=='leucosphyrus_group_2012']<- 'leucosphyrus_group'

# get the population raster
human_pop <- covs_current[[which(names(covs_current)=='human_pop')]]

# ~~~~~~~~~~~~~~~~~
# load occurrence data

# occurrence data with polygons incorporated 
occ_mbs <- read.csv('data/clean/occurrence/polygon_data_mbs.csv', stringsAsFactors = FALSE)

# correct column entries for mosquito host
occ_mbs[occ_mbs$Host=='mosquitoes with sporozoites',]$Host <- 'mosquito'

# move 5 points in Singapore onto land
# find index of points outside the mask
outside_mask <- which(is.na(extract(human_pop, occ_mbs[,c('Longitude', 'Latitude')])))
outside_points <- occ_mbs[outside_mask,]
outside_points <- outside_points[, c('Longitude', 'Latitude')]
land_points <- nearestLand(outside_points, human_pop, 10000)

# replace all outside_mask points with lats/longs for land points
for (i in 1:length(outside_mask)) {
  occ_mbs[outside_mask[[i]], c('Longitude', 'Latitude')] <- land_points[i,]
}

# change ID column name to Unique_ID
colnames(occ_mbs)[colnames(occ_mbs)=='ID'] <- 'Unique_ID'

# get extent of covs_current
ext <- extent(covs_current)

# find index of all points falling outside the extent
outside_ext_idx <- which((occ_mbs$Latitude < ext[3]) 
                         |(occ_mbs$Latitude > ext[4]) 
                         |(occ_mbs$Longitude < ext[1]) 
                         |(occ_mbs$Longitude > ext[2]))

stopifnot(length(outside_ext_idx)==0)

# subset the presence and absence points 
absence <- subset(occ_mbs, occ_mbs$Presence==0)
occ <- subset(occ_mbs, occ_mbs$Presence==1)

# ~~~~~~~~~~~~~~~~
# sort background data
# load occurrence data without coords for polygons 

occ_raw <- read.csv('data/clean/occurrence/parasite_data_mbs_bigpolys_excluded.csv', stringsAsFactors = FALSE)
occ_raw[occ_raw$Host=='mosquitoes with sporozoites',]$Host <- 'mosquito'

# subset the presence points 
presence <- subset(occ_raw, occ_raw$Presence==1)
absences2 <- subset(occ_raw, occ_raw$Presence==0)

# find out number of infections in humans, mosquitos and monkeys in occurrence dataset
human_pres <- nrow(subset(presence, presence$Host=='human'))
#human_abs <- nrow(subset(absences2, absences2$Host=='human'))
vector_pres <- nrow(subset(presence, presence$Host=='mosquito'))
#vector_abs <- nrow(subset(absences2, absences2$Host=='mosquito'))
monkey_pres <- nrow(subset(presence, presence$Host=='monkey'))
#monkey_abs <- nrow(subset(absences2, absences2$Host=='monkey'))
total_pres <- human_pres + vector_pres + monkey_pres
#total_abs <- human_abs + vector_abs + monkey_abs

# set total number of background points
bg <- 6000

# calculate number of human, vector and monkey points based on total number of background points
human_points <- round(human_pres/total_pres*bg)
vector_points <-round(vector_pres/total_pres*bg)
monkey_points <- round(monkey_pres/total_pres*bg)

# load background datasets
fascicularis <- read.csv('../knowlesi_reservoirs/data/clean/occurrence/fascicularis_presence.csv', stringsAsFactors=FALSE)
nemestrina <- read.csv('../knowlesi_reservoirs/data/clean/occurrence/nemestrina_presence.csv', stringsAsFactors = FALSE)
monkey_bias <- read.csv('../knowlesi_reservoirs/data/clean/occurrence/mammals-bias.csv', stringsAsFactors = FALSE)

# remove extra columns in fascicularis dataset 
fascicularis <- subset(fascicularis, select= -c(9:19))

# change column names of mammal-bias dataset
colnames(monkey_bias)[colnames(monkey_bias)=='decimalLatitude'] <- 'Latitude'
colnames(monkey_bias)[colnames(monkey_bias)=='decimalLongitude'] <- 'Longitude'
colnames(monkey_bias)[colnames(monkey_bias)=='year'] <- 'Year'

# add geometry-type column
monkey_bias$Geometry_type <- 'point'

# use presence data only
fascicularis <- subset(fascicularis, fascicularis$Presence==1)
nemestrina <- subset(nemestrina, nemestrina$Presence==1)

# make column names match
colnames(fascicularis)[colnames(fascicularis)=='Geometry_t'] <- 'Geometry_type'
colnames(nemestrina)[colnames(nemestrina)=='Geometry_t'] <- 'Geometry_type'

# combine monkey datasets
monkeys <- rbind(fascicularis[, c('Longitude', 'Latitude', 'Year', 'Geometry_type')], 
                 nemestrina[, c('Longitude', 'Latitude', 'Year', 'Geometry_type')],
                 monkey_bias[, c('Longitude', 'Latitude', 'Year', 'Geometry_type')]) 

# remove any points that don't fall on covariate mask
template <- covs_current[[1]]
monkeys <- insideMask(monkeys, template)

# sample background points from monkey and vector datasets and combine
monkey_bg <- monkeys[sample(1:nrow(monkeys), monkey_points, replace=FALSE),]
monkey_bg$Host <- 'monkey'

# generate pseudo-absence points for human and vector infection points
human_bg <- bgSample (human_pop, n=human_points, prob=TRUE, replace= TRUE, spatial=FALSE)
vector_bg <- bgSample (human_pop, n=vector_points, prob=TRUE, replace=TRUE, spatial=FALSE)

# convert to a dataframe and add relevant columns 
human_background <- data.frame(human_bg, stringsAsFactors = FALSE)
colnames(human_background)[colnames(human_background)=='x'] <- 'Longitude'
colnames(human_background)[colnames(human_background)=='y'] <- 'Latitude'
human_background$Year <- 2012
human_background$Geometry_type <- 'point'
human_background$Host <- 'human'

vector_background <- data.frame(vector_bg, stringsAsFactors = FALSE)
colnames(vector_background)[colnames(vector_background)=='x'] <- 'Longitude'
colnames(vector_background)[colnames(vector_background)=='y'] <- 'Latitude'
vector_background$Year <- 2012
vector_background$Geometry_type <- 'point'
vector_background$Host <- 'mosquito'

background <- rbind(human_background, vector_background, monkey_bg)

# add Unique_ID column for later to deal with polygon data in occurrence data set
background$Unique_ID <- NA

# combine the occurrence and background records, exclude "true" absence records
dat <- rbind(cbind(PA = rep(1, nrow(occ)),
                   wt = rep(1, nrow(occ)),
                   occ[, c('Unique_ID', 'Longitude', 'Latitude', 'Year', 'Geometry_type', 'Host')]),
             cbind(PA = rep(0, nrow(background)),
                   wt = 0.5,
                   background[ ,c('Unique_ID', 'Longitude', 'Latitude', 'Year','Geometry_type', 'Host')]))

# ~~~~~~~~~~~~~~~~~

# extract covariate values for each data point

# create an empty dataframe
dat_covs <- extract(covs_current, dat[, c('Longitude', 'Latitude')])

dat_covs[] <- NA

# loop through extracting covs for all other years
for (year in 2001:2012) {
  
  # truncate years to the ones we have covariates for
  years <- dat$Year
  years <- pmax(years, 2001)
  years <- pmin(years, 2012)
  
  # index for this year's data
  idx <- which(years == year)
  
  # covs for this year
  
  # get index for temporal covariates for relevant year
  
  idx2 <- which (substr(names(covs_temporal),
                        nchar(names(covs_temporal)) - 3,
                        nchar(names(covs_temporal)))==year)
  
  
  covs_year <- covs_temporal[[idx2]]
  names(covs_year)[names(covs_year)==paste0('forest_intact_', year)]<-'forest_intact'
  names(covs_year)[names(covs_year)==paste0('forest_disturbed_', year)]<-'forest_disturbed'
  names(covs_year)[names(covs_year)==paste0('open_shrublands_', year)]<-'open_shrublands'
  names(covs_year)[names(covs_year)==paste0('woody_savannas_', year)]<-'woody_savannas'
  names(covs_year)[names(covs_year)==paste0('savannas_', year)]<-'savannas'
  names(covs_year)[names(covs_year)==paste0('grasslands_', year)]<-'grasslands'
  names(covs_year)[names(covs_year)==paste0('permanent_wetlands_', year)]<-'permanent_wetlands'
  names(covs_year)[names(covs_year)==paste0('croplands_', year)]<-'croplands'
  names(covs_year)[names(covs_year)==paste0('cropland_natural_vegetation_mosaic_', year)]<-'cropland_natural_vegetation_mosaic'
  names(covs_year)[names(covs_year)==paste0('fascicularis_', year)]<-'fascicularis'
  names(covs_year)[names(covs_year)==paste0('nemestrina_', year)]<-'nemestrina'
  names(covs_year)[names(covs_year)==paste0('leucosphyrus_group_', year)]<-'leucosphyrus_group'
  
  # add nontemporal covariates
  covs_year <- addLayer(covs_year, covs_nontemporal)
  
  # extract data
  covs_year_extract <- extract(covs_year, dat[idx, c('Longitude', 'Latitude')])
  
  # check they're all there
  stopifnot(all(colnames(dat_covs) %in% colnames(covs_year_extract)))
  
  # match up the column names so they're in the right order
  match <- match(colnames(dat_covs), colnames(covs_year_extract))
  
  # extract covariates for all points
  dat_covs[idx, ] <- covs_year_extract[, match]
  
}                                                                               


# combine covariates with the other info
dat_all <- cbind(dat, dat_covs)

# add host_species column and make it numeric
dat_all$Host_species <- NA

dat_all$Host_species[dat_all$Host=='mosquito']<- 1
dat_all$Host_species[dat_all$Host=='monkey']<- 2
dat_all$Host_species[dat_all$Host=='human']<- 3

# let runBRT know that the host_species column is a discrete variable 
dat_all$Host_species <- factor(dat_all$Host_species)


# prepare dummy prediction rasters for host, vector and monkey host species
mosquito_ras <- raster('data/clean/raster/mbs_mask.grd')
monkey_ras <- mosquito_ras + 1
human_ras <- mosquito_ras + 2

# add each layer to covs_current to make a covariate set for each species for model prediction
mosquito <- addLayer(covs_current, mosquito_ras)
names(mosquito)[names(mosquito)=='layer'] <- 'Host_species'
monkey <- addLayer(covs_current, monkey_ras)
names(monkey)[names(monkey)=='layer'] <- 'Host_species'
human <- addLayer(covs_current, human_ras)
names(human)[names(human)=='layer'] <- 'Host_species'                                                                      

# check which rows in dat_covs contain an NA i.e. missing covariates (outside mask)
outside_idx <- attr(na.omit(dat_covs), 'na.action')

stopifnot(is.null(outside_idx))

# ~~~~~~~~~~~~~~~~~
# run bootstrapped BRT models

ncpu <- 50
nboot <- ncpu*10

# get random bootstraps of the data (minimum 10 pres/10 abs)
data_list <- replicate(nboot,
                       subsamplePolys(dat_all,
                                      minimum = c(10, 10),
                                      replace = TRUE),
                       simplify = FALSE)

data_list <- lapply(data_list,
                    balanceWeights2)

# initialize the cluster
sfInit(parallel = TRUE, cpus = ncpu)
sfLibrary(seegSDM)

model_list <- sfLapply(data_list,
                       runBRT,
                       gbm.x = 9:ncol(data_list[[1]]),
                       gbm.y = 1,
                       n.folds = 10,
                       gbm.coords = 4:5,
                       wt = 2)

# get cv statistics in parallel
stat_lis <- sfLapply(model_list, getStats)

# convert the stats list into a matrix using the do.call function
stats <- do.call("rbind", stat_lis)

# save them
write.csv(stats,
          paste0(outpath,
                 'stats.csv'))

# save the relative influence scores
relinf <- getRelInf(model_list)

write.csv(relinf,
          file = paste0(outpath,
                        'relative_influence.csv'))

# stop the cluster
sfStop()


# ~~~~~~~~~~~~~~
# make predictions

# loop through each set of covariates, making model predictions and plotting  
for (host in c('human','monkey', 'mosquito')){
  
  # get rasterbrick for these covs
  prediction_covs <- get(host)
  
  # get predictions in parallel
  preds_list <- sfLapply(model_list,
                         makePreds,
                         pred_covs = prediction_covs)
  
  # summarise all the ensembles
  preds <- stack(preds_list)
  
  # summarise the predictions 
  preds_sry <- combinePreds(preds, parallel = FALSE)
  
  names(preds_sry) <- c('mean',
                        'median',
                        'lowerCI',
                        'upperCI')
  
  # save the prediction summary
  writeRaster(preds_sry,
              file = paste0(outpath, 
                            'parasite_', 
                            host),
              format = 'GTiff',
              overwrite = TRUE)
  
  # plot the risk map
  png(paste0(outpath,
             host,
             '_prediction_mean.png'),
      width = 2000,
      height = 2000,
      pointsize = 30)
  
  par(oma = rep(0, 4),
      mar = c(0, 0, 0, 2))
  
  plot(preds_sry[[1]],
       axes = FALSE,
       box = FALSE)
  
  dev.off()
  
}


# make prediction to rest of SE Asia, using human host raster

# get rasterbrick for SE Asia to predict to
seasia_covs <- brick('data/clean/raster/SEAsia_covs.grd')

# drop correlated layers
seasia_covs <- dropLayer(seasia_covs, c('EVI_mean', 'EVI_SD', 'TCB_mean'))

# prepare dummy prediction raster for human host
seasia_extent <- raster('data/clean/raster/SEAsia_extent.grd')
seasia_human_ras <- seasia_extent + 3
names(seasia_human_ras)[names(seasia_human_ras)=='layer'] <- 'Host_species'

# add human dummy raster to pred_covs
seasia_covs <- addLayer(seasia_covs, seasia_human_ras)

# initialize the cluster
sfInit(parallel = TRUE, cpus = ncpu)
sfLibrary(seegSDM)

# get predictions in parallel
preds_list_seasia <- sfLapply(model_list,
                              makePreds,
                              pred_covs = seasia_covs)

# summarise all the ensembles
preds_seasia <- stack(preds_list_seasia)

sfStop()

# summarise the predictions in parallel
preds_sry_seasia <- combinePreds(preds_seasia, parallel = FALSE)

names(preds_sry_seasia) <- c('mean',
                             'median',
                             'lowerCI',
                             'upperCI')


# save the prediction summary
writeRaster(preds_sry_seasia,
            file = paste0(outpath, 
                          'SEAsia'),
            format = 'GTiff',
            overwrite = TRUE)

# plot the risk map
png(paste0(outpath,
           'SEAsia_prediction_mean_human.png'),
    width = 2000,
    height = 2000,
    pointsize = 30)

par(oma = rep(0, 4),
    mar = c(0, 0, 0, 2))

plot(preds_sry_seasia[[1]],
     axes = FALSE,
     box = FALSE)

dev.off()


# ~~~~~~~~~~~~~~~
# plot effect curves

# plot the conditional effect curves
#   png(paste0(outpath,
#              'human_conditional_effects.png'),
#       width = 2000,
#       height = 2500, 
#       pointsize = 30)
#   
#  # par(mfrow = n2mfrow(length(model_list[[1]]$effects)))
#  
#  effects <- getConditionalEffectPlots(model_list,
#                                       plot = TRUE, hold = 14, value = 3)
#  
#  dev.off()

# # plot the marginal effect curves
# png(paste0(outpath,
#            'effects.png'),
#     width = 2000,
#     height = 2500,
#     pointsize = 30)
# 
# par(mfrow = n2mfrow(length(model_list[[1]]$effects)))
# 
# effects <- getEffectPlots(model_list,
#                           plot = TRUE)
# 
# dev.off()

# plot marginal effect curves for covariates where RI > 100/no. of covs
RI_threshold <- 100/nlayers(covs_current)

# find number of covs with RI > than threshold and get names
key_covs <- which(relinf[1:nlayers(prediction_covs)] > RI_threshold)
names <- rownames(relinf[1:length(key_covs),])

# get index of key covs
key_covs_idx <- which(names(prediction_covs) %in% names)

# create raster stack of covs with RI > threshold
key_covs_ras <- prediction_covs[[key_covs_idx]]

# get the order of plots (all except relinf)!
order <- match(names, names(key_covs_ras))

# set up x axis labels and titles (in order they appear in prediction_covs, enter manually)
short_names <- c('Human population',
                 'Elevation',
                 'Urban accessibility',
                 'Croplands',
                 'Macaca nemestrina')

units <- c('people',
           'metres',
           'time',
           'proportion',
           'suitability index') # check correct axis labels!!!

effect <- getEffectPlots(model_list, plot=FALSE)

# subset effect to only include those from key covs
key_effect <- effect[key_covs_idx] 

# set up device
png(paste0(outpath,
           'effects_figure.png'),
    width = 3000,
    height = 3000,
    pointsize = 60)

# set up multi panels and margins
par(mfrow = c(3,3),
    mar = c(5, 2, 4, 2),
    oma = c(0, 3, 0, 0))

# loop through plots
for (i in 1:length(key_covs)) {
  
  # extract summary stats
  df <- key_effect[[order[i]]][, 1:4]
  
  # pick y axis
  if (i %% 3 == 1) {
    ylab = 'marginal effect'
  } else {
    ylab = ''
  }
  
  # set up empty plotting region
  plot(df[, 2] ~ df[, 1],
       type = 'n',
       ylim = c(min(c(df[, 3], rev(df[, 4]))), max(c(df[, 3], rev(df[, 4])))+0.5),
       ylab = '',
       xlab = '')
  
  # add the 95% CIs
  polygon(x = c(df[, 1], rev(df[, 1])),
          y = c(df[, 3], rev(df[, 4])),
          border = NA,
          col = grey(0.7))
  
  # add the mean line
  lines(df[, 2] ~ df[, 1],
        lwd = 5,
        col = grey(0.2))
  
  # y axis lable (only on left hand column)
  title(ylab = ylab,
        cex.lab = 1.2,
        col.lab = grey(0.3),
        xpd = NA,
        line = 2.5)
  
  # x-axis label
  title(xlab = units[order[i]],
        cex.lab = 1.2,
        col.lab = grey(0.3),
        line = 2.5)
  
  # title
  title(main = short_names[order[i]],
        line = 1.5,
        cex.main = 1.2)
  
  # relative contribution inset
  mtext(text = round(relinf[i, 1] / 100, 2),
        side = 3,
        line = -1.8,
        adj = 0.07,
        col = grey(0.5))
  
  # get x values for data distribution plot
  #x_vals <- unique(dat_all[, which(colnames(data_list[[1]])==names[i])])
  #y_vals <- rep(1.4, length(x_vals))
  
  #points(x_vals, y_vals, cex=0.4, bg='black', pch=21)
  
}

dev.off()

# generate mess map
# drop host species layer
seasia_covs_mss <- dropLayer(seasia_covs, 'Host_species')

# calculate mess map
mss <- mess(seasia_covs_mss, dat_covs, full=TRUE)

# save rmess layer and the raster stack
writeRaster(mss[['rmess']], 
            file='output/mess_raw', 
            format='GTiff', 
            overwrite=TRUE)

writeRaster(mss, 
            file='output/mess_maps',
            overwrite=TRUE)

# make binary, interpolation/extrapolation map and save
tmp <- mss[['rmess']] >= 0

tmp_masked <- mask(tmp, seasia_extent)

writeRaster(tmp_masked, 
            file='output/mess_binary', 
            format='GTiff', 
            overwrite=TRUE)


# calculating standard deviation value for each pixel
# predict, type='response', which returns probabilities for bernoulli input data

# load required package
library(logit)

# to get back to log odds scale need to use logit function
preds_list_seasia_tf <- lapply(preds_list_seasia, logit)

# summarise all the ensembles
preds_seasia_tf <- stack(preds_list_seasia_tf)

# function to calculate the mean and standard deviation for each pixel across bootstraps
combine2 <- function(x){
  
  ans <- c(mean=mean(x), 
           sd=sd(x))
  
  return(ans)
}

preds_sry_tf <- calc(preds_seasia_tf, fun=combine2)

# save standard deviation plot
writeRaster(preds_sry_tf[[2]],
            file='output/sd_raster.tif',
            format='GTiff',
            overwrite=TRUE)





