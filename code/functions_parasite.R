### function for sampling data containing polygon, sampling from the unique_ID column

subsamplePolys <- function (data, ...) {
  # given a presence-background dataset, with multiple rows for some of the
  # occurrence records, subset it so that there's only one randomly selected
  # point from each polygon and then take a bootstrap it using `subsample`.
  # Dots argument is passed to subsample.  
  
  # subset to get polygon data only, in the parasite dataset this might include 
  # presence/absence point data (if the 6km jitter around infection sites is used)
  # but not background data
  #poly_dat <- subset(data, data$true==1) - if jitter is used
  
  poly_dat <- subset(data, data$Geometry_type!='point')
  
  # get the different IDs
  u <- unique(poly_dat$Unique_ID)
  
  # remove the NA values
  u <- u[!is.na(u)]
  
  # loop through, picking an index for each based on the number available
  data_idx <- sapply(u,
                     function (identifier, data) {
                       idx <- which(data$Unique_ID == identifier)
                       sample(idx, 1)
                     },
                     data)
  
  # get the subsetted dataset
  dat <- data[data_idx, ]
  
  # get the point data, in this case the background data points
  dat_point <- subset(data, data$Geometry_type=='point')
  
  # append the subsetted polygon data to the point data
  dat_all <- rbind(dat, dat_point)
  
  # randomly subsample the dataset
  ans <- subsample(dat_all,
                   n = nrow(dat_all),
                   ...)
  
  return (ans)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~

insideMask <- function(points, mask){
  
  # checks whether lats/longs fall on non-missing pixels of a raster 
  # takes two arguments: points, a dataframe containing columns named 'Longitude' and 'Latitude' 
  # mask is a raster 
  # returns a dataframe containing only those rows with points falling on non-missing pixels
  # if all points fall on missing pixels, the function throws an error
  
  # check whether arguments are in the correct format
  stopifnot(inherits(mask, 'Raster'))
  stopifnot(inherits(points, 'data.frame'))
  
  # get indexes of points which fall inside the mask
  inside_idx <- which(!(is.na(extract(mask, points[,c('Longitude', 'Latitude')]))))
  
  # subset these points
  inside_points <- points[inside_idx, ]
  
  # if raster values exist for one point or more, return the point/s, otherwise throw an error 
  if (nrow(inside_points)==0) {stop('all points outside mask')}
  
  return (inside_points)
   
}

balanceWeights <- function(data) {
  # given a mixed spatial data frame of presence, true absence and pseudo-absence data, ensure the 
  # sum of the weights of the presences and true absences plus pseudo-absences are balanced
  
  # subset 'true' presences and absences
  true_dat <- data[data$true==1,]
  
  # subset pseudo-absences
  pseudo <- data[data$true==0,]
  
  # loop through hosts caculating background weights 
  for (host in c('human', 'mosquito', 'monkey')) {
    
    # get true data for this host 
    host_dat <- true_dat[true_dat$Host==host,]
    
    presence_total <- sum(host_dat[host_dat$PA==1,]$wt)
    absence_total <- sum(host_dat[host_dat$PA==0,]$wt) + sum(pseudo[pseudo$Host==host,]$wt)
    
    pseudo[pseudo$Host==host,]$wt <- pseudo[pseudo$Host==host,]$wt * (presence_total/absence_total)
    
    # calculate pseudo-absence weight for this host 
    
    #weight <- (nrow(host_dat[host_dat$PA==1,]) - nrow(host_dat[host_dat$PA==0,])) / nrow(pseudo[pseudo$Host==host,])
    
    # put weight into pseudo-absences data subset  
    #pseudo[pseudo$Host==host,]$wt <- weight
    
  }
  
  # re-combine 'true' presences and absences and pseudo-absences
  result <- rbind(true_dat, pseudo)
  
  return(result)
  
}

balanceWeights2 <- function(data) {
  # given a mixed spatial data frame of presence and pseudo-absence data, ensure the 
  # sum of the weights of the presences and pseudo-absences are balanced
  
  # subset presence data to add to
  dat_presence <- data[data$PA==1,]
  
  # loop through hosts caculating background weights 
  for (host in c('human', 'mosquito', 'monkey')) {
    
    # subset data for this host 
    host_dat <- data[data$Host==host,]
    
    presence <- host_dat$PA == 1
    absence <- host_dat$PA == 0
    presence_total <- sum(host_dat[presence, ]$wt)
    absence_total <- sum(host_dat[absence, ]$wt)
    
    # update weights in host data
    host_dat[absence, 'wt'] <- host_dat[absence, ]$wt * (presence_total / absence_total)
    
    # updated host absence data
    new_absence <- host_dat[absence,]
    
    # bind update host absence data with presence data
    dat_presence <- rbind(dat_presence, new_absence)
    
    
  }
  
  return(dat_presence)
  
}
