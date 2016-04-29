# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

#read in csvfile and change into rgeostats db object
filename = "combined_aggregate_15_min.csv"
rain.csv <- read.csv(file=paste("C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\precipitation_processing\\", filename, sep=""), head=TRUE, sep=",")

#subset hrsd for only the study area
hrsd_stations_in_study_area = c("MMPS-171","MMPS-185","MMPS-163","MMPS-255","MMPS-146","MMPS-004","MMPS-256","MMPS-140","MMPS-160","MMPS-144", "MMPS-036","MMPS-093-2")
levels(rain.csv$site_name) = trim.trailing(levels(rain.csv$site_name))
rain.sub.hrsd <- rain.csv[rain.csv$site_name %in% hrsd_stations_in_study_area,]
rain.csv <- subset(rain.csv, src == 'wu' | src == 'vab')
rain.csv <- rbind(rain.csv, rain.sub.hrsd)

dates=c('2013-07-02', '2013-10-09', '2014-01-11', '2014-02-13', '2014-04-15', '2014-04-25', '2014-07-10', '2014-08-18', '2014-09-08', '2014-09-09', '2014-09-13', '2014-11-26', '2014-12-24', '2015-04-14', '2015-06-02', '2015-06-24', '2015-08-07', '2015-08-20', '2015-09-30', '2015-10-02')

outliers <- data.frame()

# for (i in 1:length(dates)){
for (i in 1:length(dates[1])){
    print(paste(i, dates[i], sep=": "))  

    #select only data for given date
    r.date <- rain.csv[as.Date(rain.csv$date) %in% as.Date(dates[i]), c(1,2,6)]
    
    #check what the outliers are
    rain_vtr <- db.extract(rain.db, "precip_mm")
    upperq <- quantile(rain_vtr, na.rm = TRUE)[4]
    lowerq <- quantile(rain_vtr, na.rm = TRUE)[2]
    upperq <- unname(upperq)
    lowerq <- unname (lowerq)
    inter_quantile_range  <- upperq - lowerq
    range <- 1.5
    u.thresh <- upperq + (inter_quantile_range*range)
    l.thresh <- lowerq - (inter_quantile_range*range)
    l.thresh <- ifelse(l.thresh<0, 0, l.thresh)
    outlier_db <- db.sel(rain.db, precip_mm > u.thresh | precip_mm < l.thresh, combine="and")
    outlier_ranks = db.extract(outlier_db, c("rank"))
    if(length(outlier_ranks>0)){
      for (j in 1:length(outlier_ranks)){
        outliers <- rbind(outliers, c(outlier_ranks[j], upperq, lowerq, inter_quantile_range))  
      }
    }
    
    

    #exclude outliers
    rain.db <- db.sel(rain.db, (precip_mm < u.thresh & precip_mm > l.thresh), combine="and")

}

