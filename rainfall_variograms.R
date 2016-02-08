library(RGeostats)

get.sum.stats <- function(data, date, lag, nlag, outliers){
  data.vario <- vario.calc(data,lag=lag,nlag=20)
  data.model <- model.auto(data.vario,struct=c("Spherical","Exponential"),title="Modelling omni-directional variogram", draw = F)
  model.y <- data.vario$vardirs[[1]]$gg
  single.sum.stats <- data.frame("date" = dates[i],
                                 "lag" = lag, 
                                 "sill" = as.numeric(data.model$basics[[1]]$sill), 
                                 "range" = as.numeric(data.model$basics[[1]]$range),
                                 "modeltype" = as.character(data.model$basics[[1]]$vartype),
                                 "variance" = var(model.y, na.rm = T),
                                 "mean" = mean(model.y, na.rm = T),
                                 "numoutliers_excluded" = length(outliers)
  )
  return(single.sum.stats)
  
}




#read in csvfile and change into rgeostats db object
rain.csv <- read.csv(file="C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\precipitation_processing\\combined_aggregate_15_min.csv", head=TRUE, sep=",")
rain.db <- db.create(rain.csv,ndim=2,autoname=F,flag.grid=F)

dates=levels(unique(rain.csv$datetime))


#set precip to the z variable
rain.db <- db.locate(rain.db, "sitename")
rain.db <- db.locate(rain.db, "precip_mm","z")
#change which sources you would like to consider
type <- "wu_vab_filt_15_min"
rain.db <- db.sel(rain.db,src=="wu" | src=="vab")

#lag and number of lags for variograms
lag <- 500
nlag <- 20


orig.db <- rain.db
sum.stats <- data.frame()

for (i in 1:length(dates)){
    print(i)  
    rain.db <- orig.db
    
    #select only data for given date
    rain.db <- db.sel(rain.db,datetime==dates[i], combine="and")
    
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
    
    #add single with outliers stats to sum.stats dataframe
    ss = get.sum.stats(rain.db, dates[i], lag, nlag, c())
    sum.stats <- rbind(sum.stats, ss)
                                                        
    
    #exclude outliers
    rain.db <- db.sel(rain.db, (precip_mm < u.thresh & precip_mm > l.thresh), combine="and")
    
    #add single without outliers stats to sum.stats dataframe
    if(length(outlier_ranks)>0){
      ss <- get.sum.stats(rain.db, dates[i], lag, nlag, outlier_ranks)
      sum.stats <- rbind(sum.stats, ss)
    }

    #create experimental semivariogram and model
    data.vario <- vario.calc(rain.db,lag=lag,nlag=nlag)
    data.model <- model.auto(data.vario,struct=c("Spherical","Exponential"),title="Modelling omni-directional variogram", draw = F)

    #save experimental variogram and variogram model as jpeg
    filename <- paste(dates[i], "variogram_", type, sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\R\\figures\\variograms\\'
    jpeg(paste(path,filename,fileext,sep=""))
    plot(data.vario,npairdw=TRUE,npairpt=1, title=paste(type, "variogram for", dates[i], sep=" "), xlab = "lag (m)", ylab = "Variance")
    plot(data.model, add=TRUE)
    dev.off()
  
    #clean up
    rain.db <- db.delete(rain.db, seq(8,11))
}

write.table(sum.stats, file="variogram_summary_stats_15_min.csv", sep = ",", row.names = F)


