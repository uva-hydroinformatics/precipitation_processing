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
filename = "combined_aggregate_filt.csv"
rain.csv <- read.csv(file=paste("C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\data\\", filename, sep=""), head=TRUE, sep=",")

rain.db <- db.create(rain.csv,ndim=2,autoname=F,flag.grid=F)

dates=levels(unique(rain.csv$datetime))


#set precip to the z variable
rain.db <- db.locate(rain.db, "sitename")
rain.db <- db.locate(rain.db, "precip_mm","z")

#change which sources you would like to consider
type <- "relevant_hrsd_filt_2000_13"
#rain.db <- db.sel(rain.db,src=="hrsd" | src=="vab")

#lag and number of lags for variograms
lag <- 1000
nlag <- 20


orig.db <- rain.db
sum.stats <- data.frame()

outliers <- matrix(c("Index", "High quantile", "low quantile", "interquantile range"), ncol = 4, byrow = 1)


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
    if(length(outlier_ranks>0)){
      for (j in 1:length(outlier_ranks)){
        outliers <- rbind(outliers, c(outlier_ranks[j], upperq, lowerq, inter_quantile_range))  
      }
    }
    
    
    #add single with outliers stats to sum.stats dataframe
    ss = get.sum.stats(rain.db, dates[i], lag, nlag, c())
    sum.stats <- rbind(sum.stats, ss)
                                                        
    
    #exclude outliers
    #rain.db <- db.sel(rain.db, (precip_mm < u.thresh & precip_mm > l.thresh), combine="and")
    
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
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\figures\\R\\figures\\variograms\\rel_hrsd_all\\'
    #jpeg(paste(path,filename,fileext,sep=""))
    plot(data.vario,npairdw=TRUE,npairpt=0, title=paste("", "variogram for", dates[i], sep=" "), xlab = "lag (m)", ylab = "Variance (mm^2/d)", pin = c(1,2), varline=F)
    plot(data.model, add=TRUE, pin=c(1,1))
    #dev.off()
  
    #clean up
    rain.db <- db.delete(rain.db, seq(8,11))
}

write.table(sum.stats, file=paste(type, "variogram_summary_stats.csv", sep="_"), sep = ",", row.names = F)

#compiles all of the outliers computes the number of iqrs from the quantile and writes them to a csv file
outlier_table <- matrix(c("site_name", "src", "datetime", "precip", "low/high", "quantile", "number of iqr's away"), ncol = 7, byrow = 1)
for (j in 2:nrow(outliers)){
  out_df = rain.db[as.numeric(outliers[j,1])]
  out_df <- as.vector(t(out_df))
  x <- as.numeric(out_df[7])
  site_name <- out_df[4]
  src <- out_df[5]
  datetime <- out_df[6]
  upperq = as.numeric(outliers[j,2])
  lowerq = as.numeric(outliers[j,3])
  inter_quantile_range = as.numeric(outliers[j,4])
  if(x>upperq){
    out_of_range <- (x-upperq)/inter_quantile_range
    outlier_table <- rbind(outlier_table,c(site_name, src, datetime, x, "high", upperq, out_of_range))
  }
  else{
    out_of_range <- (lowerq-x)/inter_quantile_range
    outlier_table <- rbind(outlier_table, c(site_name, src, datetime, x, "low", lowerq, out_of_range))
  }
}
write.table(outlier_table, file=paste(type, "outliers.csv", sep = "_"), col.names = NA, sep = ",")


