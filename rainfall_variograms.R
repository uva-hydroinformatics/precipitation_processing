library(RGeostats)

dates=c('2013-07-02', '2013-10-09', '2014-01-11', '2014-02-13', '2014-04-15', '2014-04-25', '2014-07-10', '2014-08-18', '2014-09-08', '2014-09-09', '2014-09-13', '2014-11-26', '2014-12-24', '2015-04-14', '2015-06-02', '2015-06-24', '2015-08-07', '2015-08-20', '2015-09-30', '2015-10-02')

#read in csvfile and change into rgeostats db object
rain.csv <- read.csv(file="C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\precipitation_processing\\combined_aggregate.csv", head=TRUE, sep=",")
rain.db <- db.create(rain.csv,ndim=2,autoname=F,flag.grid=F)

#set precip to the z variable
rain.db <- db.locate(rain.db, "sitename")
rain.db <- db.locate(rain.db, "precip_mm","z")
#change which sources you would like to consider
type <- "wu_vab_filt"
rain.db <- db.sel(rain.db,src=="wu" | src=="vab")

orig.db <- rain.db

#create neighbor objects for interpolation
unique.neigh <- neigh.init(ndim = 2, type = 0)
moving.neigh <- neigh.init(ndim = 2, type = 2, nmini = 1, nmaxi = 8, nsect = 8, nsmax = 2, flag.sector=TRUE, dmax = 10000)

outliers <- matrix(c("Index", "High quantile", "low quantile", "interquantile range"), ncol = 4, byrow = 1)

for (i in 1:length(dates)){
    print(i)  
    rain.db <- orig.db
    
    #select only data for given date
    rain.db <- db.sel(rain.db,datetime==dates[i], combine="and")
    
    #check what the outliers are
    rain_vtr = db.extract(rain.db, "precip_mm")
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

    #create semivariogram
    data.vario <- vario.calc(rain.db,lag=500,nlag=40)
    data.model <- model.auto(data.vario,struct=c("Spherical","Exponential"),title="Modelling omni-directional variogram")

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





