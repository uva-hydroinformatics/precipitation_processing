library(RGeostats)

dates=c('2013-07-02', '2013-10-09', '2014-01-11', '2014-02-13', '2014-04-15', '2014-04-25', '2014-07-10', '2014-08-18', '2014-09-08', '2014-09-09', '2014-09-13', '2014-11-26', '2014-12-24', '2015-04-14', '2015-06-02', '2015-06-24', '2015-08-07', '2015-08-20', '2015-09-30', '2015-10-02')

#read in csvfile and change into rgeostats db object
rain.csv <- read.csv(file="C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\precipitation_processing\\combined_aggregate.csv", head=TRUE, sep=",")
rain.db <- db.create(rain.csv,ndim=2,autoname=F,flag.grid=F)

#set precip to the z variable
rain.db <- db.locate(rain.db, "sitename")
rain.db <- db.locate(rain.db, "precip_mm","z")
orig.db <- rain.db

#create neighbor objects for interpolation
unique.neigh <- neigh.init(ndim = 2, type = 0)
moving.neigh <- neigh.init(ndim = 2, type = 2, nmini = 1, nmaxi = 8, nsect = 8, nsmax = 2, flag.sector=TRUE, dmax = 10000)

outliers <- matrix(c("Index", "High quantile", "low quantile", "interquantile range"), ncol = 4, byrow = 1)

for (i in 1:20){
    print(i)  
    rain.db <- orig.db
    
    #select only data for given date
    rain.db <- db.sel(rain.db,datetime==dates[i])
    #change which sources you would like to consider
    type <- "wu_vab"
    rain.db <- db.sel(rain.db,src=="wu" | src=="vab", combine="and")
    
    #plot boxplots 
    filename <- paste(dates[i], "boxplot_", type, sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\R\\figures\\boxplots\\'
    jpeg(paste(path,filename,fileext,sep=""))
    rain_vtr = db.extract(rain.db, "precip_mm")
    boxplot(rain_vtr)
    dev.off()
   
    #check what the outliers are
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

    #create semivariogram including 4 directional one to check for anisotropy
    data.vario <- vario.calc(rain.db,lag=500,nlag=40)
    plot(data.vario,npairdw=TRUE,npairpt=1)
    data.4dir.vario <- vario.calc(rain.db,lag=1000,nlag=20,dir=c(0,45,90,135))
    plot(data.4dir.vario,title="Directional variograms")
    data.model <- model.auto(data.vario,struct=c("Spherical","Exponential"),title="Modelling omni-directional variogram")
    
    #save experimental variogram and variogram model as jpeg
    filename <- paste(dates[i], "variogram_", type, sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\R\\figures\\variograms\\'
    jpeg(paste(path,filename,fileext,sep=""))
    plot(data.vario,npairdw=TRUE,npairpt=1, title=paste(type, "variogram for", dates[i], sep=" "), xlab = "lag (m)", ylab = "Variance")
    plot(data.model, add=TRUE)
    dev.off()
    
    #create unique and moving neighborhood objects and validate with moving neighborhood
    rain.db <- xvalid(rain.db, data.model, moving.neigh)
    hist(db.extract(rain.db,"Xvalid.precip_mm.esterr"), nclass=30, main="Histogram of precip error" , xlab="Cross-validation error", col="blue")
    
    #create grid and perform kriging
    grid.db <- db.grid.init(rain.db,nodes=c(100,100))
    rain.db <- db.locate(rain.db,seq(11,12))
    rain.db <- db.locate(rain.db,"precip_mm","z")
    grid.db <- kriging(rain.db,grid.db,data.model,unique.neigh,radix="ku")
    
    #plot kriging interpolation
    ramp <- colorRampPalette(c("white","blue"))
    color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
        scale = (length(lut)-1)/(max-min)
        #dev.new(width=1.75, height=5)
        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='Precipitation_mm', main=title)
        axis(2, ticks, las=1)
        for (i in 1:(length(lut)-1)) {
         y = (i-1)/scale + min
         rect(0,y,10,y+1/scale, col=lut[i], border=NA)
        }
    }
    filename <- paste(type, 'kriging', dates[i], sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\R\\figures\\kriging\\'
    jpeg(paste(path,filename,fileext,sep=""))
    plot(grid.db, name.image="ku.precip_mm.estim", col=ramp(100), title="Estimation of precipitation for 20130703")
    plot(rain.db, pch=21, bg=1, add=T)
    par(fig=c(0.8,1,0,1),new=T)
    rain.est <- db.extract(grid.db, "ku.precip_mm.estim")
    color.bar(ramp(50),min(rain.est),max(rain.est))
    dev.off()
    
    #clean up
    rain.db <- db.delete(rain.db, seq(8,11))
}



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
write.table(outlier_table, file="outliers.csv", col.names = NA, sep = ",")



