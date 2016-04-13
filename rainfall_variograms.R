library(RGeostats)

#read in csvfile and change into rgeostats db object
filename = "daily_tots.csv"
rain.csv <- read.csv(file=paste("C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\data\\", filename, sep=""), head=TRUE, sep=",")
non_zero_columns = which(colSums(rain.csv[-1:-4], na.rm=T) !=0) + 4
rain.csv = rain.csv[, c(1:4, non_zero_columns)]
rain.filt = rain.csv[,c(2,3,non_zero_columns)]
rain.db <- db.create(rain.filt,ndim=2,autoname=F,flag.grid=F)

#change which sources you would like to consider
type <- "all data"
#rain.db <- db.sel(rain.db,src=="hrsd" | src=="vab")

#lag and number of lags for variograms
lag <- 1000
nlag <- 20


orig.db <- rain.db

unique.neigh <- neigh.init(ndim = 2, type = 0)
moving.neigh <- neigh.init(ndim = 2, type = 2, nmini = 1, nmaxi = 8, nsect = 8, nsmax = 2, flag.sector=TRUE, dmax = 10000)


for (i in 4:(length(rain.db@items))){
    rain.db <- orig.db
    
    #select only data for given date
    rain.db <- db.locate(rain.db, i, "z")
    
    #create experimental semivariogram and model
    data.vario <- vario.calc(rain.db,lag=lag,nlag=nlag)
    data.model <- model.auto(data.vario,title="Modelling omni-directional variogram", draw = F)

    #save experimental variogram and variogram model as jpeg
    filename <- paste(i, "variogram_", type, sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\figures\\R\\figures\\variograms\\rel_hrsd_all\\'
    #jpeg(paste(path,filename,fileext,sep=""))
    col_name = colnames(rain.filt)[i-1]
    date = strsplit(col_name, "X")[[1]][2]
    plot(data.vario,npairdw=TRUE,npairpt=0, title=paste("", "variogram for", date, sep=" "), xlab = "lag (m)", ylab = "Variance (mm^2/d)", pin = c(1,2), varline=F)
    plot(data.model, add=TRUE, pin=c(1,1))
    #dev.off()
    
    # Kriging ####################
    rain.db <- xvalid(rain.db, data.model, moving.neigh)
    hist(db.extract(rain.db, paste("Xvalid.", col_name, ".esterr", sep="")), nclass=30, main="Histogram of precip error" , xlab="Cross-validation error", col="blue")
    
    #create grid and perform kriging
    grid.db <- db.grid.init(rain.db,nodes=c(100,100))
    rain.db <- db.locate(rain.db,seq(24,25))
    rain.db <- db.locate(rain.db, i, "z")
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
    filename <- paste(type, 'kriging', date, sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\R\\figures\\kriging\\'
    #jpeg(paste(path,filename,fileext,sep=""))
    rain.est <- db.extract(grid.db, paste("ku.", col_name, ".estim", sep=""))
    plot(grid.db, name.image=paste("ku.", col_name, ".estim", sep=""), col=ramp(100), title="Estimation of precipitation for 20130703")
    legend.image(zlim=c(min(rain.est),max(rain.est)), col=ramp(100), position="topright")
    plot(rain.db, pch=21, bg=1, add=T)

    
      
    #clean up
    rain.db <- db.delete(rain.db, seq(24,25))
}



