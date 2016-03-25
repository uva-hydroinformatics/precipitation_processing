library(RGeostats)

#read in csvfile and change into rgeostats db object
filename = "2014-07-10-fifteen_min.csv"
rain.csv <- read.csv(file=paste("C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\data\\all data\\", filename, sep=""), head=TRUE, sep=",")
non_zero_columns = which(colSums(rain.csv[-1:-4], na.rm=T) !=0)
rain.csv = rain.csv[, c(1:4, non_zero_columns+4)]
rain.db <- db.create(rain.csv,ndim=2,autoname=F,flag.grid=F)


#set precip to the z variable
rain.db <- db.locate(rain.db, "sitename")
rain.db <- db.locate(rain.db, "x", "x")
rain.db <- db.locate(rain.db, "y", "x")

#change which sources you would like to consider
type <- "all data"
#rain.db <- db.sel(rain.db,src=="hrsd" | src=="vab")

#lag and number of lags for variograms
lag <- 1000
nlag <- 20


orig.db <- rain.db

for (i in 5:(length(rain.db@items)-1)){
    rain.db <- orig.db
    
    #select only data for given date
    rain.db <- db.locate(rain.db, i, "z")
    
    #create experimental semivariogram and model
    data.vario <- vario.calc(rain.db,lag=lag,nlag=nlag)
    data.model <- model.auto(data.vario,struct=c("Spherical","Exponential"),title="Modelling omni-directional variogram", draw = F)

    #save experimental variogram and variogram model as jpeg
    filename <- paste(i, "variogram_", type, sep="_")
    fileext <- ".jpg"
    path <- 'C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\figures\\R\\figures\\variograms\\rel_hrsd_all\\'
    #jpeg(paste(path,filename,fileext,sep=""))
    date = colnames(rain.csv)[i]
    date = strsplit(date, "X")[[1]][2]
    plot(data.vario,npairdw=TRUE,npairpt=0, title=paste("", "variogram for", date, sep=" "), xlab = "lag (m)", ylab = "Variance (mm^2/d)", pin = c(1,2), varline=F)
    plot(data.model, add=TRUE, pin=c(1,1))
    #dev.off()
  
    #clean up
    rain.db <- db.delete(rain.db, seq(8,11))
}



