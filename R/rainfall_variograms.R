library(RGeostats)
library(XLConnect)
library(DBI)
library(RSQLite)

#read in csvfile and change into rgeostats db object
data_dir = "../../Data/"
con = dbConnect(RSQLite::SQLite(), dbname=paste(data_dir,"master.sqlite", sep=""))
table = "fif"
rain.csv <- dbGetQuery(con, paste("select * from ", table, sep=""))
non_zero_columns = which(colSums(rain.csv[-1:-4], na.rm=T) !=0) + 4

#change which sources you would like to consider
type <- "all data"
#rain.db <- db.sel(rain.db,src=="hrsd" | src=="vab")

#lag and number of lags for variograms
lag <- 1000
nlag <- 20

unique.neigh <- neigh.init(ndim = 2, type = 0)
moving.neigh <- neigh.init(ndim = 2, type = 2, nmini = 1, nmaxi = 8, nsect = 8, nsmax = 2, flag.sector=TRUE, dmax = 10000)

model.attr = data.frame()
l = length(non_zero_columns)
for (i in 1:l){
    print(paste("doing ", i, " of ", l, sep=""))
    
    rain.filt = rain.csv[,c(2,3,non_zero_columns[i])]
    rain.db <- db.create(rain.filt, ndim=2,autoname=F,flag.grid=F)
  
    
    date = colnames(rain.filt)[3]

    rain.db <- db.locate(rain.db, col_name, "z")
    
    #create experimental semivariogram and model
    data.vario <- vario.calc(rain.db,lag=lag,nlag=nlag)
    data.model <- model.auto(data.vario, struct = c("Spherical"), title="Modelling omni-directional variogram", draw = F)
    
    # model characteristics
    sill = as.numeric(data.model$basics[[1]]$sill)
    range = as.numeric(data.model$basics[[1]]$range)
    modeltype = as.character(data.model$basics[[1]]$vartype)
    m = data.frame('date' = date, 'sill' = sill, 'range' = range, 'type' = modeltype)
    model.attr = rbind(model.attr, m)

    #jpeg(paste(path,filename,fileext,sep=""))
#     plot.new()
#     par(pin=c(.1,.2))
#     split.screen(c(1, 2), erase = FALSE)
#     screen(1)
#     plot(data.vario,npairdw=TRUE,npairpt=0, title=paste("", "variogram for", date, sep=" "), xlab = "lag (m)", ylab = "Variance (mm^2/d)", pin = c(1,2), varline=F, reset=FALSE)
#     plot(data.model, add=TRUE, pin=c(1,1))
    #dev.off()
    
    # Kriging ####################
#     rain.db <- xvalid(rain.db, data.model, moving.neigh)
    #hist(db.extract(rain.db, paste("Xvalid.", col_name, ".esterr", sep="")), nclass=30, main="Histogram of precip error" , xlab="Cross-validation error", col="blue")
    
    #create grid and perform kriging
#     grid.db <- db.grid.init(rain.db,nodes=c(100,100))
#     rain.db <- db.locate(rain.db,seq(l+1, l+2))
#     rain.db <- db.locate(rain.db, i, "z")
#     grid.db <- kriging(rain.db,grid.db,data.model,unique.neigh,radix="ku")
    
    #plot kriging interpolation
#     screen(2)
#     ramp <- colorRampPalette(c("white","blue"))
#     filename <- paste(type, 'kriging', date, sep="_")
#     fileext <- ".jpg"
#     #jpeg(paste(path,filename,fileext,sep=""))
#     rain.est <- db.extract(grid.db, paste("ku.", col_name, ".estim", sep=""))
#     plot(grid.db, name.image=paste("ku.", col_name, ".estim", sep=""), col=ramp(100), title=paste("Estimation of precipitation for ", date, sep=""), reset=FALSE)
#     legend.image(zlim=c(min(rain.est),max(rain.est)), col=ramp(100), position="topright")
#     plot(rain.db, pch=21, bg=1, add=T)
    #dev.off()
    
}
dbWriteTable(con, paste(table, "_model_params", sep=""), model.attr, overwrite=TRUE)



