library(RGeostats)
library(SpatialTools)


#read in csvfile and change into rgeostats db object
filename = "combined_aggregate_filt.csv"
rain.csv = read.csv(file=paste("C:\\Users\\jeff_dsktp\\Box Sync\\Sadler_1stPaper\\rainfall\\data\\", filename, sep=""), head=TRUE, sep=",")

#with wu points
pts = c(rain.csv$x, rain.csv$y)
pts_m = unique(matrix(pts, nrow=length(pts)/2, ncol=2, byrow = F))
d = dist1(pts_m)
m = mean(d)

#plot
plot(pts_m)


#without wu points
rain.csv1 = subset(rain.csv, src != "wu")
pts1 = c(rain.csv1$x, rain.csv1$y)
pts_m1 = unique(matrix(pts1, nrow=length(pts)/2, ncol=2, byrow = F))
d1 = dist1(pts_m1)
m1 = mean(d1)

#plot
plot(pts_m1)

print(paste("percent reduction in average distance between gauges: ", (m1-m)/m1*100, "%"))