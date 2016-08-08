### Kriging

### Step 1: import data and shapefile for coordinates

library(maptools)

ehec <- read.csv("C:/Users/April/Desktop/MSc project/ireland complete 2014-04-24.csv")

library(maptools)

ehec.town <- readShapeSpatial("C:/Users/April/Desktop/MSc project/Data/Towns.shp") 
windows()
plot(ehec.town)
ehec.map <- readShapeSpatial("C:/Users/April/Desktop/Msc project/Data/city_towns.shp")
plot(ehec.map)

library(spdep)

ehec$rmort <- EBest(ehec$n_dead, ehec$n_pop)[,"raw"]        # raw and smoothed 
ehec$smort <- EBest(ehec$n_dead, ehec$n_pop)[,"estmm"]      # mortality 

ehec$east <- coordinates(ehec.map)[,1]
ehec$north <- coordinates(ehec.map)[,2]

library(geoR)
ehec.geo <- as.geodata(cbind(ehec$east, ehec$north, ehec$smort)) 
windows()
plot.geodata(ehec.geo)



####################### STEP 3: Kriging

### define grid for prediction locations = pixles
bound.shp <- readShapeSpatial("C:/Users/April/Desktop/Msc project/Data/city_towns.shp")

ehec.bdr <- read.table("C:/Users/April/Desktop/Stats for Health Sciences (POPM 6290)/Take home project/repom6290takehomeproject/EHEC-THP/EHEC-THP-2013/LSaxony.txt", header=T)
ehec.bdr

ehec.ply <- as.matrix(cbind(ehec$east, ehec$north))
krige.grid <- expand.grid(seq(min(ehec$east), max(ehec$east), l=100), # 100 lines
                          seq(min(ehec$north), max(ehec$north), l=100))

# all observations are weighted by their spatial dependence

windows() 
plot(ehec.map, type="l")              # plot ontario boundary
points(krige.grid, pch="+")          # plot grid location
points(ehec$east, ehec$north, pch=16, col="blue")   # plot data locations, i.e. centroids

### geostatistical kriging predictions at grid locations

ehec.uk <- krige.conv(ehec.geo, 
                     krige = krige.control(
                       type.krige = "ok",   		
                       trend.d = "cte", trend.l = "cte",  # "cte" means constant trend (ie. no trend)  
                       obj.model = ehec.beta), # model is here
                     locations = krige.grid, borders=ehec.ply) # use polygon to say just inside for prediction

### Mapping the predicted values as percentages

win.graph()
junk <- ehec.uk
junk$predict <- -junk$predict # the predicted values are the risks, so between 0 and 1, so if you are using heat colours, then have to switch them like this
junk$predict
image(junk, krige.grid,
      col = heat.colors(10), 
      xlab=c("easting"), ylab=c("northing"))
# plotting the negative predicted value to indicate high risk by red 
junk <- ehec.uk
junk$predict <- 100*junk$predict
contour(junk, labcex=1.2, method = "edge", ad=T)

lines(wnv.ply, col="white", lwd=5)    # make boundary nice
lines(wnv.ply, col="black", lwd=2)


### Last step: Relative Risk Map (risk of exposed vs. risk of unexposed)

# Idea scale the risk map by the baseline risk
# define baseline risk as the risk outside the disease cluster
# i.e. as the risk of the less exposed population

# Step 1: apply scan test to identify cluster
# Step 2: estimate baseline risk
# Step 3: scale risk map to Relative Risk Map

my.dat <- as.data.frame(coordinates(lsax.regions.shp))   # extract centroids
names(my.dat) <- c("x","y")

my.dat$n <- ehec$n98
my.dat$m <- ehec$m99

my.dat$raw <- 100*(ehec$m99 / ehec$n98)

library(SpatialEpi)
my.geo <- my.dat[c("x","y")]

my.expect <- expected(cases=my.dat$m, population=my.dat$n ,1)
my.cluster <- kulldorff(geo=my.geo, cases=my.dat$m, 
                         population=my.dat$n, expected.cases=my.expect,
                         pop.upper.bound=0.30, alpha.level=0.05, n.simulations=999, plot=F)

my.cluster$most.likely.cluster[c(1:5,8)]
my.cluster$secondary.clusters[[1]]$location.IDs.included

n.unexp <- (sum(ehec$n) # total number of birds, minus birds in one cluster, minus birds in the other cluster
            - my.cluster$most.likely.cluster$population 
            - my.cluster$secondary.clusters[[1]]$population)
m.unexp <- (sum(ehec$m) # total number of cases, minus cases in one cluster, minus cases in the other cluster
            - my.cluster$most.likely.cluster$number.of.cases 
            - my.cluster$secondary.clusters[[1]]$number.of.cases)
base.risk <- m.unexp / n.unexp # this is the base risk! this is the risk outside of the clusters

# plotting the relative risk = risk / base-risk 
windows()
junk <- ehec.uk # object with all these pixel values
junk$predict <- -junk$predict / base.risk
?image
image(junk, krige.grid,
      col = heat.colors(10), 
      xlab=c("easting"), ylab=c("northing"), main="Relative risk of EHEC in Lower Saxony")
junk <- ehec.uk
junk$predict <- junk$predict / base.risk
contour(junk, labcex=1.2, method = "edge", ad=T)

lines(ehec.ply, col="white", lwd=5)    # make boundary nice
lines(ehec.ply, col="black", lwd=2)

# draw a circle

my.geo[30,]
spDistsN1(as.matrix(my.geo[7,]), as.matrix(my.geo[7,]))

# Add circular scanning windows 
# plot a circle for each significant cluster
# need to know the center of the cluster and its diameter
# is provided by SaTScan
# will have to see how to extract this from R Kulldorff() results
# cluster center is the first region in the list
# diameter is distance to last region in cluster list

cluz <- seq(0, 2*pi, length=1000)
clux <- sin(cluz)
cluy <- cos(cluz)

polygon((clux*48067.85)+ 7.2226e+006,(cluy*48067.85)+ 936753, col=4, dens=0, lwd=5)
polygon((clux*63198.42)+ 7.02266e+006,(cluy*63198.42)+ 750660, col=4, dens=0, lwd=5)
