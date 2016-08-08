##### PLOT COORDINATES AND LINES TEST #####
library(sp)
data(meuse)
str(meuse) # what is the structure
coordinates(meuse) <- c("x", "y")
windows(); plot(meuse)
title("points")
cc <- coordinates(meuse)
m.sl <- SpatialLines(list(Lines(list(Line(cc)), "1")))
plot(m.sl)
title("lines")

##### CREATE A SPATIALPOINTS OBJECT #####
library(sp)
data(meuse)
coords <- SpatialPoints(meuse[, c("x", "y")])
summary(coords)

##### CREATE A SPATIALPOINTS DATAFRAME OBJECT #####
meuse1 <- SpatialPointsDataFrame(coords, meuse)
names(meuse1)

##### PLOT THE SPATIALPOINTS DATAFRAME OBJECT #####
plot(as(meuse1, "Spatial"), axes = TRUE)
plot(meuse1, add = TRUE)
plot(meuse1[meuse1$ffreq == 1, ], col = "green", add = TRUE)
meuse1$ffreq1 <- as.numeric(meuse1$ffreq)
plot(meuse1, col = meuse1$ffreq1, pch = 19)
labs <- c("annual", "every 2-5 years", "> 5 years")
cols <- 1:nlevels(meuse1$ffreq)
legend("topleft", legend = labs, col = cols, pch = 19, bty = "n")

###### CHOOSE CLASS INTERVALS AND COLOURS #######
library(classInt)
library(RColorBrewer)
q5 <- classIntervals(meuse1$zinc, n = 5, style = "quantile")
q5
pal <- brewer.pal(3, "Blues")
pal
plot(q5, pal = pal)
q5Colours <- findColours(q5, pal)
plot(meuse1, col = q5Colours, pch = 19)
legend("topleft", fill = attr(q5Colours, "palette"), legend = names(attr(q5Colours, "table")), bty = "n")

###### SSPLOT() ######
meuseZincSP <- spplot(meuse1, "zinc", at = q5, col.regions = pal, legendEntries = c("under 186", "186-246", "246-439", "439-737", "over 737"))
print(meuseZincSP)

##### BUBBLE PLOTS ######
library(lattice)
bubble1 <- bubble(meuse1, "zinc", maxsize = 2, key.entries = 100 * 2^(0:4))
print(bubble1)

#### CONTOUR LINES AND GRIDDED DATA ####
library(maptools)
volcano_sl <- ContourLines2SLDF(contourLines(volcano))
volcano_sl$level1 <- as.numeric(volcano_sl$level)
pal <- terrain.colors(nlevels(volcano_sl$level))
plot(volcano_sl, bg = "grey70", col = pal[volcano_sl$level1], lwd = 3)
data(meuse.grid)
coords <- SpatialPixels(SpatialPoints(meuse.grid[, c("x", "y")]))
meuseg1 <- SpatialPixelsDataFrame(coords, meuse.grid)
meuseg1$ffreq1 <- as.numeric(meuseg1$ffreq)
image(meuseg1, "ffreq1")
bpal <- colorRampPalette(pal)(41)
gridPlotSP <- spplot(meuseg1, "dist", col.regions = bpal, cuts = 40)
print(gridPlotSP)

###### Raster data ######
library(rgdal)
sp27gtiff <- readGDAL("/Users/charlie/Dropbox/asdarHowTo/asdarExamples/SP27GTIF.TIF")
sp27gtiff <- readGDAL("C:/Users/April/Downloads/hdimg_57540_0000_1770_svar_unkn_1_tf.TIF")

##### Using maptools to read data #####
library(maptools)
list.files()
scot<-readShapePoly("scot_BNG.shp")
