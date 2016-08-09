# Code to analyze the disease clustering of cholera in Ireland during the 
#  cholera outbreak of 1848 - 1850. Disease clustering is measured using 
#  Moran's I correlation coefficients which test the strength of the spatial 
#  clustering of disease. A correlogram of Moran's I coefficients visualizes 
#  how far cholera observations are dependent or clustering.

#### Initial configuration ####
# R version 3.1.3 (2015-03-09)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 8 x64 (build 9200)

# Create a new comment

#### Working directory ####
# Clone or download the repository and set the working directory
#  with setwd to the folder where the repository is located.

setwd("C:/Users/April/Desktop/2014_MSc_Cholera-Ireland-1848-1850/data")

#### Libraries ####

library(maptools) # version 0.8-34
library(spdep) # version 0.5-88

#### Constants ####

data.file.path <- "Analysis-data_Counties_Cholera-Ireland-1848-1850.txt"
shapefile.path <- "Ireland-shapefile_County"

##### Analysis of disease clustering with regional data #####

#### Create neighbourhood structures ####

county.map <- readShapeSpatial(fn=shapefile.path)

# queen = FALSE selects 1st order neighbourhood structure by county borderline
county.neighbours.borderline <- poly2nb(county.map, queen = FALSE)

# shows the neighbourhood structure of the object
str(county.neighbours.borderline)

# creates a matrix object of the county centroids
coords <- as.data.frame(coordinates(county.map))
names(coords) <- c("x", "y")

# county centroid data added to county cholera dataset
county.cholera.data <- read.table(data.file.path, header=T)
county.cholera.data$x <- coords$x
county.cholera.data$y <- coords$y

# county incidence data added to county cholera dataset
county.cholera.data$inc <- county.cholera.data$CASES / county.cholera.data$POP

#### Moran's I correlation coefficients ####
# tests the strength of the spatial clustering of disease with regional data

moran.test(
  county.cholera.data$inc, 
  nb2listw(county.neighbours.borderline), 
  alternative = "two.sided"
)

##### Correlogram of Moran's I coefficients
# The correlogram visualizes how far observations are dependent or clustering

morans.i.corr <- sp.correlogram(
  county.neighbours.borderline, 
  county.cholera.data$inc, 
  order = 6, 
  method = "I", 
  zero.policy = TRUE
)
print(morans.i.corr)
plot(morans.i.corr)
