setwd("C:/Users/April/Desktop/MSc project/github/spacetime-vis")

library(lattice)
library(ggplot2)
library(latticeExtra)
library(gridBase)
library(gridExtra)


myTheme <- custom.theme.2(pch=19, cex=0.7,
                          region=rev(brewer.pal(9, 'YlOrRd')),
                          symbol = brewer.pal(n=8, name = "Dark2"))
myTheme$strip.background$col='transparent'
myTheme$strip.shingle$col='transparent'
myTheme$strip.border$col='transparent'

xscale.components.custom <- function(...){
  ans <- xscale.components.default(...)
  ans$top=FALSE
  ans}
yscale.components.custom <- function(...){
  ans <- yscale.components.default(...)
  ans$right=FALSE
  ans}
myArgs <- list(as.table=TRUE,
               between=list(x=0.5, y=0.2),
               xscale.components = xscale.components.custom,
               yscale.components = yscale.components.custom)
defaultArgs <- lattice.options()$default.args

lattice.options(default.theme = myTheme,
                default.args = modifyList(defaultArgs, myArgs))

##################################################################
## Choropleth maps
##################################################################


##################################################################
## Administrative boundaries
##################################################################

library(sp)
library(maptools)

#old <- setwd(tempdir())

IreMap <- readShapePoly(fn="Ireland choropleth")
Encoding(levels(IreMap$NAME_1)) <- "latin1"


##################################################################
## Map
##################################################################

library(colorspace)  

pdf(file="figs/population.pdf")
quantPal <- rev(heat_hcl(16))
spplot(IreMap["POP"], col='transparent', col.regions=quantPal)
dev.off()

pdf(file="figs/cases.pdf")
quantPal <- rev(heat_hcl(16))
spplot(IreMap["CASES"], col='transparent', col.regions=quantPal)
dev.off()

pdf(file="figs/cases per pop.pdf")
quantPal <- rev(heat_hcl(16))
spplot(IreMap["Cases.Pop"], col='transparent', col.regions=quantPal)
dev.off()

pdf(file="figs/case fatality.pdf")
quantPal <- rev(heat_hcl(16))
spplot(IreMap["Case.fatal"], col='transparent', col.regions=quantPal)
dev.off()
