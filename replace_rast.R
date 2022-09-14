# replacing raster with terra 


# PDF slide 84 ------------------------------------------------------------

library(terra)
# Creating a raster from a matrix
r1 <- rast(matrix(rnorm(19*13), nrow = 19), crs =  "EPSG:4326")
# define extent 
ext(r1)<-c(xmin=5, xmax=15, ymin=-5, ymax=10)

r1
class       : SpatRaster 
dimensions  : 19, 13, 1  (nrow, ncol, nlyr)
resolution  : 0.7692308, 0.7894737  (x, y)
extent      : 5, 15, -5, 10  (xmin, xmax, ymin, ymax)
coord. ref. : lon/lat WGS 84 (EPSG:4326) 
source      : memory 
name        :     lyr.1 
min value   : -2.777259 
max value   :  2.850702 

#plot 
plot(r1, main = "Raster made from a matrix")
# Plot the center of the pixels
points(crds(r1), pch=3, cex=0.5)


# PDF slide 85 ------------------------------------------------------------
# Run these 4 lines in this order to install the "hires" version of "rnaturalearth"
#install.packages("Rtools")
#install.packages("devtools")
devtools::install_github("ropenscilabs/rnaturalearth")
devtools::install_github("ropenscilabs/rnaturalearthhires")

setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial/spatial/")
de_dem <- rast("deutschland_dgm.asc")
crs(de_dem) <-  "ESRI:31494"
print(de_dem)
class       : SpatRaster 
dimensions  : 910, 720, 1  (nrow, ncol, nlyr)
resolution  : 1000, 1000  (x, y)
extent      : 4030000, 4750000, 5230000, 6140000  (xmin, xmax, ymin, ymax)
coord. ref. : Germany_Zone_4 (ESRI:31494) 
source      : deutschland_dgm.asc 
name        : deutschland_dgm 

setwd("/home/ahmed")

# PDF slide 86 ------------------------------------------------------------

global(de_dem, 'range', na.rm=TRUE) # min and max
                range     max
deutschland_dgm -178.46 2770.35
global(de_dem, 'mean', na.rm=TRUE)
                  mean
deutschland_dgm 312.5505
# if #1 didnot work use #2
global(de_dem, fun='median', na.rm=TRUE) #1
median(values(de_dem), na.rm = TRUE)#2
[1] 256.21

de_dem <- setMinMax(de_dem) # add range permanently to SpatRaster
print(de_dem)
class       : SpatRaster 
dimensions  : 910, 720, 1  (nrow, ncol, nlyr)
resolution  : 1000, 1000  (x, y)
extent      : 4030000, 4750000, 5230000, 6140000  (xmin, xmax, ymin, ymax)
coord. ref. : Germany_Zone_4 (ESRI:31494) 
source      : deutschland_dgm.asc 
name        : deutschland_dgm 
min value   :         -178.46 
max value   :         2770.35 