# replacing raster with terra 

# install order -----------------------------------------------------------
#install.packages("Rtools")
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libproj-dev 

install.packages("devtools")
devtools::install_github("ropenscilabs/rnaturalearth")
devtools::install_github("ropenscilabs/rnaturalearthhires")

pkgs<-utils::available.packages()

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

setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial")
# save plots 
png("images/matrix_raster_terra.png", width = 800,
    height= 800, res = 150)
#plot 
plot(r1, main = "Raster made from a matrix")
# Plot the center of the pixels
points(crds(r1), pch=3, cex=0.5)

dev.off()

# PDF slide 85 ------------------------------------------------------------
# Run these 4 lines in this order to install the "hires" version of "rnaturalearth"
install.packages("Rtools")
install.packages("devtools")
devtools::install_github("ropenscilabs/rnaturalearth")
devtools::install_github("ropenscilabs/rnaturalearthhires")

setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial/spatial/")
de_dem <- rast("deutschland_dgm.asc")
crs(de_dem)<-"+proj=tmerc +lat_0=0 +lon_0=12 +k=1 +x_0=4500000 +y_0=0 +ellps=bessel +units=m +no_defs +type=crs"
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

# PDF slide 87 ------------------------------------------------------------

sqrt(de_dem)
class       : SpatRaster 
dimensions  : 910, 720, 1  (nrow, ncol, nlyr)
resolution  : 1000, 1000  (x, y)
extent      : 4030000, 4750000, 5230000, 6140000  (xmin, xmax, ymin, ymax)
coord. ref. : Germany_Zone_4 (ESRI:31494) 
source      : memory 
name        : deutschland_dgm 
min value   :         0.00000 
max value   :        52.63412 
de_dem + de_dem*4 # Need to have same dimensions
class       : SpatRaster 
dimensions  : 910, 720, 1  (nrow, ncol, nlyr)
resolution  : 1000, 1000  (x, y)
extent      : 4030000, 4750000, 5230000, 6140000  (xmin, xmax, ymin, ymax)
coord. ref. : Germany_Zone_4 (ESRI:31494) 
source      : memory 
name        : deutschland_dgm 
min value   :         -892.30 
max value   :        13851.75 

# PDF slide 88------------------------------------------------------------
par(mfrow=c(1,3))
terra::hist(de_dem, main="Distribution of elevation \n values",
             breaks=40,maxcell=1000000)
terra::boxplot(de_dem, ylab= "Elevation", main = "Boxplot")
terra::plot(de_dem, main = "Basic plot",
            col = RColorBrewer::brewer.pal(7, "BrBG"), 
            range = c(0, 2500),
            ylim = ext(de_dem)[c(3,4)])



png("../images/histbox_dem_terra.png", width = 900,
    height= 450, res = 150)

par(mfrow=c(1,3))
terra::hist(de_dem, main="Distribution of elevation \n values",
            breaks=40,maxcell=1000000)
terra::boxplot(de_dem, ylab= "Elevation", main = "Boxplot")
terra::plot(de_dem, main = "Basic plot",
            col = RColorBrewer::brewer.pal(7, "BrBG"), 
            range = c(0, 2500),
            ylim = ext(de_dem)[c(3,4)])

dev.off()

# PDF slide 89------------------------------------------------------------

dem_repro <- terra::project(de_dem,
                           "+proj=longlat +datum=WGS84")
dem_repro

class       : SpatRaster 
dimensions  : 732, 901, 1  (nrow, ncol, nlyr)
resolution  : 0.01127346, 0.01128598  (x, y)
extent      : 5.499419, 15.6568, 47.03692, 55.29826  (xmin, xmax, ymin, ymax)
coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
source      : memory 
name        : deutschland_dgm 
min value   :        -138.226 
max value   :        2689.770 

png("../images/reproj_dem_terra.png", width = 800,
    height= 800, res = 150)

terra::plot(dem_repro, col= terrain.colors(12))
dev.off()
# PDF slide 90------------------------------------------------------------

writeRaster(x = dem_repro, "dem_repro_terra.tif")

plot(rast("dem_repro.tif"));plot(rast("dem_repro_terra.tif"))
# PDF slide 91------------------------------------------------------------
terrain_all <- terrain(dem_repro, unit='degrees',
                       v=c("slope", "aspect", "TPI",
                             "TRI", "roughness", "flowdir"))
terrain_all
class(terrain_all)

png("../images/terrain_terra.png", width = 1500,
    height= 1000, res = 150)

plot(terrain_all, cex=1.5)
dev.off()

terrain_all$TRI
# PDF slide 93------------------------------------------------------------
library(rnaturalearth)
bundes <- ne_states(country="germany") # Obtain borders
plot(terrain_all$TRI)
plot(bundes, add=TRUE)
class(bundes) # Notice the class of the object
[1] "SpatialPolygonsDataFrame"
attr(,"package")
[1] "sp"
# SpatRaster can also be created:
c(terrain_all$roughness, terrain_all$TPI)

png("../images/tri_bundes_terra.png", width = 800,
    height= 800, res = 150)
plot(terrain_all$TRI)
plot(bundes, add=TRUE)
dev.off()
# PDF slide 94------------------------------------------------------------
ext(dem_repro)
SpatExtent : 5.4994194913659, 15.6568047584716, 47.0369192127676, 55.2982592300136 (xmin, xmax, ymin, ymax)

par(mfrow=c(2,1))
crop_extent <- ext(c(8,12,50,54))
cropped_dem <- crop(dem_repro, crop_extent)
plot(cropped_dem, main= "Cropped to extent")
plot(bundes, add=TRUE)
masked_dem <- mask(dem_repro, vect(bundes))
plot(masked_dem, main= "Masked to polygon")


png("../images/crop_mask_terra.png", width = 800,
    height= 1500, res = 150)
par(mfrow=c(2,1))
plot(cropped_dem, main= "Cropped to extent", cex=1.5)
plot(bundes, add=TRUE)

plot(masked_dem, main= "Masked to polygon", cex=1.5)
dev.off()


# PDF slide 95------------------------------------------------------------
setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial/spatial/")
library(terra)
kreis_ogr <- vect("kreis.gpkg")
class(kreis_ogr)
[1] "SpatialPolygonsDataFrame"
attr(,"package")
[1] "sp"

png("../images/kreis_ogr_terra.png", width = 800,
    height= 800, res = 150)
plot(kreis_ogr, main = "Default terra plot")
dev.off()

plot(kreis_ogr, main = "Default terra plot")
# PDF slide 96------------------------------------------------------------
library(tidyverse)
kreis_ogrT <- project(kreis_ogr,"EPSG:4326")

plot(dem_repro, xlim = c(11.5,15.5),
     ylim=c(50,52))
plot(kreis_ogrT, add=TRUE)

png("../images/dem_kreis_ogr_terra.png", width = 850,
    height= 600, res = 150)
plot(dem_repro, xlim = c(11.5,15.5),
     ylim=c(50,52))
plot(kreis_ogrT, add=TRUE)
dev.off()

library(sf)
kreis_sf_2 <- st_as_sf(kreis_ogr)
class(kreis_sf_2)
kreis_sf <- read_sf("kreis.gpkg")
kreis_sf == kreis_sf_2

# PDF slide 97------------------------------------------------------------
kreis_ogrSub <- kreis_ogrT[grep("Kreisfreie",  kreis_ogrT$KREIS)]
                           
plot(dem_repro, col= terrain.colors(12),
     xlim = c(11.5,15.5), ylim=c(50,52),
     main = "Main cities in Sachsen from terra")
plot(kreis_ogrSub, add=TRUE)

png("../images/dem_kreis_ogr_sub_terra.png", width = 750,
    height= 500, res = 150)
plot(dem_repro, col= terrain.colors(12),
     xlim = c(11.5,15.5), ylim=c(50,52),
     main = "Main cities in Sachsen from terra", cex = 1.15)
plot(kreis_ogrSub, add=TRUE)
dev.off()



# Slide 98 ----------------------------------------------------------------
library(stringr)
library(dplyr)
kreis_sfSub <- kreis_sfT %>%
  dplyr::filter(str_detect(KREIS, "Kreisfreie"))


kreis_sfT <- st_transform(kreis_sf,
                          "EPSG:4326")

cropped <- crop(dem_repro, kreis_sfT)
masked_dem_sn <- mask(cropped, kreis_sfT)

masked.spdf<- terra::as.data.frame(masked_dem_sn,
                                   xy =TRUE) %>% 
  dplyr::rename(elev = deutschland_dgm)

raster_gg <- ggplot(masked.spdf) +
  geom_tile(aes(fill=elev, x=x, y=y)) +
  geom_sf(data = kreis_sfT, fill=NA,
          colour="black", size = 0.5) +
  geom_sf_label(data = kreis_sfSub, aes(label=KREIS),
                fill=NA, color= "red2", label.size = 0) +
  coord_sf() +
  labs(x=NULL, y=NULL, fill="m.a.s.l.",
       title = "Raster with different vectors") +
  theme_light(base_size = 11) +
  scale_fill_gradientn(colours = terrain.colors(12))

# vect --------------------------------------------------------------------

library(terra)

### SpatVector from a geom matrix
x1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60))
x1<-vect(x1, "polygons")
plot(x1)


x2 <- rbind(c(-180,-20), c(-140,55))
x2<-vect(z[,1:4], "lines")
plot(x2)

x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
hole <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))
z <- rbind(cbind(object=1, part=1, x1, hole=0), cbind(object=2, part=1, x3, hole=0),
           cbind(object=3, part=1, x2, hole=0), cbind(object=3, part=1, hole, hole=1))
colnames(z)[3:4] <- c('x', 'y')

p <- vect(z, "polygons")

plot(p)
# Q -----------------------------------------------------------------------
setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial/spatial/")

library(terra)
library(raster)

# raster 
de_dem1 <- raster::raster("deutschland_dgm.asc")

raster::crs(de_dem1) <- "+proj=tmerc +lat_0=0 +lon_0=12 +k=1 +x_0=4500000 +y_0=0 +ellps=bessel +units=m +no_defs +type=crs"
print(de_dem1)

de_dem_proj_raster<-raster::projectRaster(de_dem1,
                                          crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")



# terra
de_dem2 <- terra::rast("deutschland_dgm.asc")
terra::crs(de_dem2) <- "+proj=tmerc +lat_0=0 +lon_0=12 +k=1 +x_0=4500000 +y_0=0 +ellps=bessel +units=m +no_defs +type=crs"
print(de_dem2)

de_dem_proj_terra<-terra::project(de_dem2,
                                          "+proj=longlat +datum=WGS84 +no_defs +type=crs")

print(de_dem_proj_raster)
print(de_dem_proj_terra)


par(mfrow=c(1,2))

raster::plot(de_dem_proj_raster)
terra::plot(de_dem_proj_terra)
library(rbenchmark)
benchmark("raster" = {raster::projectRaster(de_dem1,
                                           crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")}, 
          "terra" = {terra::project(de_dem2,
                                   "+proj=longlat +datum=WGS84 +no_defs +type=crs")})
