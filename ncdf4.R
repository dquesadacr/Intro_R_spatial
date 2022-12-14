# NetCDF (network Common Data Form) is a file format for 
# storing multidimensional scientific data (variables) such as temperature,
# humidity, pressure, wind speed, and direction. 

library(ncdf4)
library(ncdf4.helpers)

# file name 
library(ncdf4)
library(ncdf4.helpers)
setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial/spatial/")
f<-"example.nc"

ncin<-nc_open(f)

#https://pjbartlein.github.io/REarthSysSci/netCDF.html

print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
head(lon)
[1]  0  2  4  6  8 10

# get latitude
lat <- ncvar_get(ncin,"lat")
head(lat)
[1]  -89 -87 -85 -83 -81 -79
# get time
time <- ncvar_get(ncin,"time")
head(time)
[1] 46020.5 46021.5 46022.5 46023.5
# convert time 
time_converted<-ncdf4.helpers::nc.get.time.series(f = ncin)
head(time_converted)
[1] "1976-01-01 12:00:00" "1976-01-02 12:00:00" "1976-01-03 12:00:00" "1976-01-04 12:00:00"

# get precipitation 
tmp_array <- ncvar_get(ncin,"pr")
# convert to mm/day 
tmp_array <-tmp_array*86400

ncatt_get(ncin,"pr","long_name")
 
ncatt_get(ncin,"pr","units")
 
ncatt_get(ncin,"pr","_FillValue")

# get dim
dim(tmp_array)
[1] 180  90   4
# take one layer
tmp_slice <- tmp_array[,,1]

library(lattice)
library(RColorBrewer)
# quick map
image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))
# another option 
library(fields)
image.plot(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))

png("../images/ncdf4_plot_fields.png", width = 800,
    height= 800, res = 150)
image.plot(lon,lat,tmp_slice)

dev.off()
# terra -------------------------------------------------------------------
library(terra)
r<-rast(f)

# meta data
r

# convert units
r<-r*86400
r1<-r

ext(r1)<-c(-180,180,-90,90)
png("../images/netCDf_terra.png", width = 800,
    height= 800, res = 150)
plot(r)
dev.off()

plot(r[[1]])


# convert to df
r.df<-terra::as.data.frame(r, xy=T)

# change order of lon
r.df$x<-ifelse(r.df$x>180, r.df$x -360, r.df$x)

library(ggplot2)
library(rnaturalearth)
library(sf)

# to stay PC and avoid boarders issue
spdf_world <- ne_coastline(returnclass = "sf")


ggplot()+
  geom_tile(r.df,
              mapping =  aes(x,y,
                             fill=pr_1))+
  geom_sf(spdf_world, 
          mapping = aes(), 
          fill=NA, color ="black")+
  
  scale_fill_viridis_c(begin = 0,
                       end = 1)+
  theme_light(base_size = 8) +
  scale_x_continuous(limits = c(-180,180), 
                     expand = c(0,0),
                     breaks = seq(-180,180,10))



#
library(robustbase) 
library(readr)

currentDataset <- read_csv("https://statsnotebook.io/blog/data_management/example_data/stars.csv")
res <- lmrob(light ~ temperature,
             data=currentDataset)
summary(res)
cbind(coef(res),confint(res, level = 0.95))
