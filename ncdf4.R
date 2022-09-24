# NetCDF (network Common Data Form) is a file format for 
# storing multidimensional scientific data (variables) such as temperature,
# humidity, pressure, wind speed, and direction. 

library(ncdf4)
library(ncdf4.helpers)

# file name 
setwd("/media/ahmed/Volume/MA-Betreuung/R-Kurse/RKurs-Github/Intro_R_spatial/")
f<-"spatial/example.nc"

ncin<-nc_open(f)

#https://pjbartlein.github.io/REarthSysSci/netCDF.html

print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
head(lat)

# get time
time <- ncvar_get(ncin,"time")
head(time)

time_okay<-ncdf4.helpers::nc.get.time.series(f = ncin)

# get temperature
tmp_array <- ncvar_get(ncin,"tas")

 ncatt_get(ncin,"tas","long_name")
 ncatt_get(ncin,"tas","units")
 ncatt_get(ncin,"tas","_FillValue")
dim(tmp_array)


tmp_slice <- tmp_array[,,5]

library(lattice)
library(RColorBrewer)
# quick map
image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))
