# this is a function written for the students to acquire
# a time-series for a specific location from a NetCDF file

get_ts_as_csv <- function(file_name,
                          lon,
                          lat,
                          output_name) {
  library(terra)

  message(paste0("Reading the file: ", file_name))

  file_rast <- rast(file_name)

  file_crds <- crds(file_rast, df = FALSE, na.rm = TRUE)

  file_lat <- sort(unique(file_crds[, 1]), decreasing = T)

  file_lon <- sort(unique(file_crds[, 2]))


  lon_index <- which(abs(file_lon - lon) == min(abs(file_lon - lon)))

  message(paste0("The closest longitude to the one you provided is ", file_lon[lon_index]))


  lat_index <- which(abs(file_lat - lat) == min(abs(file_lat - lat)))

  message(paste0("The closest latitude to the one you provided is ", file_lat[lat_index]))

  file_value <- t(as.matrix(file_rast[lon_index, lat_index, , drop = TRUE]))

  ts_df <- data.frame(
    TIME = terra::time(file_rast),
    VALUE = file_value
  )

  names(ts_df) <- c("TIME", varnames(file_rast))

  message(paste0("Writing csv file: ", output_name))

  write.csv(
    x = ts_df,
    file = output_name,
    row.names = FALSE
  )

  message("DONE ;-)")
  rm(list = ls())
  gc()
}
# test
# file_name <- "/media/ahmed/Volume/CORDEX/tas_AFR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_19700101-19701230_corrected.nc"
# lon <- 0
# lat <- 0
# output_name <- "/media/ahmed/Volume/CORDEX/tas_AFR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_19700101-19701230_corrected.csv"
# 
# get_ts_as_csv(
#   file_name,
#   lon,
#   lat,
#   output_name
# )
