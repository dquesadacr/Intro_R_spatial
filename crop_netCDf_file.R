
crop_netCDf_file <- function(input_file,
                             shape_file,
                             output_file) {
  require(terra)

  r <- rast(input_file)

  v <- vect(shape_file)

  r_out <- crop(x = r, y = v, mask = TRUE, snap = "in")

  writeCDF(r_out,
    filename = output_file,
    varname = varnames(r),
    longname = longnames(r),
    compression = 9,
    overwrite = TRUE
  )

  message("DONE ;-)")
  rm(list = ls())
  gc()
}

# test
input_file <- "/media/ahmed/Volume/CORDEX/tas_AFR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_19700101-19701230_corrected.nc"

shape_file <- "sdn_adm_cbs_nic_ssa_20200831_shp/sdn_admbnda_adm0_cbs_nic_ssa_20200831.shp"

output_file <- "/media/ahmed/Volume/CORDEX/tas_AFR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_19700101-19701230_cropped.nc"
