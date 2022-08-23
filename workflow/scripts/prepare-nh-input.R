#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(arrow)
library(ncdf4)
library(RNetCDF)
library(rnrfa)
library(RcppRoll)
library(yaml)

options(dplyr.summarise.inform = FALSE)

## Extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config <- read_yaml(args[1])
  aggr_period = args[2]
  outputroot = args[3]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "external/R/utils.R")) # TODO eventually put utils in package
config = parse_config(config)

## First we handle the time series data
outputdir <- file.path(outputroot, "analysis", aggr_period, "nh-input", "time_series")
dir.create(outputdir, recursive = TRUE)

ds <- open_dataset(file.path(outputroot, "analysis", aggr_period, "input")) %>% collect()
station_ids <- ds$ID %>% unique()
n_stations <- length(station_ids)
## subsets <- x$subset %>% unique()

for (i in 1:n_stations) {
  id <- station_ids[i]
  xx <- ds %>% filter(ID %in% id & subset %in% "best_n")
  yrs <- xx$year + xx$lead_time - 1
  dates <- paste(yrs, "12-01", sep = "-") %>% as.POSIXct(tz = "UTC", format = "%Y-%m-%d")
  tunits <- "days since 1970-01-01"
  nc_dates <- RNetCDF::utinvcal.nc(tunits, dates)

  ## Define dimensions
  ## From neuralhydrology docs: "The netcdf has to have one
  ## coordinate called date, containing the datetime index."
  timedim <- ncdim_def("time", tunits, as.double(nc_dates))

  ## Define variables
  fillvalue <- 1e32
  dlname <- "North Atlantic Oscillation index"
  nao_def <- ncvar_def("NAO", "hPa", timedim, fillvalue, dlname)
  dlname <- "East Atlantic index"
  ea_def <- ncvar_def("EA", "hPa", timedim, fillvalue, dlname)
  dlname <- "Atlantic Multidecadal Variability"
  amv_def <- ncvar_def("AMV", "K", timedim, fillvalue, dlname)
  dlname <- "European precipitation"
  prec_def <- ncvar_def("P", "mm d-1", timedim, fillvalue, dlname)
  dlname <- "UK temperature"
  tmp_def <- ncvar_def("T", "K", timedim, fillvalue, dlname)

  ## Create netCDF4 file and put data
  nc_filename <- file.path(outputdir, paste0(id, ".nc"))
  ncout <- nc_create(nc_filename, list(nao_def, ea_def, amv_def, prec_def, tmp_def), force_v4 = TRUE)
  ncvar_put(ncout, nao_def, xx$nao)
  ncvar_put(ncout, ea_def, xx$ea)
  ncvar_put(ncout, amv_def, xx$amv)
  ncvar_put(ncout, prec_def, xx$european_precip)
  ncvar_put(ncout, tmp_def, xx$uk_temp)
  nc_close(ncout)
}

## Now we add the static attributes
camels_attr_fs <- c("CAMELS_GB_climatic_attributes.csv",
                    "CAMELS_GB_humaninfluence_attributes.csv",
                    "CAMELS_GB_hydrogeology_attributes.csv",
                    "CAMELS_GB_hydrologic_attributes.csv",
                    ## "CAMELS_GB_hydrometry_attributes.csv",
                    "CAMELS_GB_landcover_attributes.csv",
                    "CAMELS_GB_soil_attributes.csv",
                    "CAMELS_GB_topographic_attributes.csv")

camels_datadir <- file.path(config$aux_data$camels, "data")
destdir <- file.path(outputroot, "analysis", aggr_period, "nh-input", "attributes")
dir.create(destdir, recursive = TRUE)

for (i in 1:length(camels_attr_fs)) {
  filepath <- file.path(camels_datadir, camels_attr_fs[i])
  file.copy(filepath, destdir, overwrite = TRUE)
}
