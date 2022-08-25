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

## Create output directories
outputdir <- file.path(outputroot, "analysis", aggr_period, "nh-input") #, "time_series")
dir.create(outputdir, recursive = TRUE)
dir.create(file.path(outputdir, "time_series"))
dir.create(file.path(outputdir, "attributes"))

## ####################################################### ##
## 1 - Time series data
## ####################################################### ##

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
  timedim <- ncdim_def("date", tunits, as.double(nc_dates))

  ## Define variables
  fillvalue <- 1e32
  dlname <- "Q95"
  q95_def <- ncvar_def("Q95", "m3 d-1", timedim, fillvalue, dlname)
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
  nc_filename <- file.path(outputdir, "time_series", paste0(id, ".nc"))
  ncout <- nc_create(nc_filename, list(q95_def, nao_def, ea_def, amv_def, prec_def, tmp_def), force_v4 = TRUE)
  ncvar_put(ncout, q95_def, xx$Q_95)
  ncvar_put(ncout, nao_def, xx$nao)
  ncvar_put(ncout, ea_def, xx$ea)
  ncvar_put(ncout, amv_def, xx$amv)
  ncvar_put(ncout, prec_def, xx$european_precip)
  ncvar_put(ncout, tmp_def, xx$uk_temp)
  nc_close(ncout)
}

## ####################################################### ##
## 2 - Static attributes
## ####################################################### ##

camels_attr_fs <- c("CAMELS_GB_climatic_attributes.csv",
                    "CAMELS_GB_humaninfluence_attributes.csv",
                    "CAMELS_GB_hydrogeology_attributes.csv",
                    "CAMELS_GB_hydrologic_attributes.csv",
                    ## "CAMELS_GB_hydrometry_attributes.csv",
                    "CAMELS_GB_landcover_attributes.csv",
                    "CAMELS_GB_soil_attributes.csv",
                    "CAMELS_GB_topographic_attributes.csv")

camels_datadir <- file.path(config$aux_data$camels, "data")
destdir <- file.path(outputdir, "attributes")
dir.create(destdir, recursive = TRUE)
for (i in 1:length(camels_attr_fs)) {
  filepath <- file.path(camels_datadir, camels_attr_fs[i])
  file.copy(filepath, destdir, overwrite = TRUE)
}


## ####################################################### ##
## 3 - Basin list
## ####################################################### ##

conn <- file(file.path(outputdir, "basins.txt"))
writeLines(as.character(station_ids[1:5]), conn)
close(conn)

## ## ####################################################### ##
## ## 4 - YAML configuration file
## ## ####################################################### ##

## template <- read_yaml("resources/nh-config-template.yml")
## template$experiment_name <- paste0("cudalstm_", n_stations, "_basins_", aggr_period)
## template$train_basin_file <- "basins.txt"
## template$validation_basin_file <- "basins.txt"
## template$test_basin_file <- "basins.txt"
## template$train_start_date <- "01/12/1961"
## template$train_end_date <- "01/12/2006"
## template$validation_start_date <- "01/12/1961"
## template$validation_end_date <- "01/12/2006"
## template$test_start_date <- "01/12/1961"
## template$test_end_date <- "01/12/2006"
## ## template$device <- "cuda:0"
## template$device <- "cpu"
## template$validate_every <- as.integer(3)
## template$validate_n_random_basins <- as.integer(1)
## template$metrics <- "NSE"
## template$model <- "cudalstm"
## template$head <- "regression"
## template$output_activation <- "linear"
## template$hidden_size <- as.integer(20)
## template$initial_forget_bias <- as.integer(3)
## template$output_dropout <- 0.4
## template$optimizer <- "Adam"
## template$loss <- "MSE"
## template$learning_rate <- "-0: 1e-2\n-30: 5e-3\n-40: 1e-3" #list("0" = "1e-2", "30"= "5e-3", "40"= "1e-3")
## ## template$learning_rate <- list("0" = "1e-2", "30" = "5e-3", "40" = "1e-3")
## template$batch_size <- as.integer(256)
## template$epochs <- as.integer(50)
## template$clip_gradient_norm <- as.integer(1)
## template$predict_last_n <- as.integer(1)
## template$seq_length <- as.integer(8)
## template$num_workers <- as.integer(8)
## template$log_interval <- as.integer(5)
## template$log_tensorboard <- as.integer(0) # FIXME should be True or False, but can't remove quotes [0/1 should be correctly interpreted by Python yaml loader]
## template$log_n_figures <- as.integer(0)
## template$save_weights_every <- as.integer(1)
## template$dataset <- "generic"
## template$data_dir <- noquote(".")
## template$dynamic_inputs <- c("NAO", "EA", "AMV", "P", "T")
## template$target_variables <- "Q95"
## template$clip_targets_to_zero <- "Q95"

## write_yaml(template, file.path(outputdir, "basins.yml"), fileEncoding = "UTF-8")

## test <- read_yaml(file.path(outputdir, "basins.yml"))
## current <- read_yaml("basin.yml")

## items <- names(test)
## for (i in 1:length(items)) {
##   nm <- items[i]
##   x <- test[[nm]]
##   y <- current[[nm]]
##   if (!isTRUE(all.equal(x, y))) {
##     cat(sprintf("Mismatch for item: %s\n", nm))
##     cat(sprintf("New: %s\n", x))
##     cat(sprintf("Ref: %s\n", y))
##   }
## }
