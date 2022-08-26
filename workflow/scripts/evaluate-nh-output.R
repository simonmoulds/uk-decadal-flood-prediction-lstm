#!/usr/bin/env Rscript

## Author : Simon Moulds
## Date   : Jan 2022

library(tidyverse)
library(lubridate)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(scales)
library(arrow)
library(sf)
library(yaml)

options(dplyr.summarise.inform = FALSE)

## Extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config = read_yaml(args[1])
  outputroot = args[2]
  ## args = commandArgs()
  ## m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  ## cwd <- dirname(regmatches(args, m))
}

## Text size for plot labels
axis_title_size_large = 9
axis_title_size = 8
axis_title_size_small = 7
axis_label_size_large = 7
axis_label_size = 6
axis_label_size_small = 5
legend_label_size = 6
legend_title_size = 8
tag_label_size = 8
strip_label_size = 8

mean_square_error_skill_score = function(obs, exp, na.rm = TRUE) {
  if (na.rm) {
    idx <- is.na(obs) | is.na(exp)
    obs <- obs[!idx]
    exp <- exp[!idx]
  }
  ## mse = mean((exp - obs) ^ 2)
  ## mse_ref = mean((mean(obs) - obs) ^ 2)
  ## msss = 1 - (mse / mse_ref)
  ## correlation
  r = cor(obs, exp, method = "pearson")
  ## potential skill [= coefficient of determination]
  ps = r ^ 2
  ## slope reliability
  srel = (r - (sd(exp) / sd(obs))) ^ 2
  ## standardized mean error
  sme = ((mean(exp) - mean(obs)) / sd(obs)) ^ 2
  msss = ps - srel - sme
  list(msss = msss, ps = ps, srel = srel, sme = sme)
}

## For testing:
outputroot <- "results"

nh_inputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "nh-output", "time_series")

gamlss_inputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "gamlss-output", "prediction")
metadata <- read_parquet(file.path(outputroot, "nrfa-metadata.parquet"))

y <- open_dataset(gamlss_inputdir) %>% collect()
station_ids <- y$ID %>% unique()
n_stations <- length(station_ids)

msss_list <- list()
for (i in 1:n_stations) {
  stn <- station_ids[i]
  catchment_area <- metadata %>% filter(id %in% stn) %>% `$`(catchment_area)

  ## Read GAMLSS model prediction
  yy <- y %>% filter(ID %in% stn) %>% rename(gamlss_exp = exp) %>% dplyr::select(year, obs, gamlss_exp)
  yy$obs <- yy$obs * 86400 / catchment_area * 1000 / 1000 / 1000
  yy$gamlss_exp <- yy$gamlss_exp * 86400 / catchment_area * 1000 / 1000 / 1000

  ## Read LSTM prediction
  xx <- read_csv(file.path(nh_inputdir, paste0(stn, ".csv")), show_col_types = FALSE)
  xx <- xx %>%
    mutate(year = lubridate::year(date)) %>%
    rename(obs = Q95_obs, nh_exp = Q95_sim) %>%
    dplyr::select(year, obs, nh_exp, -date, -time_step) %>%
    mutate(year = year - 1) # FIXME
    ## gather(-season_year, key = "key", value = "value")

  ## Compute MSSS
  idx <- is.na(yy$gamlss_exp) | is.na(xx$nh_exp)
  gamlss_msss <- mean_square_error_skill_score(yy$obs[!idx], yy$gamlss_exp[!idx])$msss
  lstm_msss <- mean_square_error_skill_score(xx$obs[!idx], xx$nh_exp[!idx])$msss
  msss_list[[i]] <- tibble(ID = stn, GAMLSS=gamlss_msss, LSTM=lstm_msss)

  ## Make plot
  xx <- xx %>%
    left_join(yy %>% dplyr::select(-obs), by = c("year")) %>%
    gather(-year, key = "key", value = "value") %>%
    filter(!key %in% "obs")

  p <- ggplot() +
    theme_bw() +
    geom_line(
      aes(y=value, x=year, colour = key), data=xx
    ) +
    ## scale_fill_manual(values = cbbPalette) +
    ## scale_color_manual(values = cbbPalette) +
    ## scale_fill_discrete(values = cbbPalette) +
    ## scale_colour_discrete(values = cbbPalette) +
    ## facet_wrap(. ~ ID, ncol = 1) + #, labeller = label_parsed) +
    ylab(expression(Streamflow~(mm~d^{-1}))) +
    xlab("") +
    ## scale_x_continuous(breaks = pretty_breaks()) +
    ## ## scale_y_continuous(expression(Streamflow~(m^{3}~s^{-1})), breaks = pretty_breaks()) +
    ## ## N.B. use alpha to create another legend
    ## ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
    ## ## geom_point(
    geom_line(
      aes(y=obs, x=year), #, alpha="Observed"),
      color = "black",
      data=yy,
      size = 1 #0.2
    ) +
    ## scale_alpha_manual(name=NULL, values=1, breaks="Observed") +
    ## ggtitle(sprintf("ID = %d", stn_id)) +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.title = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = strip_label_size),
          ## legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_label_size),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = axis_label_size_small),
          axis.text.x = element_text(size = axis_label_size_small))
}

msss <- do.call("rbind", msss_list)
