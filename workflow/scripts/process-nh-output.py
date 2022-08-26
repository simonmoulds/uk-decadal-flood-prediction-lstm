#!/usr/bin/env python3

import os
import sys
import shutil
import datetime
import pickle
from pathlib import Path
from tqdm import tqdm

import pyarrow
import pandas as pd
import xarray as xr
from ruamel.yaml import YAML

# Get command line arguments
config = sys.argv[1]
aggr_period = sys.argv[2]
outputroot = sys.argv[3]

# For testing:
# config = 'config/config.yml'
# outputroot = 'results'
# aggr_period = 'yr2to9_lag'
run_dir = 'cudalstm_499_basins_yr2to9_lag_2608_112005'

inputdir = os.path.join(outputroot, 'analysis', aggr_period, 'nh-input')
outputdir = os.path.join(outputroot, 'analysis', aggr_period, 'nh-output')
try:
    os.makedirs(outputdir)
except FileExistsError:
    pass

try:
    os.makedirs(os.path.join(outputdir, 'time_series'))
except FileExistsError:
    pass

with open(os.path.join(inputdir, 'runs', run_dir, 'test', 'model_epoch075', 'test_results.p'), 'rb') as fp:
    results = pickle.load(fp)

station_ids = list(results.keys())
n_stations = len(station_ids)
freq = '1AS-DEC'

for i in range(n_stations):
    stn = station_ids[i]
    df = results[stn][freq]['xr'].to_dataframe()
    csv_filename = str(stn) + '.csv'
    df.to_csv(os.path.join(outputdir, 'time_series', csv_filename))
