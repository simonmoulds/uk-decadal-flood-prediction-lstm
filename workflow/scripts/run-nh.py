#!/usr/bin/env python3

import os
import sys
# import shutil
# import datetime
from pathlib import Path

# import matplotlib.pyplot as plt
# import torch
from neuralhydrology.evaluation import metrics
from neuralhydrology.nh_run import start_run, eval_run

# Get command line arguments
config = sys.argv[1]
nh_config = sys.argv[2]
# aggr_period = sys.argv[2]
# outputroot = sys.argv[3]

start_run(config_file = Path(nh_config))

# if torch.cuda.is_available():
#     start_run(config_file = Path('531_basins_arc.yml'))
# else:
#     raise OSError("No GPU available!")
#     # start_run(config_file = Path("531_basins_arc.yml"), gpu = -1)

# Neural Hydrology assigns a unique name to the output run
# directory. We find this name by identifying the most recent
# directory, then evaluate the model output.
run_dirs = [d for d in filter(os.path.isdir('runs'), os.listdir('.'))]
run_dirs.sort(key = lambda x: os.path.getmtime(x))
run_dir = run_dirs[-1]
eval_run(run_dir, period = 'test')
