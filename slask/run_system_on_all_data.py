from glob import glob
import numpy as np
import scipy as sp
import json
import os
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import sys

if os.path.split(os.getcwd())[-1] != 'structure-from-sound-python':
    os.chdir("../")

sys.path.append('.')
from src.detectors import gcc_phat_detector
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
import src.tdoa_datasets_module as tdoa_datasets_module
import src.system_function_stuff

datasets, experiments_dict = tdoa_datasets_module.get_data_paths("./data")

config = json.load(open("config.json","r"))

for dataset in datasets:
    for experiment_path in experiments_dict[dataset]:
        print(experiment_path)
        src.system_function_stuff.run_system(experiment_path, config["matlab"] )