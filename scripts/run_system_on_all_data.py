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
import src.system_general 
datasets, experiments_dict = tdoa_datasets_module.get_data_paths("./data")
config = json.load(open("config.json","r"))



# settings
run_concurrently = True
n_workers = 10 # number of parallell processes to use, will probably depend on your computers RAM

if __name__ == '__main__':
    if run_concurrently:
        from multiprocessing import pool
        experiment_paths = []
        for dataset in datasets:
            for experiment_path in experiments_dict[dataset]:
                experiment_paths.append(experiment_path)

        pool = pool.Pool(n_workers)
        pool.map(src.system_general.run_system, experiment_paths)
        pool.close()

    else:
        for dataset in datasets:
            for experiment_path in experiments_dict[dataset]:
                print(experiment_path)
                src.system_general.run_system(experiment_path)