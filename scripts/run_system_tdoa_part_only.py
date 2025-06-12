import os,sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np
import pandas as pd
import shutil
from detectors import gcc_phat_detector
from tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
import system_settings

def run_system(experiment_path):
    # Paths -----------------------

    # construct paths
    input_folder = experiment_path
    path_list = experiment_path.split(os.sep)
    path_list[-4] = os.path.join("results", system_settings.system_name)
    output_folder = os.sep.join(path_list)

    # run tdoa_vector_to_positions
    system_settings.tdoa_vector_to_position_function(output_folder)



if __name__ == "__main__":
    run_system(sys.argv[1]) 