import numpy as np
import os
import pandas as pd
from pathlib import Path

from src.detectors import gcc_phat_detector
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
import system_settings

def run_system(experiment_path, matlab_path):
    # Paths -----------------------

    # construct paths
    input_folder = experiment_path
    path_list = experiment_path.split(os.sep)
    path_list[-4] = os.path.join("results", system_settings.system_name)
    output_folder = os.sep.join(path_list)

    # create output folder if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # run detector
    system_settings.detector(input_folder,output_folder=output_folder)

    # run detections_to_tdoa_vector
    system_settings.tdoa_matrix_to_tdoa_vector_function(output_folder)

    # run tdoa_vector_to_positions
    system_settings.tdoa_vector_to_position_function(output_folder)
