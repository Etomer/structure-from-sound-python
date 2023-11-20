import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil

from src.detectors import gcc_phat_detector
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
import system_settings

def run_system(experiment_path):
    # Paths -----------------------

    # construct paths
    input_folder = experiment_path
    path_list = experiment_path.split(os.sep)
    path_list[-4] = os.path.join("results", system_settings.system_name)
    output_folder = os.sep.join(path_list)

    # create output folder if it doesn't exist
    head = os.path.abspath(output_folder)
    correct_folder = False
    while len(head) > 0:
        head, tail = os.path.split(head)
        if tail == "structure-from-sound-python":
            correct_folder = True
            break
    if not correct_folder:
        raise Exception("Path of results folder is outside of this project, OR the name of this project has been named to something other than structure-from-sound-python")

    shutil.rmtree(output_folder)
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # run detector
    system_settings.detector(input_folder,output_folder=output_folder)

    # run detections_to_tdoa_vector
    system_settings.tdoa_matrix_to_tdoa_vector_function(output_folder)

    # run tdoa_vector_to_positions
    system_settings.tdoa_vector_to_position_function(output_folder)
