import numpy as np
import os
import pandas as pd
from pathlib import Path

from src.detectors import gcc_phat_detector
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector

def run_system(experiment_path, matlab_path):
    # Paths -----------------------
    # Given relative to base of direcotry

    input_folder = experiment_path
    path_parts = []
    (head, tail) = os.path.split(experiment_path)
    path_parts.append(tail)
    while tail != "structure-from-sound-python":
        (head, tail) = os.path.split(head)
        path_parts.append(tail)

    output_folder = os.path.join(head, tail, "results", *path_parts[1:])
    print(output_folder)

    # create output folder if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    gcc_phat_detector(input_folder, output_folder=output_folder,
                  window_length=10000, speed_of_movement=1)

    tdoa_matrix_to_tdoa_vector(
    output_folder, output_folder=output_folder, cutoff_fraction_of_all_measuremnets=1/2)
    tdoav = np.load(output_folder + "tdoa_vectors.npy")
    df = pd.DataFrame(tdoav)
    df.to_csv(output_folder + "tdoa_vectors_to_matlab.csv")

    matlab = matlab_path

    os.system(matlab + "\" addpath('./matlab/matlab');tdoa_for_python('" +
            output_folder + "');exit\"")

