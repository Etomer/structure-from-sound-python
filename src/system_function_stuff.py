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

    path_list = experiment_path.split(os.sep)
    path_list[-4] = "results"
    output_folder = os.sep.join(path_list)

    # create output folder if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    gcc_phat_detector(input_folder, output_folder=output_folder,
                  window_length=10000, speed_of_movement=1)

    tdoa_matrix_to_tdoa_vector(
    output_folder, output_folder=output_folder, cutoff_fraction_of_all_measuremnets=1/2)
    tdoav = np.load(os.path.join(output_folder,"tdoa_vectors.npy"))
    df = pd.DataFrame(tdoav)
    df.to_csv(os.path.join(output_folder ,"tdoa_vectors_to_matlab.csv"))

    matlab = matlab_path

    os.system(matlab + " -nodesktop -nosplash -nodisplay -r " + "\" addpath('./matlab/matlab');tdoa_for_python('" +
            output_folder + "');exit\"")

