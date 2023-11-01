import os
import numpy as np
import pandas as pd
import json

config = json.load(open("config.json","r"))

def matlab_tdoa_vector_to_positions(input_folder, output_folder=None):
    if output_folder == None:
        output_folder = input_folder

    # resaves the tdoa_vectors to csv so matlab can read them
    tdoav = np.load(os.path.join(input_folder,"tdoa_vectors.npy"))
    df = pd.DataFrame(tdoav)
    df.to_csv(os.path.join(output_folder ,"tdoa_vectors_to_matlab.csv"))

    # calls the matlab code
    os.system(config["matlab"] + " -nodesktop -nosplash -nodisplay -r " + "\" addpath('./matlab/matlab');tdoa_for_python('" +
            output_folder + "');exit\"")
