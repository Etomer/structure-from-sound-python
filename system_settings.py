# Settings file for system

# imports
from functools import partial
import sys
import os
#if os.path.split(os.getcwd())[-1] != 'structure-from-sound-python':
#    os.chdir("../")
#sys.path.append('.')
import src.detectors 
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
from src.tdoa_vector_to_position import python_tdoa_vector_to_positions

# General settings
system_name = "python_solver" #"gcc_long"#"learned_detector2" # name to store the results under
#system_name = "gcc_long"#"learned_detector2" # name to store the results under
#system_name = "chirp_detector" # name to store the results under

# functions for system to run 
detector = partial(src.detectors.learned_detector, #src.detectors.chirp_detector, , # detector function
                   # detector settings 
                    window_length=10000, 
                    #chirp_gt_path = "data/tdoa_20201016/meta/reference sound/Chirp 2 Hz.wav"
)

tdoa_matrix_to_tdoa_vector_function = partial(tdoa_matrix_to_tdoa_vector, # detections_to_tdoa_vector function
                                        cutoff_fraction_of_all_measuremnets=1/2
)

tdoa_vector_to_position_function = partial(python_tdoa_vector_to_positions, # tdoa_vector_to_positions function
                                           
)


