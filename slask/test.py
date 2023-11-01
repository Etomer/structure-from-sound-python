import os
import sys
import json
sys.path.append('.')
import src.system_general


config = json.load(open("config.json","r"))

src.system_general.run_system("./data/tdoa_20201016/data/chirp_0001", config["matlab"] )

# print(os.getcwd())
# dataset_name = "tdoa_20201016"
# experiment_name = "music_0015"
# output_folder = "./results/" + dataset_name + "/data/" + experiment_name + "/"
# # tdoa_matrix_to_tdoa_vector(
# #     output_folder, output_folder=output_folder, cutoff_fraction_of_all_measuremnets=1/2)
