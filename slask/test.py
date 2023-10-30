import os
import sys
import json

os.chdir(os.path.dirname(os.path.abspath(__file__)))

os.chdir("..")

sys.path.append('.')
import src.system_function_stuff

config = json.load(open("config.json", "r"))
os.chdir("..")
print("1 ", os.getcwd())

src.system_function_stuff.run_system("/home/terus_work/research-code/structure-from-sound-python/data/tdoa_20201016/data/music_0015", config["matlab"] )

# print(os.getcwd())
# dataset_name = "tdoa_20201016"
# experiment_name = "music_0015"
# output_folder = "./results/" + dataset_name + "/data/" + experiment_name + "/"
# # tdoa_matrix_to_tdoa_vector(
# #     output_folder, output_folder=output_folder, cutoff_fraction_of_all_measuremnets=1/2)
