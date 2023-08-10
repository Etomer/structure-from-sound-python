import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
import functools
import sys

sys.path.append('../src')
import tdoa_datasets_module as tdoa

# Settings
# sizes of windows to run detections methods on
window_sizes = [2000, 6000, 10000, 30000]

original_gcc_phat = functools.partial(tdoa.gcc_phat, high_freq_cutoff=False)
methods = {"org_gcc_phat": original_gcc_phat, "gcc_phat": tdoa.gcc_phat}

# runs methods on all data with ground truth
experiments, recordings = tdoa.get_data_paths(
    "../data/", with_ground_truth=True)
results_path = "data_for_plotting"


# create folders
if not os.path.exists(results_path):
    os.mkdir(results_path)
for method_name in methods:
    if not os.path.exists(results_path + "/" + method_name):
        os.mkdir(results_path + "/" + method_name)

for window_size in window_sizes:
    for method_name in methods:
        for experiment in experiments:
            for recording_folder in recordings[experiments[0]]:
                tdoa_chunk_estimation, tdoa_chunk_gt = tdoa.evaluate_tdoa_estimator_on_recording(
                    methods[method_name], recording_folder, chunk_length=window_size)
                np.save(results_path + "/" + method_name + "/" + recording_folder.split("/")
                        [-1] + "_" + str(window_size), tdoa_chunk_estimation)
                np.save(results_path + "/" + method_name + "/" + recording_folder.split("/")
                        [-1] + "_" + str(window_size)+"_gt", tdoa_chunk_gt)
