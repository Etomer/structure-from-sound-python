from glob import glob
import numpy as np
import scipy as sp
import json
import os
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import sys

if os.path.split(os.getcwd())[-1] != 'structure-from-sound-python':
    os.chdir("../")

sys.path.append('.')

from src.detectors import gcc_phat_detector
from src.tdoa_matrix_to_tdoa_vector import tdoa_matrix_to_tdoa_vector
from src import tdoa_datasets_module

# import importlib
# importlib.reload(tdoa_datasets_module)
import system_settings

def plotTdoaMatrixInliers(experiment_names, tdoa_matrix_inlier_fractions):
    plt.figure()
    plt.scatter(experiment_names, tdoa_matrix_inlier_fractions)
    plt.xticks(rotation='vertical')
    # plt.xticks(rotation=45)
    plt.ylim(0, 1)
    plt.title('Fraction of inliers in tdoa matrix')
    plt.xlabel("Experiment")
    plt.ylabel("Fraction")
    # plt.show()

def plotSenderErrors(experiment_names, positions_sender_errors):
    plt.figure()
    plt.scatter(experiment_names, positions_sender_errors)
    plt.xticks(rotation='vertical')
    plt.ylim(0, 1)
    plt.title('RMS error for senders')
    plt.xlabel("Experiment")
    plt.ylabel("Rms error (m)")

if __name__=="__main__":
    approx_room_size = 10  # meters
    datasets, experiments_dict = tdoa_datasets_module.get_data_paths(os.path.join(".", "results", system_settings.system_name))
    # datasets = glob(os.path.join(".", "results", system_settings.system_name) + "/*/")
    
    #re_experiments = [exp.split("/")[-1] for exp in experiments]
    # experiments_dict = {exp:[] for exp in datasets}
    # for i,exp in enumerate(datasets):
    #     temp = glob(os.path.join(exp[:-1] + "/data/") + "*")
    #     for j in temp:
    #         experiments_dict[exp].append(j)

    experiment_names = []
    tdoa_matrix_inlier_fractions = []
    tdoa_matrix_missing_fractions = []
    positions_receiver_errors = []
    positions_sender_errors = []
    sender_inlier_fractions = []
    
    for dataset_result_path in datasets:
        for experiment_result_path in experiments_dict[dataset_result_path]:
            head, dataset_name = os.path.split(dataset_result_path)
            if len(dataset_name) == 0:
                _, dataset_name = os.path.split(head)
            _, experiment_name = os.path.split(experiment_result_path)
            data_folder = os.path.join(".", "data", dataset_name, "data", experiment_name)
            result_folder = os.path.join(
                ".", "results", system_settings.system_name, dataset_name, "data", experiment_name)
            
            try:
                missing_ratio_tdoa_matrix, inlier_ratio_tdoa_matrix, std_tdoa_inliers = tdoa_datasets_module.compute_stats_tdoa_matrix(data_folder, result_folder)
                distances_receivers, distances_senders, rms_receivers, rms_senders, sender_inliers, sender_inlier_fraction = tdoa_datasets_module.compute_stats_positions(data_folder, result_folder, 0.2)
                
                experiment_names.append(dataset_name+"\\" + experiment_name)
                tdoa_matrix_inlier_fractions.append(inlier_ratio_tdoa_matrix)
                tdoa_matrix_missing_fractions.append(missing_ratio_tdoa_matrix)
                positions_receiver_errors.append(rms_receivers)
                positions_sender_errors.append(rms_senders)
                sender_inlier_fractions.append(sender_inlier_fraction)
            except FileNotFoundError as e:
                print("Failed for:", dataset_name+"\\" + experiment_name)
                print(e)
            except Exception as e:
                print("Failed for:", dataset_name+"\\" + experiment_name)
                print(e)
    
    plotTdoaMatrixInliers(experiment_names, tdoa_matrix_inlier_fractions)
    plotSenderErrors(experiment_names, positions_sender_errors)

    plt.show()





