import os
import matplotlib.pyplot as plt
import sys
import json

if os.path.split(os.getcwd())[-1] != 'structure-from-sound-python':
    os.chdir("../")

sys.path.append('.')

from src import tdoa_datasets_module

import system_settings

if __name__=="__main__":
    approx_room_size = 10  # meters
    datasets, experiments_dict = tdoa_datasets_module.get_data_paths(os.path.join(".", "results", system_settings.system_name))

    experiment_names = []
    tdoa_matrix_inlier_fractions = []
    tdoa_vector_inlier_fractions = []
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
            info = json.load(open(os.path.join(data_folder, "info.json"), 'r'))
            if not info["has_ground_truth"]:
                continue

            try:
                missing_ratio_tdoa_matrix, inlier_ratio_tdoa_matrix, std_tdoa_inliers = tdoa_datasets_module.compute_stats_tdoa_matrix(data_folder, result_folder)
                residuals_tdoa_vector, inlier_ratio_tdoa_vector, std_tdoa_vector_inliers = tdoa_datasets_module.compute_stats_tdoa_vector(data_folder, result_folder)
                distances_receivers, distances_senders, rms_receivers, rms_senders, sender_inliers, sender_inlier_fraction = tdoa_datasets_module.compute_stats_positions(data_folder, result_folder, 0.2)
                
                experiment_names.append(dataset_name+"\\" + experiment_name)
                tdoa_matrix_inlier_fractions.append(inlier_ratio_tdoa_matrix)
                tdoa_matrix_missing_fractions.append(missing_ratio_tdoa_matrix)
                tdoa_vector_inlier_fractions.append(inlier_ratio_tdoa_vector)
                positions_receiver_errors.append(rms_receivers)
                positions_sender_errors.append(rms_senders)
                sender_inlier_fractions.append(sender_inlier_fraction)
            except Exception as e:
                print("Failed for:", dataset_name+"\\" + experiment_name)
                print(e)
    

    fig, ax = plt.subplots(2, 1)
    plt.axes(ax[1])
    plt.scatter(experiment_names, tdoa_matrix_inlier_fractions, label="Tdoa m", alpha=0.8)
    plt.scatter(experiment_names, tdoa_vector_inlier_fractions, label="Tdoa v", alpha=0.8)
    plt.scatter(experiment_names, sender_inlier_fractions, label="Senders", alpha=0.8)
    plt.legend()

    plt.xticks(rotation='vertical')
    plt.ylim(0, 1)
    plt.title("Inlier fraction")
    plt.xlabel("Experiment")
    plt.ylabel("fraction")
    plt.subplots_adjust(bottom=0.4)

    plt.axes(ax[0])
    plt.scatter(experiment_names, positions_sender_errors, label="Senders", alpha=0.8)
    plt.scatter(experiment_names, positions_receiver_errors, label="Receivers", alpha=0.8)
    plt.legend()

    plt.tick_params(
        labelbottom=False
    )
    plt.ylim(0, 1)
    plt.title("Position error")
    plt.ylabel("Rms error (m)")
    plt.subplots_adjust(bottom=0.1)


    plt.show()





