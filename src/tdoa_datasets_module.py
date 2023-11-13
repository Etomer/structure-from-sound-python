import scipy as sp
from glob import glob
import json
import numpy as np
import pandas as pd
import time
import os

# default constant values
CHUNK_LENGTH = 5000 # samples
SPEED_OF_SOUND = 343 # m/s
 
#def plot_gcc_phat_image(sounds, fs, tdoa_gt, times_gt, chunk_length=5000)

def interpolate_gt_at_times(recording_folder, interpolation_times):
    position_gt, _, times_gt = read_tdoa_sound_ground_truth(recording_folder)
    speakerPositions = position_gt["speaker"]
    interpolatedGt = [np.interp(interpolation_times, times_gt, speakerPositions[i]) for i in range(3)]
    return np.array(interpolatedGt)



def evaluate_tdoa_estimator_on_recording(tdoa_estimator, recording_folder, print_status=True, chunk_length=CHUNK_LENGTH):
    if print_status:
        print(recording_folder)
        prev_status_update = time.time()
    fs, sounds = read_tdoa_sound(recording_folder)
    position_gt, tdoa_gt, times_gt = read_tdoa_sound_ground_truth(
        recording_folder)

    problems = divide_into_tdoa_chunks(sounds, fs, tdoa_gt, times_gt, chunk_length=chunk_length)

    nmics = problems[0][0].shape[0]
    nproblems = len(problems)

    tdoa_chunk_gt = np.zeros((nmics, nmics, nproblems))
    tdoa_chunk_estimation = np.zeros((nmics, nmics, nproblems))
    for i, problem in enumerate(problems):
        if print_status:
            nowtime = time.time()
            if prev_status_update + 1 < nowtime:
                print(str(i) + "/" + str(len(problems)))
                prev_status_update = nowtime
        tdoa_chunk_estimation[:, :, i] = tdoa_estimator(problem)
        tdoa_chunk_gt[:, :, i] = problem[1]
    return tdoa_chunk_estimation, tdoa_chunk_gt

def gcc_phat(problem, keep_dim = False, high_freq_cutoff=300):
    c =  sp.fft.fft(problem[0])
    nmics = problem[0].shape[0]

    if keep_dim:
        result = np.empty((nmics,nmics,problem[0].shape[1]))
    else:
        result = np.empty((nmics,nmics))
    for mic1 in range(problem[0].shape[0]):
        for mic2 in range(mic1+1, problem[0].shape[0]):
            temp = c[mic1,:]*np.conj(c[mic2,:])
            temp /= np.abs(temp) + 1e-10
            temp[high_freq_cutoff:-high_freq_cutoff] = 0
            res = sp.fft.ifft(temp)

            if keep_dim:
                result[mic1,mic2] = np.real(res)
            else: 
                bi = np.argmax(res)
                result[mic1,mic2] = bi if bi < problem[0].shape[1]/2 else bi - problem[0].shape[1]
            result[mic2,mic1] = -result[mic1,mic2]

    return result


def divide_into_tdoa_chunks(sounds, fs , tdoa_gt, times_gt, chunk_length=CHUNK_LENGTH, keep_missing_gt_chunks=False, speed_of_sound = SPEED_OF_SOUND):
    """
    Divides sound into chunks of a given length, and matches with ground truth for tdoa

    Input : 
        sounds, numpy_array of size=(nmics, sound_length)
        fs, int; sampling frequency
        tdoa_gt, numpy_array of size=(nmics, nmics, gt_len); containing time-difference groudn truth between pair of mics (in the unit samples)
        times_gt, numpy_array of size=(gt_len); contains times of ground_truth sampling, t=0 is when sound recording started

    Output : 
        list((
            sounds, numpy_array of size=(nmics, chunk_length)
            tdoa_gt, numpy_array of size=(nmics, nmics); 
                i.e. this assumes the chunk_length is short enough so that tdoa can be considered constant in the chunk
            time, int
                time since start of the original recording
            ))
    """

    problems = []
    starts = np.arange(0,sounds.shape[-1],chunk_length)[:-1]

    for start in starts:

        sound = sounds[:,start:(start+chunk_length)]

        gt_sample_in_period = (times_gt < (start + chunk_length)/fs) & (times_gt > start/fs)
        if (gt_sample_in_period.sum() == 0) and not keep_missing_gt_chunks:
            continue 

        gt_index = np.where(gt_sample_in_period)[0][0]

        problems.append((sound,tdoa_gt[:,:,gt_index]*fs/speed_of_sound, times_gt[gt_index]))
    return problems



def get_data_paths(data_folder, with_ground_truth=False):
    """
    returns paths to (experiments, recordings), where experiments are collection of recordings

    experiments : list(str) of experiment paths
    recordings : dict(str -> list(str)) with experiment paths as keys and list of recording paths as values
    """
    experiments = glob(data_folder + "/*/")
    
    #re_experiments = [exp.split("/")[-1] for exp in experiments]
    recordings = {exp:[] for exp in experiments}
    for i,exp in enumerate(experiments):
        temp = glob(os.path.join(exp[:-1] + "/data/") + "*")
        for j in temp:
            recordings[exp].append(j)

    if with_ground_truth:
        gt_recordings = {}
        for exp in experiments:
            temp = []
            for recording in recordings[exp]:
                info = json.load(open(recording + "/info.json",'r'))
                if info["has_ground_truth"]:
                    temp.append(recording)
            gt_recordings[exp] = temp
    
        recordings = gt_recordings


    return (experiments, recordings)


def read_tdoa_sound(case_folder):
    """
    Reads and packages a recording into a tuple containing sample rate and a numpy tensor with the shape (number_of_recordings, length_of_recordings)    
    """
    #sound_files = glob(case_folder + "/*.wav")
    
    info = json.load(open(case_folder + "/info.json",'r'))
    n = info["number_of_mics"]

    sounds = []
    for i in range(n):
        fs,s = sp.io.wavfile.read(case_folder + "/Track " + str(i+1)+".wav")
        sounds.append(s)

    return_obj = np.zeros((n,len(sounds[0])))
    for i in range(n):
        return_obj[i] = sounds[i]

    return (fs,return_obj)
    
    
def read_tdoa_sound_ground_truth(case_folder):
    """
    Reads and packages the ground truth to a recording
    
    RETURNS
    tuple : (positions, tdoa, time),
        positions : tuple {"mics":size=(number_of_mics, 3, timesteps), "speaker": size=(3, timesteps)}
        
        tdoa : size=(number_of_mics, numer_of_mics, timesteps),
            with tdoa(i,j) = toa_j - toa_i and tdoa(i,i) = nan, (in the units meters)
            
        time : (timesteps),
            time since recording started, for each of the timesteps
    """
    df_tdoa = pd.read_csv(os.path.join(case_folder, "gt_tdoa.csv"))
    df_pos  = pd.read_csv(os.path.join(case_folder, "gt_positions.csv"))
    info = json.load(open(os.path.join(case_folder, "info.json"),'r'))
    n = info["number_of_mics"]
    t = df_tdoa.shape[0]

    positions = {"mics": np.zeros((n,3,t)), "speaker": np.zeros((3,t))}
    
    tdoa = np.empty((n,n,t))
    time = np.array(df_tdoa['time'])

    dims = ['x','y','z']
    for i in range(n):
        for j,d in enumerate(dims):
            positions["mics"][i,j,:] = df_pos['mic' + str(i+1) + "_" + d]
    for j,d in enumerate(dims):
        positions["speaker"][j,:] = df_pos['speaker' + "_" + d]

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            else:
                tdoa[i,j,:] = df_tdoa['mic' + str(i+1) + "-mic" + str(j+1)]
                tdoa[j,i,:] = -df_tdoa['mic' + str(i+1) + "-mic" + str(j+1)]
    return (positions, tdoa, time)

def procrustes(points_to_map_from,points_to_map_to, tol=0.5, n_iters=100):
    """
    finds robust euclidian transform on the form:
    points_to_map_to = points_to_map_from @ R + t
    using RANSAC (with selecting 3 points)

    Based on : Horn, Berthold K. P. (1987-04-01). "Closed-form solution of absolute orientation using unit quaternions"
    
    Input:
    points_to_map_from : np.array (n x 3),
    points_to_map_to : np.array (n x 3),
    
    Output:
    transform : tuple (R,t)
    """
    most_inliers = 0

    for i in range(n_iters):
        ind = np.random.permutation(points_to_map_from.shape[0])[:3]
        p1 = points_to_map_from[ind]
        p2 = points_to_map_to[ind]

        p1n = p1 - np.mean(p1,axis=0)
        p2n = p2 - np.mean(p2,axis=0)
        u,s,v = np.linalg.svd(p1n.T@p2n)
        R = u@v

        t = np.mean(p2,axis=0) - np.mean(p1@R,axis=0)

        res = np.linalg.norm(points_to_map_from@R + t - points_to_map_to, axis=1)
        inliers = np.sum(res < tol)
        
        if most_inliers < inliers:
            most_inliers = inliers
            best_set = res < tol
    # re-estimate with inlier set
    ind = best_set
    p1 = points_to_map_from[ind]
    p2 = points_to_map_to[ind]

    p1n = p1 - np.mean(p1,axis=0)
    p2n = p2 - np.mean(p2,axis=0)
    u,s,v = np.linalg.svd(p1n.T@p2n)
    R = u@v

    t = np.mean(p2,axis=0) - np.mean(p1@R,axis=0) 
    return (R,t)


def get_detections_with_gt(data_folder, result_folder):
    """
    Returns the detections along with time-synced ground truth 

    Input:
        data_folder : string, Path to folder where ground truth is stored as file "gt_positions.csv"
        result_folder : string, Path to where detections are stored as "detections.npy" and detection_times as "detections_times.npy"

    Output:
        detections : np.array(n_windows, n_mics, n_mics), detections

        tdoa_gt : np.array(n_windows, n_mics, n_mics), ground truth

        detection_times : np.array(n_windows), times for detections 
    """

    detections = np.load(os.path.join(result_folder, "detections.npy"))
    detection_times = np.load(os.path.join(result_folder, "detection_times.npy"))


    speaker_gt_pos = interpolate_gt_at_times(
        data_folder, detection_times).T
    positions, tdoa, time = read_tdoa_sound_ground_truth(
        data_folder)
    receiver_gt_positions = np.nanmedian(positions["mics"], axis=2)
    nmics = receiver_gt_positions.shape[0]

    tdoa_gt = np.zeros((speaker_gt_pos.shape[0], nmics, nmics))

    for i in range(nmics):
        for j in range(i+1, nmics):
            # check sign here
            tdoa_gt[:, i, j] = (np.linalg.norm(speaker_gt_pos - receiver_gt_positions[j, :],
                                            axis=1) - np.linalg.norm(speaker_gt_pos - receiver_gt_positions[i, :], axis=1))
            tdoa_gt[:, j, i] = -tdoa_gt[:, i, j]
    return detections, tdoa_gt, detection_times

def get_positions_with_gt(data_folder, result_folder):
    """
    Returns the detections along with time-synced ground truth 

    Input:
        data_folder : string, Path to folder where ground truth is stored as file "gt_positions.csv"
        result_folder : string, Path to where positions are stored as "sender_positions.csv","receiver_positions.csv" and times as "tdoa_vector_times.npy"

    Output:
        receiver_positions 
        sender_positions 
        receiver_gt_positions
        sender_gt_positions
    """
    tdoa_vector_times = np.load(os.path.join(
        result_folder, "tdoa_vector_times.npy"))

    sender_gt_positions = interpolate_gt_at_times(
        data_folder, tdoa_vector_times).T
    positions, tdoa, time = read_tdoa_sound_ground_truth(
        data_folder)
    receiver_gt_positions = np.nanmedian(positions["mics"], axis=2)
    sender_positions = pd.read_csv(os.path.join(
        result_folder, "sender_positions.csv"), header=None).to_numpy()
    receiver_positions = pd.read_csv(os.path.join(
        result_folder, "receiver_positions.csv"), header=None).to_numpy()

    R, t = procrustes(
        receiver_positions.T, receiver_gt_positions, tol=0.5, n_iters=1000)
    receiver_positions = receiver_positions.T @ R + t

    sender_positions = pd.read_csv(os.path.join(
        result_folder, "sender_positions.csv"), header=None).to_numpy()
    sender_positions = sender_positions.T @ R + t
    return receiver_positions, sender_positions, receiver_gt_positions, sender_gt_positions

def compute_stats_tdoa_matrix(data_folder, result_folder, tol=0.3):
    detections, tdoa_gt, _ = get_detections_with_gt(data_folder, result_folder)
    res = detections - tdoa_gt
    missing_ratio_tdoa_matrix = np.sum(np.isnan(res))/np.size(res)
    res = res[np.logical_not(np.isnan(res))]
    inlier_ratio_tdoa_matrix = sum(abs(res) < tol)/np.size(res)
    tdoa_inlier_vals = res[abs(res) < tol]
    std_tdoa_inliers = np.std(tdoa_inlier_vals)
    return missing_ratio_tdoa_matrix, inlier_ratio_tdoa_matrix, std_tdoa_inliers

def compute_stats_tdoa_vector(data_folder, result_folder, tol=0.3):
    tdoav = np.load(os.path.join(result_folder,"tdoa_vectors.npy")) # load estimation

    receiver_positions, sender_positions, receiver_gt_positions, sender_gt_positions = get_positions_with_gt(data_folder,result_folder)
    tdoav_gt = sp.spatial.distance.cdist(sender_gt_positions,receiver_gt_positions)
    def sync_tdoa_to_gt(tdoav, tdoav_gt, tol = tol):
        re_tdoav = np.zeros(tdoav.shape)
        for i in range(tdoav.shape[0]):
            t = tdoav[i]
            tgt = tdoav_gt[i]
            
            most_inliers = -1
            for j in range(t.shape[0]):
                offset = tgt[j] - t[j]
                res = t + offset - tgt
                inliers = np.abs(res) < tol
                if np.sum(inliers) > most_inliers:
                    best_inliers = inliers
                    most_inliers = np.sum(inliers)
                    best_offset = offset
            A = np.ones((most_inliers,1))
            b = tgt[best_inliers] - t[best_inliers]    
            best_offset = np.linalg.lstsq(A,b,rcond=None)
            re_tdoav[i] = t + best_offset[0]
        return re_tdoav

    synced_tdoav = sync_tdoa_to_gt(tdoav, tdoav_gt)
    residuals = synced_tdoav - tdoav_gt
    
    residuals_flatten = np.ndarray.flatten(residuals)
    temp = np.ndarray.flatten(residuals < tol)
    inlier_ratio_tdoa_matrix = np.nansum(temp)/temp.shape[0]
    std_tdoa_inliers = np.std(residuals_flatten[temp])

    return residuals, inlier_ratio_tdoa_matrix, std_tdoa_inliers

def compute_stats_positions(data_folder, result_folder, inlier_tol = 0.3):
    receiver_positions, sender_positions, receiver_gt_positions, sender_gt_positions = get_positions_with_gt(
    data_folder, result_folder)
    distances_receivers = np.linalg.norm(receiver_gt_positions-receiver_positions, axis=1)
    distances_senders = np.linalg.norm(sender_gt_positions-sender_positions, axis=1)
    inlier_senders = distances_senders<inlier_tol

    fraction_inlier_senders = sum(inlier_senders)/np.size(inlier_senders)
    rms_receivers = np.sqrt(np.mean(np.square(distances_receivers[~np.isnan(distances_receivers)])))
    rms_senders = np.sqrt(np.mean(np.square(distances_senders[inlier_senders])))

    return distances_receivers, distances_senders, rms_receivers, rms_senders, inlier_senders, fraction_inlier_senders
    
    